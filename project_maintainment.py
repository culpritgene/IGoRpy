#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 20:50:34 2018

@author: lunar
"""
import os
os.chdir('/home/lunar/Desktop/IgoR_crawling/CODE_Project/')
from load_data import *
import project_results as res_pj
from utils import *


warnings.simplefilter('ignore')

### we should put the contents of navigational dictionary to new class - much simpler to hold and remodel
class project_navigation():
    def __init__(self, Core_path=None, primers=None, Inputs=None, Outputs=None, species=None, ALL=None ):
        if ALL: 
            Core_path = ALL['Core_path']
#            primers = ALL['Primers']
            Inputs = ALL['Inputs']
            Outputs = ALL['Outputs']
            species = ALL['Species']
        else: 
            pass
        
        self.Core_path = Core_path
        if species and not Outputs:
            [sub_run('mkdir {0}Outputs/{1}'.format(self.Core_path, x)) for x in species]
        self.Inputs =  {gr0: gr1+'/' for gr1, gr0 in [re.match('(^.*/(.*$))', x).groups() for x in glob.glob(Core_path+'Data/*')]}
        self.Outputs = {gr0: gr1+'/' for gr1, gr0 in [re.match('(^.*/(.*$))', x).groups() for x in glob.glob(Core_path+'Outputs/*')]}
        self.species = list(self.Inputs.keys())
        self.Primers = primers
        if not primers:
            print("Attention, no base structure for Igor's mandatory files!")
            
        def space_out_new():
            pass
               
        
        
class new_project():
    def __init__(self, species, name, project_navigation, init_model):
        self.__dict__.update(project_navigation.__dict__)
        self.name = name
        self.species = species
        self.init_model = init_model
        self.subsets = {sp:{} for sp in self.species}
        self.batches = {sp:[] for sp in self.species}
        self.base_models = {sp:None for sp in self.species}
        self.derived_models = {sp:{} for sp in self.species}
        self.derived_models_paths = {sp:{} for sp in self.species}

#####
    def store_meta(self):
        print('storing Meta in Core folder')
        json.dump( {'Species':self.species, 'Inputs':self.Inputs, 'Outputs':self.Outputs},
                      open(self.Core_path+self.name+'_meta.json', 'w'))

#####--Make-A-New
    def initialize_for_multiple_samples(self):
        for sp in self.species:
            full_sample_paths = glob.glob(self.Inputs[sp]+'*tsv')
            for f in full_sample_paths:
                batch, subset_full_name = make_subset_from_sample_data(sp, f)
                self.batches[sp].append(batch)
                self.subsets[sp].update( {batch: subset_full_name })
    
    def initialize_base_models(self):
        for sp in self.species:
            new_base_model_name = self.Outputs[sp]+sp+'_base_model.txt'
            with open(new_base_model_name, 'w') as f:
                F2 = gen_model_parms(1, self.Primers)
                f.write(F2)
            self.base_models[sp] = new_base_model_name

#####--Load----------------------------
                
    def retrive_batches(self):
        for sp in self.species:
            Temp = [re.match('.*/(.*)_inference', x) for x in glob.glob(self.Outputs[sp]+'*')]
            [self.batches[sp].append(x.group(1)) for x in Temp if x]
            
    def upload_samples(self):
        for sp in self.species:
            for batch in self.batches[sp]:
                Temp = [re.findall('.*{}.*_SEQ.txt'.format(batch), x) for x in glob.glob(self.Inputs[sp]+'processed/*')]
                self.subsets[sp].update( {batch: [x for x in Temp if x][0][0]} )   
    
    def upload_base_model(self):
        for sp in self.species:
            self.base_models[sp] = glob.glob(self.Outputs[sp]+'*.txt')[0]
    
    def load_compiled_model_with_files(self):
        self.retrive_batches()
        self.upload_base_model()
        self.upload_samples()
        
####-Execute--------------------------------
    def align_ig(self, sp, batch,):
        command_string = ('{0} -set_wd {1} -batch {2} -set_genomic --V {3} -set_genomic --J {4} -set_genomic --D ' 
        '{5} -read_seqs {6} -align --all'.format(self.Primers['Igor'], self.Outputs[sp], 
                        batch, self.Primers['Segments']['Vs'], self.Primers['Segments']['Js'], 
                                self.Primers['Segments']['Ds'], self.subsets[sp][batch], self.base_models[sp]))
        # make a call
        print(command_string)
        print(f'aligning {batch}')
        p = subprocess.run('{}'.format(command_string), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print(f'alligned successfully: {batch}')


    def infer_ig(self, sp, batch, N_iters=5):
        command_string = ('{0} -set_wd {1} -batch {2} -set_custom_model {3} -infer --N_iter {4}'.format(
            self.Primers['Igor'], self.Outputs[sp], batch, self.base_models[sp], N_iters))
        print(command_string)
        print(f'infering {batch} with N_iters {N_iters}')
        p = subprocess.run('{}'.format(command_string), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print(f'infered successfully: {batch}')
    
    def generate_ig(self, sp, batch, Number_of_seqs, name_flag, additional_flags, return_CDR3=True):
        anchors = ''
        CDR3= ''
        if return_CDR3:
            anchors = '-set_CDR3_anchors --V {} --J {}'.format(self.Primers['Anchors']['Vs'], 
                                                                               self.Primers['Anchors']['Js'])
            CDR3 = '--CDR3 '
        command_string = ('{0} -set_wd {1} -batch {2} -set_custom_model {3} {4} {5} -generate {6} {7}--name {8} {9}'.format(
            self.Primers['Igor'], self.Outputs[sp], batch, self.derived_models_paths[sp][batch][0], self.derived_models_paths[sp][batch][1], anchors,
                                 Number_of_seqs, CDR3 ,name_flag, additional_flags))
        print(command_string)
        print(f'generating from {batch} with Number of sequences {Number_of_seqs} and {additional_flags}')
        p = subprocess.run('{}'.format(command_string), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print(f'generated successfully! (from {batch})')
 
    def run_for_all(self):
        ## Create decorator that makes above functions run against whole sample of batches
        pass

####-Load_Results-----------------------------    
       
    
    def upload_derived_models(self):
        for sp in self.species:
            for batch in self.batches[sp]:
                path_to_pars = self.Outputs[sp]+batch+'_inference/final_parms.txt'
                path_to_margs = self.Outputs[sp]+batch+'_inference/final_marginals.txt'

                self.derived_models_paths[sp][batch] = (path_to_pars, path_to_margs)
                self.derived_models[sp][batch] = genmod.GenModel( path_to_pars, path_to_margs)



class Multi():
    """ Class wrapper for basic dict. this dictionary contains all projects uploaded from project names. """
    def __init__(self, core, primers, projects):
        self.core = core
        self.primers = primers
        self.project_names = projects
        self.projects = {pr: {'meta': json.load(open(core+pr+'/'+pr+'_meta.json', 'r')) } for pr in projects}
     
    def load_frames(self):
        for pr in self.project_names:
            print(f'loading frame of project {pr}' )
            self.projects[pr]['frame']= new_project(self.projects[pr]['meta']['Species'], pr, 
                         project_navigation(primers=self.primers, ALL=self.projects[pr]['meta']), 1)
        ### load all files and all models
            self.projects[pr]['frame'].load_compiled_model_with_files()
            self.projects[pr]['frame'].upload_derived_models()
            
    def load_results(self, dist_F=JSD,  ):
         for pr in self.project_names:
            self.projects[pr]['res']= res_pj.project_models_results(self.projects[pr]['frame'])
            self.projects[pr]['res'].get_events_probs()
            self.projects[pr]['res'].store_unstacked_results()
            self.projects[pr]['res'].store_unstacked_results_2(dist_f=dist_F)
            
        
        
        
        
        
        
        
        
        
        




                                                                       