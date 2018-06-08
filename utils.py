#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 20:49:59 2018

@author: lunar
"""

import os
import json

os.chdir('/home/lunar/Desktop/IgoR_crawling/CODE_Project/')
from load_data import *

def make_subset_from_sample_data(species, path_to_file, name=None, frame='Out', Quantile=0.99):
    N = 18000
    ## also use date, just for tracebility
    now = datetime.datetime.now()
    date = '-'.join( map(str, [now.year, now.month, now.day]))
    
    r = re.match('(.*/)(.*).tsv', path_to_file)
    path = r.group(1)
    if not name:
        name = r.group(2)
    # we will use central 80% by length only for investigating type of frame distribution
    sample = pd.read_csv(path_to_file , sep='\t')
    
    # select by frame type (alway use out)
    outs = sample[sample['frame_type']==frame]
    ## strip by length to get central 80%
    outs_central = outs[ (outs['cdr3_length']>outs['cdr3_length'].quantile(q=Quantile))&(outs['cdr3_length']<outs['cdr3_length'].quantile(q=Quantile)) ]
    ### Choose only a fraction of all sequences randomly, here fraction is a set number e.g 18k sequences
    N0_c = outs_central.shape[0]
    try:
        assert N0_c>=N, "there is not enough unique sequences in this file: {} < {}".format(N0_c, N)
    except AssertionError:
        print('Warning: sample from {} has only {} unique sequences from 84% interval'.format(name, N0_c,) )
        if N0_c/N > 0.8:
            print('We will use what we have got, though')
            N=N0_c
        else:
            print('Probably we should discard this sample, or investigate it more thorougly')
#     subset_out = outs_central.iloc[np.random.choice(range(1, outs_central.shape[0]), N)]
    subset_out = outs_central.iloc[random.sample(range(N0_c), N)]
    ## NO NEED TO sort by index, why do you even need to restore order here?
    #subset_out = subset_out.sort_index()
    
    ### make batchname, this batchname will be used to create separate folder with this models data!
    batch = '{}_{}_{}'.format(name, frame ,N)
    ### make subset name, this name will be used to access subset
    subset_full_name = path+'processed/'+batch+'_subset_{}_SEQ.txt'.format(date)
    subset_out.to_csv(subset_full_name[:-8]+'_FILE.txt', index=False)
    subset_out.rearrangement.to_csv(subset_full_name, index=False)
    
    ### add to logs of current project
#     with open(Mouse_logs, 'a') as logs:
#         logs.write('Created new subset ' + subset_full_name +' \n')
#         logs.write('Batch for this subset is ' + batch +' \n')
    
    return batch, subset_full_name

############### MAKE NEW BASE MODELS ##################

def read_template_files(self, templates):
### Read V,J,D genomic templates from template files
        fragments = OrderedDict()
        for f, g in templates.items():
            Genes = {}
            with open(g, 'r') as file:
                while True:
                    l = file.readline()
                    if not l:
                        break
                    m = re.match( r'>(.*)\n', l)
                    if m:
                        Genes.update({m.group(1):file.readline().rstrip('\n')})
                    else:
                        pass
            
            fragments.update({f:Genes})
        return fragments

def gen_model_parms(self, primers):
### Generate initial model parameters from base model with reference templates
        fragments = read_template_files(1, primers['Segments'])
        ## here we just vivisect example model, load new sequences inside it and read what is produced as normal model
        ### for simplicity we use pandas to store interim file
        with open(primers['Example_model'], 'r') as file:
            F = file.readlines()
        indices = [i for i,val in enumerate(F) if val.startswith('#')][:4]
        VDG = []
        shuffle_=  {g:shuffle( range(0, len(list(fragments[g].items())) )) for g in fragments.keys()}
        
        VDG += [ ['%'+';'.join(x)+';{}\n'.format(i) for x,i in zip(  list(fragments[g].items()),
                                                                  shuffle_[g] ) ] for g in fragments.keys()]
        assert len(VDG)==len(indices)-1 
        F1 = F[0:indices[0]+1] + VDG[0] + [F[indices[1]]] + VDG[1] + [F[indices[2]]] + VDG[2] + F[indices[3]:]
#--------------------------
## Substitute dimentions in edges
        F2 = re.sub(r'(GeneChoice_([V|D|J])_gene_Undefined_side_prio[6|7]_size)(\d+)', lambda x:
                              r'{}{}'.format( x.group(1), len(fragments[x.group(2)])), ''.join(F1))
        return F2


############ SIMPLE UTILS #########################

### Some Util functions
def JSD(dist1,dist2, offset=10e-13):
        dist1 = dist1 + offset
        dist2 = dist2 + offset
        return np.sum(dist1*np.log(2*dist1/(dist1+dist2))+dist2*np.log(2*dist2/(dist1+dist2)))
    
def NestedDictValues(d):
    all_vals = []
    for v in d.values():
        if isinstance(v, dict):
            all_vals += NestedDictValues(v)
        else:
            all_vals += v
    return all_vals

def NestedDictKeys(d):
    inner_keys = []
    if all(isinstance(v, dict) for v in d.values()):
        for v in d.values():            
            inner_keys += NestedDictKeys(v)
    else:
        inner_keys += d.keys()
    return inner_keys

def unNestDict(d):
    inner_items = {}
    if all(isinstance(v, dict) for v in d.values()):
        for v in d.values():
            inner_items.update(unNestDict(v))
    else:
        inner_items.update(d)
    return inner_items

def map_nested_dicts(ob, func):
    if isinstance(ob, Mapping):
        return {str(k): map_nested_dicts(v, func) for k, v in ob.items()}
    else:
        return func(ob)












