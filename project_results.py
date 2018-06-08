#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 20:51:30 2018

@author: lunar
"""
import os
import json

os.chdir('/home/lunar/Desktop/IgoR_crawling/CODE_Project/')
from load_data import *
from utils import *

class project_models_results():
    def __init__(self, project_class):
        self.species = project_class.species
        self.batches = project_class.batches
        self.models = project_class.derived_models
        self.Outputs = project_class.Outputs
        ## test model here is just to get familiarized with syntax.
        # marginals - for probabilities
        # events - get events realization - for indices
        self.Test_model = self.models[self.species[0]][self.batches[self.species[0]][0]]
        self.gene_events = ['v_choice', 'd_gene', 'j_choice', 'direct_d_gene']
        self.indel_events = list(self.Test_model.marginals[0].keys())[3:]
        ### events probabilities are stored as Events-> Species -> Batches -> probs, so we can quickly 
        # combine together (bundle) all probabilities for one particular event. Probabilities are stored
        # as pandas DataFrames, because they can be ordered randomly for each batch.
        self.events_probs = {}
        self.sel = None
        # self.flattened_probs = {}
    
    def get_gene_choice_v_or_j(self, i):
        ## i=0 - v_choice, i=2 - j_choice
        eve = self.gene_events[i]
        vj_choice={}
        for sp in self.species:
            vj_choice.update({sp:{}})
            for batch in self.batches[sp]:
                model = self.models[sp][batch]
                vj_choice[sp].update({ re.match('(.*_(Out|In))_\d+', batch).group(1) : pd.Series(model.marginals[0][eve], 
                                index=model.events[i].get_realization_vector()) })
            vj_choice[sp]=pd.DataFrame(vj_choice[sp])
        self.events_probs.update({eve:vj_choice})
        
    def get_gene_choice_d(self):
        i=1
        eve = self.gene_events[i]
        d_choice={}
        for sp in self.species:
            d_choice.update({sp:{}})
            for batch in self.batches[sp]:
                model = self.models[sp][batch]
                d_choice[sp].update({ re.match('(.*_(Out|In))_\d+', batch).group(1) : pd.DataFrame(model.marginals[0][eve], columns=model.events[i].get_realization_vector(),
                                 index=model.events[i+1].get_realization_vector()) })
        self.events_probs.update({eve:d_choice})
        
    def get_direct_d_gene_probs(self):
        i=3
        eve = self.gene_events[i]
        d_direct_choice={}
        for sp in self.species:
            d_direct_choice.update({sp:{}})
            for batch in self.batches[sp]:
                model = self.models[sp][batch]
                d_direct_choice[sp].update({ re.match('(.*_(Out|In))_\d+', batch).group(1) : pd.Series(model.marginals[0]['d_gene'].T.dot(model.marginals[0]['j_choice']),
                                   model.events[1].get_realization_vector()) })
            d_direct_choice[sp]=pd.DataFrame(d_direct_choice[sp])
        self.events_probs.update({eve:d_direct_choice})
    
    def get_del_probs(self, i, eve):
        ## i determines on what positions deletions are considered
        r_deletions={}
        for sp in self.species:
            r_deletions.update({sp:{}})
            r_del_b = {}
            for batch in self.batches[sp]:
                model = self.models[sp][batch]
                
                r_del_b.update({ re.match('(.*_(Out|In))_\d+', batch).group(1) : model.marginals[0][eve].flatten()})
#                 r_deletions[sp].update({ re.match('(.*_(Out|In))_\d+', batch).group(1) : pd.DataFrame( model.marginals[0][eve],
#                                    model.events[i].get_realization_vector()).T })
            r_deletions[sp]=pd.DataFrame(r_del_b)
        self.events_probs.update({eve:r_deletions})    
    
    def get_ins_dinuc(self, i):
        ## i=4 - vd_ins, i=6 - dj_ins
        ## i=5 - vd_dinuc, i=7 - dj_dinuc
        ## all are shaped as vecotr columns
        eve = self.indel_events[i]
        Event={}
        for sp in self.species:
            Event.update({sp:{}})
            for batch in self.batches[sp]:
                model = self.models[sp][batch]
                Event[sp].update({ re.match('(.*_(Out|In))_\d+', batch).group(1) : pd.Series(model.marginals[0][eve], 
                                                            np.arange(len(model.marginals[0][eve])))})
            Event[sp]=pd.DataFrame(Event[sp])
        self.events_probs.update({eve:Event})   


    def get_events_probs(self):
        self.get_gene_choice_v_or_j(0)
        self.get_direct_d_gene_probs()
        self.get_gene_choice_v_or_j(2)
        self.get_gene_choice_d()
        for i, eve in zip([0,1,1,2], self.indel_events[:4]): 
            if i!=None:
                self.get_del_probs(i, eve)
        self.get_gene_choice_v_or_j
        self.get_ins_dinuc(4)
        self.get_ins_dinuc(5)
        self.get_ins_dinuc(6)
        self.get_ins_dinuc(7)

    def bundle_together(self, event, ds=None):
        if event:
            return reduce(lambda left, right: pd.merge(left,right, left_index=True, right_index=True), list(self.events_probs[event].values()))
        else:
            if isinstance(ds[0], pd.DataFrame):
                return reduce(lambda left, right: pd.merge(left,right, left_index=True, right_index=True),ds)
            elif isinstance(ds[0], dict):
                return NestedDictValues(ds)
    
    def compute_pairwise_distance(self, df, func=JSD):
        ## Compute JSD divergence. Symmetric metric allowing us to compare probability distributions
        JSD_res = {}
        for c1 in df.columns:
            c = {}
            for c2 in df.columns:
                c.update( {c2: func(df[c1],df[c2])} )
            JSD_res.update( {c1: pd.Series(c)} )
        return pd.DataFrame(JSD_res)
    
    def store_unstacked_results(self, list_of_labels=['v_choice', 'direct_d_gene', 'j_choice', 'vd_ins', 'dj_ins', 'dj_dinucl', 'vd_dinucl']):
        res = []
        for label in list_of_labels:
            H = self.compute_pairwise_distance(self.bundle_together(label))
            H = H.mask( np.tril(np.ones(H.shape, bool),False))
            res.append(pd.DataFrame(H.stack().values, index=[':'.join(x) for x in H.stack().index], columns=[label]))
        self.combined_unstacked = reduce(lambda left, right: pd.merge(left,right, left_index=True, right_index=True), res)
    
    def compute_JSD_averaged(self, a, b):
        ## drop zero columns 
        drop_a = a.T[a.sum(axis=0)==0].index
        drop_b = b.T[b.sum(axis=0)==0].index
        drop = set(list(drop_a)+ list(drop_b))
        a = a.drop(drop, axis=1)
        b = b.drop(drop, axis=1)
        res = []
        for c1, c2 in itertools.combinations( a.columns, 2):
            res.append(JSD(a[c1],b[c2]))
        res_m = np.mean(res)
        return res_m             
    
    def store_unstacked_results_2(self, list_of_labels = ['v_3_del', 'd_5_del', 'd_3_del', 'j_5_del'], dist_f=JSD):
        reses = []
        for label in list_of_labels:
            res = {}
            D = unNestDict(self.events_probs[label])
            if dist_f==JSD:
                for c1,c2 in [(ind.split(':')) for ind in self.combined_unstacked.index]:
                    res.update({':'.join( (c1,c2) ):(dist_f( D[c1].flatten() ,D[c2].flatten() ))})
            else:    
                for c1,c2 in [(ind.split(':')) for ind in self.combined_unstacked.index]:
                    res.update({':'.join( (c1,c2) ):(dist_f(D[c1],D[c2]))})
            self.combined_unstacked=self.combined_unstacked.join(pd.DataFrame(res, index=[label]).T)
                    
    def make_distance_matrix(self):
        batches_arr = NestedDictValues(self.batches)
#         cols =  list(map( lambda x: ':'.join(x),  itertools.combinations(NestedDictValues(RES.batches),2 )))
#         df = pd.DataFrame( np.zeros( [len(batches_arr),len(batches_arr)]),  columns=batches_arr, index=batches_arr)  
        df = pd.DataFrame( np.triu(np.ones([len(batches_arr),len(batches_arr)],bool),True) , columns=batches_arr, index=batches_arr)  
        return df
        
    def prepare_gen_intersect(self, batch1, batch2, affix1, affix2=None):

        def ret_sp(batch):
            for sp in self.batches.keys():
                if batch in self.batches[sp]:
                    return sp

        if not affix2:
            affix2=affix1
        path_to_gen1 = self.Outputs[ret_sp(batch1)]+batch1+'_generated/'+batch1[:-10]+affix1
        path_to_gen2 = self.Outputs[ret_sp(batch2)]+batch2+'_generated/'+batch2[:-10]+affix2
        return path_to_gen1, path_to_gen2       
    
    @staticmethod    
    def flatten_df(df):
        if type(df)==pd.Series or type(df)==pd.DataFrame:
            return df.values.flatten()
        else:
            return df.flatten()
    
    def make_flattened(self, event_list=None):
        try:
            del self.events_probs['flatten']
        except:
            pass
        flatten = {}
        if not event_list:
            event_list=list(self.events_probs.keys())
        for sp in self.events_probs['v_choice'].keys():
                flatten[sp]={}
                fl_ = {}
                for batch in self.events_probs['v_choice'][sp]:
                    vals = [self.events_probs[event][sp][batch] for event in event_list]
                    vals = list(map(self.flatten_df, vals))
                    vals = np.concatenate(vals)
                    fl_.update({batch: vals})
                fl_ = pd.DataFrame(fl_)
                flatten[sp] = fl_
        self.events_probs.update({'flatten':flatten})


    def save_selection(self,):
        ks = list(self.sel.keys())
        for k in ks:
            print(f'storing Model {k} on work table')
            json.dump( map_nested_dicts(self.sel[k], lambda x:x), open(k+'model_dump.json', 'w'))


    def select_model(self, batch):      
        model_dict = {}
        for eve in self.events_probs.keys():
            for sp in self.species:
                try: 
                    k= self.events_probs[eve][sp][batch].to_dict()
                    model_dict.update({eve:k})
                except:
                    pass
        print(model_dict.keys())
        self.sel={batch:model_dict}
        
    def plot_single(self, event, batch):
        spp = None
        for sp, B in self.batches.items():
            for b in B: 
                if re.match(batch, b):
                    spp = sp
                    
        if event=='v_3_del':
            a = self.events_probs[event][spp][batch]
            a = pd.DataFrame(a, columns=range(-4,17), index= self.events_probs['v_choice'][spp][batch].index)
            a.T.boxplot(rot=75,  figsize=(10,5))
            plt.title(batch)

                
        
        
        
        
        
