""" Class for SBN modeling. The code is based on a previous version of libsbn. 
   To locate it, try: 

   git checkout 19898c9cdb349690b46914988178a18d933d46ff

"""

import numpy as np
from utils import logsumexp, softmax, softmax_parser, split_sum, aggregate_sbn


class SBNModel:
    EPS = 1e-40
    def __init__(self, inst):
        self.inst = inst
        self.sbn_parameters = np.array(self.inst.sbn_parameters, copy=False)
        self.num_sbn_parameters = np.size(self.sbn_parameters)
        self._sbn_parameters = np.zeros(self.num_sbn_parameters)
        index_dict, range_dict = self.inst.get_indexers()
        self.param_range = np.array(sorted(range_dict.values(), key=lambda x:x[0]))
        self.sbn_parameters[:] = softmax_parser(self._sbn_parameters, self.param_range)
        
    def sample_trees(self, count):
        self.inst.sample_trees(count)
    
    def tree_loglikelihood(self, grad=False, value_and_grad=False):
        rootsplit_index_list, subsplit_index_list = zip(*self.inst.make_indexer_representations())
        sbn_index_list = np.concatenate((np.expand_dims(np.array(rootsplit_index_list), axis=-1), np.array(subsplit_index_list)), axis=-1)
        CPDs = np.append(self.sbn_parameters, 0.0)
        prob_mat = CPDs[sbn_index_list].clip(self.EPS)
        loglikelihood = logsumexp(np.sum(np.log(prob_mat), axis=-1), axis=-1)
        if not grad:
            return loglikelihood
        
        wts = softmax(np.sum(np.log(prob_mat), axis=-1), axis=-1)
        CPDs_grad = np.zeros((self.inst.tree_count(), self.num_sbn_parameters))
        aggregate_sbn(np.expand_dims(wts, -1)/prob_mat, sbn_index_list, CPDs_grad, (self.num_sbn_parameters, ))
        
        range_sum_grad = np.stack([split_sum(self.sbn_parameters*grad, self.param_range) for grad in CPDs_grad])
        CPDs_grad = (CPDs_grad - range_sum_grad) * self.sbn_parameters
        
        if not value_and_grad:
            return CPDs_grad
        else:
            return loglikelihood, CPDs_grad