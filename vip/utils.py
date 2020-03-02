""" Useful functions for efficient numpy array operations. """

import numpy as np
import numba
from numba.typed import List
from numba import njit


def logsumexp(arrA, axis=0, keepdims=False):
    max_arrA = np.max(arrA, axis=axis, keepdims=True)
    tmp = np.log(np.sum(np.exp(arrA - max_arrA), axis=axis, keepdims=keepdims))
    if not keepdims:
        max_arrA = max_arrA.squeeze(axis=axis)
    return  tmp + max_arrA
	

def softmax(arrA, axis=0):
    max_arrA = np.max(arrA, axis=axis, keepdims=True)
    exp_arr = np.exp(arrA - max_arrA)
    return exp_arr / np.sum(exp_arr, axis=axis, keepdims=True)
        

@njit    
def split_sum(arrA, range_arr):
    output = np.empty_like(arrA)
    for i in range(len(range_arr)):
        start, end = range_arr[i]
        output[start:end] = np.sum(arrA[start:end])
    return output

@njit    
def split_max(arrA, range_arr):
    output = np.empty_like(arrA)
    for i in range(len(range_arr)):
        start, end = range_arr[i]
        output[start:end] = np.max(arrA[start:end])
    return output
    
def softmax_parser(arrA, range_arr):
    max_arrA = split_max(arrA, range_arr)
    exp_arrA = np.exp(arrA - max_arrA)
    normal_arrA = split_sum(exp_arrA, range_arr)
    return exp_arrA / normal_arrA
    

@njit
def aggregate_sbn(matA, index_mat, out, exclude):
    k, m, n = index_mat.shape
    for l in range(k):
        for i in range(m):
            for j in range(n):
                if index_mat[l,i,j] not in exclude:
                    out[l, index_mat[l,i,j]] += matA[l, i, j]
