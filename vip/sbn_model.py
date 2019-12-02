"""
Stub (unimplemented) class for SBN modeling.
"""
import numpy as np


class SBNModel:
    def __init__(self, inst):
        self.sbn_parameters = np.array(inst.sbn_parameters, copy=False)
