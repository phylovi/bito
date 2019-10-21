import abc
import numpy as np


class BranchModel(abc.ABC):
    def __init__(self, get_raw_representation, scalar_model):
        self.get_raw_representation = get_raw_representation
        self.scalar_model = scalar_model
        # TODO make the scalar model fixed


class SplitModel(BranchModel):
    def __init__(self, get_raw_representation, scalar_model):
        super().__init__(get_raw_representation, scalar_model)


def of_name(name, get_raw_representation, scalar_model):
    choices = {"split": SplitModel}
    if name in choices:
        choice = choices[name]
    else:
        raise Exception(f"BranchModel {name} not known.")
    return choice(get_raw_representation, scalar_model)
