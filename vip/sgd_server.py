import numpy as np


class SGD_Server(object):
    """A class for stochastic gradient descent.

    Methods
    -------
    SGD (with or without momentum)
    Adam
    AMSGrad
    RMSProp
    ADAGrad
    ADADelta
    """

    def __init__(
        self,
        dim_dict,
        beta_0=0.9,
        beta_1=0.999,
        beta_1_ams=0.99,
        gamma=0.9,
        epsilon=1e-08,
        decay=0.0,
        momentum=0.9,
    ):
        self.beta_0, self.beta_1 = beta_0, beta_1
        self.beta_1_ams = beta_1_ams
        self.gamma = 0.9
        self.mom = momentum
        self.decay = decay
        self.eps = epsilon
        self.vars = dim_dict.keys()
        self.mean_grad = {var: np.zeros(d) for var, d in dim_dict.items()}
        self.var_grad = {var: np.zeros(d) for var, d in dim_dict.items()}
        self.var_delta = {var: np.zeros(d) for var, d in dim_dict.items()}
        self.var_grad_max = {var: np.zeros(d) for var, d in dim_dict.items()}
        self.t = 0

    def sgd(self, stepsz_dict, param_dict, grad_dict):
        for var in self.vars:
            grad = grad_dict[var] - self.decay * param_dict[var]
            self.mean_grad[var] = (
                self.mom * self.mean_grad[var] + stepsz_dict[var] * grad
            )
        return self.mean_grad

    def adam(self, stepsz_dict, param_dict, grad_dict):
        self.t += 1
        update_dict = {}
        for var in self.vars:
            grad = grad_dict[var] - self.decay * param_dict[var]
            self.mean_grad[var] = (
                self.beta_0 * self.mean_grad[var] + (1.0 - self.beta_0) * grad
            )
            self.var_grad[var] = (
                self.beta_1 * self.var_grad[var] + (1.0 - self.beta_1) * grad ** 2
            )
            hat_mean_grad = self.mean_grad[var] / (1.0 - self.beta_0 ** self.t)
            hat_var_grad = self.var_grad[var] / (1.0 - self.beta_1 ** self.t)
            update_dict[var] = (
                stepsz_dict[var] * hat_mean_grad / (np.sqrt(hat_var_grad) + self.eps)
            )
        return update_dict

    def amsgrad(self, stepsz_dict, param_dict, grad_dict):
        self.t += 1
        update_dict = {}
        for var in self.vars:
            grad = grad_dict[var] - self.decay * param_dict[var]
            self.mean_grad[var] = (
                self.beta_0 * self.mean_grad[var] + (1.0 - self.beta_0) * grad
            )
            self.var_grad[var] = (
                self.beta_1_ams * self.var_grad[var]
                + (1.0 - self.beta_1_ams) * grad ** 2
            )
            self.var_grad_max[var] = np.maximum(
                self.var_grad_max[var], self.var_grad[var]
            )
            hat_mean_grad = self.mean_grad[var] / (1.0 - self.beta_0 ** self.t)
            hat_var_grad = self.var_grad_max[var] / (1.0 - self.beta_1_ams ** self.t)
            update_dict[var] = (
                stepsz_dict[var] * hat_mean_grad / (np.sqrt(hat_var_grad) + self.eps)
            )
        return update_dict

    def rmsprop(self, stepsz_dict, param_dict, grad_dict):
        update_dict = {}
        for var in self.vars:
            grad = grad_dict[var] - self.decay * param_dict[var]
            self.var_grad[var] = (
                self.gamma * self.var_grad[var] + (1.0 - self.gamma) * grad ** 2
            )
            update_dict[var] = (
                stepsz_dict[var] * grad / np.sqrt(self.var_grad[var] + self.eps)
            )
        return update_dict

    def adagrad(self, stepsz_dict, param_dict, grad_dict):
        update_dict = {}
        for var in self.vars:
            grad = grad_dict[var] - self.decay * param_dict[var]
            self.var_grad[var] = self.var_grad[var] + grad ** 2
            update_dict[var] = (
                stepsz_dict[var] * grad / np.sqrt(self.var_grad[var] + self.eps)
            )
        return update_dict

    def adadelta(self, stepsz_dict, param_dict, grad_dict):
        update_dict = {}
        for var in self.vars:
            grad = grad_dict[var] - self.decay * param_dict[var]
            self.var_grad[var] = (
                self.gamma * self.var_grad[var] + (1.0 - self.gamma) * grad ** 2
            )
            update_dict[var] = (
                np.sqrt(
                    (self.var_delta[var] + self.eps) / (self.var_grad[var] + self.eps)
                )
                * grad
            )
            self.var_delta[var] = (
                self.gamma * self.var_delta[var]
                + (1.0 - self.gamma) * update_dict[var] ** 2
            )
        return update_dict
