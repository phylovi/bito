import numpy as np
import pandas as pd


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


class AdaptiveStepsizeOptimizer:
    def __init__(self, model):
        self.model = model
        self.trace = []
        self.step_number = 0
        self.window_size = 5
        self.step_size = model.suggested_step_size()
        self.stepsize_increasing_rate = 1.2
        self.stepsize_decreasing_rate = 1 - 1e-2
        # The amount by which we drop the stepsize after realizing that it's gotten too
        # big.
        self.stepsize_drop_from_peak = 2
        self.stepsize_increasing = True
        self.best_elbo = -np.inf
        self.best_param_matrix = np.zeros(model.param_matrix.shape)
        self.optimizer = SGD_Server({"params": model.param_matrix.shape})

    def turn_around(self):
        """Triggered when the stepsize has gotten too big or a gradient step
        fails, restoring the previously best seen parameters."""
        np.copyto(self.model.param_matrix, self.best_param_matrix)
        self.step_size /= self.stepsize_drop_from_peak
        self.stepsize_increasing = False

    def gradient_step(self, target_log_like, grad_target_log_like):
        # Are we starting to decrease our objective function?
        if self.stepsize_increasing and self.step_number >= 2 * self.window_size:
            last_epoch = self.trace[-self.window_size :]
            prev_epoch = self.trace[-2 * self.window_size : -self.window_size]
            if np.mean(last_epoch) < np.mean(prev_epoch):
                self.turn_around()
        # Perform gradient step.
        if self.stepsize_increasing:
            self.step_size *= self.stepsize_increasing_rate
        else:
            self.step_size *= self.stepsize_decreasing_rate
        self.model.sample_and_prep_gradients()
        if not self.model.gradient_step(
            self.optimizer, self.step_size, grad_target_log_like(self.model.z)
        ):
            self.turn_around()  # Gradient step failed.
        self.trace.append(self.model.elbo_estimate(target_log_like, particle_count=500))
        if self.trace[-1] > self.best_elbo:
            self.best_elbo = self.trace[-1]
            np.copyto(self.best_param_matrix, self.model.param_matrix)
        self.step_number += 1

    def gradient_steps(self, target_log_like, grad_target_log_like, step_count):
        for _ in range(step_count):
            self.gradient_step(target_log_like, grad_target_log_like)
