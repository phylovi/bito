"""Classes that perform optimization and contain relevant state."""

import abc
import numpy as np
import vip.sgd_server
from vip.sbn_model import SBNModel
from vip.scalar_model import ScalarModel


class BaseOptimizer(abc.ABC):
    def __init__(
        self, sbn_model: SBNModel, scalar_model: ScalarModel, elbo_estimator_fun
    ):
        self.sbn_model = sbn_model
        self.scalar_model = scalar_model
        self.estimate_elbo = elbo_estimator_fun
        self.trace = []
        self.step_number = 0
        self.step_size = scalar_model.suggested_step_size()
        self.sbn_step_size = 0.001
        self.sgd_server = vip.sgd_server.SGD_Server(
            {
                "scalar_params": scalar_model.q_params.shape,
                "sbn_params": sbn_model.sbn_parameters.shape,
            }
        )

    def _simple_gradient_step(self, grad_dict, history=None):
        """Just take a simple gradient step.

        scalar_grad is the gradient with respect to the scalar variables.
        sbn_grad is the gradient with respect to the SBN variables.

        Return True if the gradient step was successful.
        """
        scalar_grad = grad_dict["scalar_params"]
        sbn_grad = grad_dict["sbn_params"]
        assert self.scalar_model.q_params.shape == scalar_grad.shape
        if not np.isfinite(np.array([scalar_grad])).all():
            return False
        assert self.sbn_model.sbn_parameters.shape == sbn_grad.shape
        update_dict = self.sgd_server.adam(
            {"scalar_params": self.step_size, "sbn_params": self.sbn_step_size},
            {
                "scalar_params": self.scalar_model.q_params,
                "sbn_params": self.sbn_model.sbn_parameters,
            },
            grad_dict,
        )
        self.scalar_model.q_params += update_dict["scalar_params"]
        self.sbn_model.sbn_parameters += update_dict["sbn_params"]
        if history is not None:
            history.append(self.scalar_model.q_params.copy())
            history.append(self.sbn_model.sbn_parameters.copy())
        return True

    def gradient_step(self, grad_dict, history=None):
        gradient_step_was_successful = self._simple_gradient_step(grad_dict)
        self.update(gradient_step_was_successful)

    @abc.abstractmethod
    def update(self, gradient_step_was_successful):
        """Update the parameters of the optimizer based on any information it
        stores."""
        pass


class SimpleOptimizer(BaseOptimizer):
    def __init__(
        self, sbn_model: SBNModel, scalar_model: ScalarModel, elbo_estimator_fun
    ):
        super().__init__(sbn_model, scalar_model, elbo_estimator_fun)
        self.stepsize_decreasing_rate = 1 - 1e-2

    def update(self, gradient_step_was_successful):
        if gradient_step_was_successful:
            self.step_size *= self.stepsize_decreasing_rate
        else:
            self.step_size /= 2
        self.step_number += 1


class BumpStepsizeOptimizer(BaseOptimizer):
    """An optimizer that increases the stepsize until it's too big, then
    decreases it."""

    def __init__(
        self, sbn_model: SBNModel, scalar_model: ScalarModel, elbo_estimator_fun
    ):
        super().__init__(sbn_model, scalar_model, elbo_estimator_fun)
        self.window_size = 5
        self.stepsize_increasing_rate = 1.2
        self.stepsize_decreasing_rate = 1 - 1e-2
        # The amount by which we drop the stepsize after realizing that it's gotten too
        # big.
        self.stepsize_drop_from_peak = 4
        self.stepsize_increasing = True
        self.best_elbo = -np.inf
        self.best_q_params = np.zeros(scalar_model.q_params.shape)

    def _turn_around(self):
        """Triggered when the stepsize has gotten too big or a gradient step
        fails, restoring the previously best seen parameters."""
        np.copyto(self.scalar_model.q_params, self.best_q_params)
        self.step_size /= self.stepsize_drop_from_peak
        self.stepsize_increasing = False

    def update(self, gradient_step_was_successful):
        if not gradient_step_was_successful:
            self._turn_around()
        # Are we starting to decrease our objective function?
        if self.stepsize_increasing and self.step_number >= 2 * self.window_size:
            last_epoch = self.trace[-self.window_size :]
            prev_epoch = self.trace[-2 * self.window_size : -self.window_size]
            if np.mean(last_epoch) < np.mean(prev_epoch):
                self._turn_around()
        # Adjust step size.
        if self.stepsize_increasing:
            self.step_size *= self.stepsize_increasing_rate
        else:
            self.step_size *= self.stepsize_decreasing_rate
        self.trace.append(self.estimate_elbo(particle_count=500))
        if self.trace[-1] > self.best_elbo:
            self.best_elbo = self.trace[-1]
            np.copyto(self.best_q_params, self.scalar_model.q_params)
        self.step_number += 1
        return np.isfinite(self.trace[-1])


def of_name(name, sbn_model: SBNModel, scalar_model: ScalarModel, elbo_estimator_fun):
    choices = {"simple": SimpleOptimizer, "bump": BumpStepsizeOptimizer}
    if name in choices:
        choice = choices[name]
    else:
        raise Exception(f"Optimizer {name} not known.")
    return choice(sbn_model, scalar_model, elbo_estimator_fun)
