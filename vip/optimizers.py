import numpy as np
import vip.sgd_server
from vip.sbn_model import SBNModel
from vip.scalar_models import ScalarModel


class BaseOptimizer:
    def __init__(self, sbn_model: SBNModel, scalar_model: ScalarModel):
        self.sbn_model = sbn_model
        self.scalar_model = scalar_model
        self.trace = []
        self.step_number = 0
        self.step_size = scalar_model.suggested_step_size()
        self.sgd_server = vip.sgd_server.SGD_Server(
            {"scalar_params": scalar_model.q_params.shape}
        )

    def _simple_gradient_step(self, vars_grad, history=None):
        """Just take a simple gradient step.

        Return True if the gradient step was successful.
        """
        assert self.scalar_model.q_params.shape == vars_grad.shape
        if not np.isfinite(np.array([vars_grad])).all():
            self.scalar_model.clear_sample()
            return False
        update_dict = self.sgd_server.adam(
            {"scalar_params": self.step_size},
            {"scalar_params": self.scalar_model.q_params},
            {"scalar_params": vars_grad},
        )
        self.scalar_model.q_params += update_dict["scalar_params"]
        self.scalar_model.clear_sample()
        if history is not None:
            history.append(self.scalar_model.q_params.copy())
        return True


class SimpleOptimizer(BaseOptimizer):
    def __init__(self, sbn_model: SBNModel, scalar_model: ScalarModel):
        super().__init__(sbn_model, scalar_model)
        self.stepsize_decreasing_rate = 1 - 1e-2

    def gradient_step(self, target_log_like, grad_target_log_like, which_variables):
        """
        Take a gradient step.

        target_log_like and grad_target_log_like are all in terms of branch lengths.
        """
        self.scalar_model.sample_and_prep_gradients(which_variables)
        # This gradient is in terms of branch lengths...
        grad_log_p_z = grad_target_log_like(self.scalar_model.brlen_sample)
        # ... so we need to move it into variable space.
        vars_grad = np.zeros(
            (self.scalar_model.variable_count, self.scalar_model.param_count)
        )
        for branch_index, variable_index in enumerate(which_variables):
            # import pdb
            # pdb.set_trace()
            for param_index in range(self.scalar_model.param_count):
                vars_grad[variable_index, param_index] += (
                    np.sum(
                        grad_log_p_z[:, branch_index]
                        * self.scalar_model.grad_z[:, variable_index, param_index],
                        axis=0,
                    )
                    - self.scalar_model.grad_sum_log_q[variable_index, param_index]
                )
        if self._simple_gradient_step(vars_grad):
            self.step_size *= self.stepsize_decreasing_rate
        else:
            self.step_size /= 2
        self.step_number += 1
        return True


class BumpStepsizeOptimizer(BaseOptimizer):
    """An optimizer that increases the stepsize until it's too big, then
    decreases it."""

    def __init__(self, sbn_model: SBNModel, scalar_model: ScalarModel):
        super().__init__(sbn_model, scalar_model)
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

    def gradient_step(self, target_log_like, grad_target_log_like, which_variables):
        raise NotImplementedError("We don't use which_variables yet.")
        # Are we starting to decrease our objective function?
        if self.stepsize_increasing and self.step_number >= 2 * self.window_size:
            last_epoch = self.trace[-self.window_size :]
            prev_epoch = self.trace[-2 * self.window_size : -self.window_size]
            if np.mean(last_epoch) < np.mean(prev_epoch):
                self._turn_around()
        # Perform gradient step.
        if self.stepsize_increasing:
            self.step_size *= self.stepsize_increasing_rate
        else:
            self.step_size *= self.stepsize_decreasing_rate
        self.scalar_model.sample_and_prep_gradients()
        if not self._simple_gradient_step(
            grad_target_log_like(self.scalar_model.brlen_sample)
        ):
            self._turn_around()  # Gradient step failed.
        self.trace.append(
            self.scalar_model.elbo_estimate(target_log_like, particle_count=500)
        )
        if self.trace[-1] > self.best_elbo:
            self.best_elbo = self.trace[-1]
            np.copyto(self.best_q_params, self.scalar_model.q_params)
        self.step_number += 1
        return np.isfinite(self.trace[-1])


def of_name(name, sbn_model: SBNModel, scalar_model: ScalarModel):
    choices = {"simple": SimpleOptimizer, "bump": BumpStepsizeOptimizer}
    if name in choices:
        choice = choices[name]
    else:
        raise Exception(f"Optimizer {name} not known.")
    return choice(sbn_model, scalar_model)
