import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp

import optimizers

import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp

tf.enable_v2_behavior()
tfd = tfp.distributions


class TFContinuousParameterModel:
    def __init__(
        self, q_factory, initial_params, branch_count, particle_count, step_size=0.01
    ):
        assert initial_params.ndim == 1
        self.q_factory = q_factory
        self.param_matrix = np.full((branch_count, len(initial_params)), initial_params)
        self.optimizer = optimizers.SGD_Server({"params": self.param_matrix.shape})
        self.particle_count = particle_count
        self.step_size = step_size
        # The current stored sample.
        self.x = None
        # The gradient of x with respect to the parameters of q.
        self.grad_x = None
        # The stochastic gradient of log sum q for x.
        self.grad_log_sum_q = None

    @property
    def branch_count(self):
        return self.param_matrix.shape[0]

    @property
    def param_count(self):
        return self.param_matrix.shape[1]

    def sample(self):
        with tf.GradientTape(persistent=True) as g:
            tf_params = tf.constant(self.param_matrix, dtype=tf.float32)
            g.watch(tf_params)
            q_distribution = self.q_factory(tf_params)
            tf_x = q_distribution.sample(self.particle_count)
            q_term = tf.math.reduce_sum(tf.math.log(q_distribution.prob(tf_x)))
        self.x = tf_x.numpy()
        # The Jacobian is laid out as particles x edges x edges x params.
        self.grad_x = np.sum(g.jacobian(tf_x, tf_params).numpy(), axis=2)
        self.grad_log_sum_q = g.gradient(q_term, tf_params).numpy()
        del g  # Should happen anyway but being explicit to remember.
        return self.x

    def clear_sample(self):
        self.x = None
        self.grad_x = None
        self.grad_log_sum_q = None

    @staticmethod
    def _chain_rule(grad_log_p_x, grad_x):
        return (
            np.tensordot(grad_log_p_x.transpose(), grad_x, axes=1)
            .diagonal()
            .transpose()
        )

    @staticmethod
    def _slow_chain_rule(grad_log_p_x, grad_x):
        particle_count, branch_count, param_count = grad_x.shape
        result = np.zeros((branch_count, param_count))
        for branch in range(branch_count):
            for param in range(param_count):
                for particle in range(particle_count):
                    result[branch, param] += (
                        grad_log_p_x[particle, branch] * grad_x[particle, branch, param]
                    )
        return result

    def linspace_one_branch(self, branch, min_x, max_x, num=50, default=0.1):
        """Fill a num x branch_count array with default, except for the branch
        column, which gets filled with a linspace."""
        a = np.full((num, self.branch_count), default)
        a[:, branch] = np.linspace(min_x, max_x, num)
        return a

    def elbo_gradient_using_current_sample(self, grad_log_p_x):
        assert self.grad_x is not None
        if not np.allclose(
            self._chain_rule(grad_log_p_x, self.grad_x),
            self._slow_chain_rule(grad_log_p_x, self.grad_x),
        ):
            print("chain rule isn't close")
            print(grad_log_p_x)
            print(self.grad_x)
            print(
                self._chain_rule(grad_log_p_x, self.grad_x)
                - self._slow_chain_rule(grad_log_p_x, self.grad_x)
            )
        unnormalized_result = (
            self._chain_rule(grad_log_p_x, self.grad_x) - self.grad_log_sum_q
        )
        return unnormalized_result / self.particle_count

    def gradient_step(self, grad_log_p_x, history=None):
        grad = self.elbo_gradient_using_current_sample(grad_log_p_x)
        update_dict = self.optimizer.adam(
            {"params": self.step_size}, {"params": self.param_matrix}, {"params": grad}
        )
        self.param_matrix += update_dict["params"]
        self.clear_sample()
        if history is not None:
            history.append(self.param_matrix.copy())

    def plot_1d(self, ax, target_log_like, which_branch, max_x=0.5):
        min_x = max_x / 100
        x_vals = self.linspace_one_branch(which_branch, min_x, max_x, 100)
        q_distribution = self.q_factory(self.param_matrix)
        df = pd.DataFrame(
            {
                "x": x_vals[:, which_branch],
                "target": sp.special.softmax(target_log_like(x_vals)),
                "fit": sp.special.softmax(
                    q_distribution.log_prob(x_vals).numpy()[:, which_branch]
                ),
            }
        )
        return df.plot(
            ax=ax,
            x="x",
            y=["target", "fit"],
            kind="line",
            title=q_distribution._name + " " + str(self.param_matrix[which_branch, :]),
        )

    def plot(self, target_log_like, max_x=0.5):
        f, axarr = plt.subplots(
            self.branch_count, sharex=True, figsize=(8, 1.5 * self.branch_count)
        )
        for which_branch in range(self.branch_count):
            self.plot_1d(axarr[which_branch], target_log_like, which_branch, max_x)
        plt.tight_layout()
