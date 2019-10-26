import abc
import numpy as np

import vip.priors
import vip.scalar_model


class BranchModel(abc.ABC):
    def __init__(self, scalar_model_name, inst):
        self.get_raw_representation = inst.get_psp_indexer_representations
        self.scalar_model = vip.scalar_model.of_name(
            scalar_model_name, variable_count=self._compute_variable_count(inst)
        )
        self.log_prior = vip.priors.log_exp_prior
        self.grad_log_prior = vip.priors.grad_log_exp_prior

    @staticmethod
    @abc.abstractmethod
    def _compute_variable_count(inst):
        """Compute the number of variables needed for a given instance.

        This is a private method that is just used in the initialization
        of the BranchModel.
        """
        pass


class SplitModel(BranchModel):
    def __init__(self, scalar_model_name, inst):
        super().__init__(scalar_model_name, inst)

    @staticmethod
    def _compute_variable_count(inst):
        return inst.psp_indexer.details()["after_rootsplits_index"]

    def px_branch_representation(self):
        """The ith entry of this array gives the index of the split
        corresponding to the ith branch."""
        return [
            np.array(representation[0])
            for representation in self.get_raw_representation()
        ]

    def mode_match(self, modes):
        self.scalar_model.mode_match(modes)

    def sample(self, particle_count, px_branch_representation):
        return self.scalar_model.sample(particle_count, px_branch_representation)

    def log_prob(self, theta_sample, px_branch_representation):
        return sum(
            self.scalar_model.log_prob(
                theta_sample[particle_idx, :], which_variables=branch_to_split
            )
            for particle_idx, branch_to_split in enumerate(px_branch_representation)
        )

    def sample_and_gradients(self, particle_count, px_branch_representation):
        return self.scalar_model.sample_and_gradients(
            particle_count, px_branch_representation
        )

    def scalar_grad(
        self, theta_sample, branch_gradients, px_branch_to_split, dg_dpsi, dlog_qg_dpsi
    ):
        """Do a gradient for the scalar parameters in terms of splits.

        See the tex for details. Comments of the form eq:XX refer to
        equations in the tex. See class-level docstring for information
        about `px_`.
        """
        # Calculate the gradient of the log unnormalized posterior.
        dlogp_dtheta = np.zeros_like(theta_sample)
        for particle_idx, (_, log_grad_raw) in enumerate(branch_gradients):
            # This :-2 is because of the two trailing zeroes that appear at the end of
            # the gradient.
            dlogp_dtheta[particle_idx, :] = np.array(log_grad_raw, copy=False)[:-2]
        dlogp_dtheta += self.grad_log_prior(theta_sample)
        # Now build up the complete gradient with respect to the variables.
        grad = np.zeros(
            (self.scalar_model.variable_count, self.scalar_model.param_count)
        )
        for particle_idx, branch_to_split in enumerate(px_branch_to_split):
            for branch_idx, variable_idx in enumerate(branch_to_split):
                grad[variable_idx, :] += (
                    # eq:dLdPsi
                    dlogp_dtheta[particle_idx, branch_idx]
                    * dg_dpsi[particle_idx, variable_idx, :]
                    - dlog_qg_dpsi[particle_idx, variable_idx, :]
                )
        return grad


def of_name(branch_model_name, scalar_model_name, inst):
    choices = {"split": SplitModel}
    if branch_model_name in choices:
        choice = choices[branch_model_name]
    else:
        raise Exception(f"BranchModel {branch_model_name} not known.")
    return choice(scalar_model_name, inst)
