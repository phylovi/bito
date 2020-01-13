"""Classes for modeling phylogenetic branch lengths."""


import abc
import numpy as np

import vip.priors
import vip.scalar_model
from vip.scalar_model import LogNormalModel


class BranchModel(abc.ABC):
    """Our abstract base class for branch lengths."""

    def __init__(self, scalar_model_name, inst):
        self.make_raw_representation = inst.make_psp_indexer_representations
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

    @abc.abstractmethod
    def mode_match(self, split_modes):
        """Do some crazy heuristics for mode matching.

        split_modes is the set of modes as indexed by the splits.
        """
        pass

    @abc.abstractmethod
    def sample(self, px_branch_representation):
        """Sample from the branch model.

        Returns a list of samples over particles, each one being a
        sample indexed by branches according to
        px_branch_representation.
        """
        pass

    @abc.abstractmethod
    def sample_all(self, particle_count):
        """Sample all of the splits from the model.

        Only makes sense for split-based models.
        """
        pass


class SplitModel(BranchModel):
    """Branch modeling parameterized by splits."""

    def __init__(self, scalar_model_name, inst):
        super().__init__(scalar_model_name, inst)

    @staticmethod
    def _compute_variable_count(inst):
        return inst.psp_indexer.details()["after_rootsplits_index"]

    def px_branch_representation(self):
        """The split-based representation of each of the trees.

        The ith entry of the resulting list is the representation of the
        ith tree. The jth entry of such a representation is the index of
        the split corresponding to the ith branch.
        """
        return [
            np.array(representation[0])
            for representation in self.make_raw_representation()
        ]

    def mode_match(self, split_modes):
        self.scalar_model.mode_match(split_modes)

    def sample(self, px_branch_representation):
        return self.scalar_model.sample(px_branch_representation)

    def sample_all(self, particle_count):
        return self.scalar_model.sample_all(particle_count)

    def log_prob(self, theta_sample, px_branch_representation):
        return sum(
            self.scalar_model.log_prob(
                theta_sample[particle_idx, :], which_variables=branch_to_split
            )
            for particle_idx, branch_to_split in enumerate(px_branch_representation)
        )

    def sample_and_gradients(self, px_branch_representation):
        return self.scalar_model.sample_and_gradients(px_branch_representation)

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


class PSPModel(BranchModel):
    """The object containing a PSP model.

    See `psp_indexer.hpp` for a description of how the PSPs are indexed.
    Especially important: the after_rootsplits_index acts as a sentinel
    for when there isn't a component of a given PSP. We have an entry
    for that sentinel so we don't have to write special-purpose code,
    and instead keep that entry at zero so that it doesn't contribute
    anything. That way we can just do the usual summation across PSP
    components.
    """

    def __init__(self, scalar_model_name, inst):
        if scalar_model_name != "lognormal":
            raise Exception("PSP only works with LogNormal.")
        super().__init__(scalar_model_name, inst)
        details = inst.psp_indexer.details()
        # This is how things should be laid out.
        assert details["rootsplit_position"] == 0
        assert details["subsplit_down_position"] == 1
        assert details["subsplit_up_position"] == 2
        self.after_rootsplits_index = details["after_rootsplits_index"]
        # Here we aren't using the scalar_model as it was designed -- we are replacing
        # most of its functionality by PSPModel methods.
        self.q_params = self.scalar_model.q_params
        # Set the sentinel parameters to be zero.
        self.q_params[-1, :] = 0.0

    @staticmethod
    def _compute_variable_count(inst):
        # Here the +1 is because we have a sentinel.
        return inst.psp_indexer.details()["first_empty_index"] + 1

    def px_branch_representation(self):
        """The PSP-based representation of each of the trees."""
        return [
            np.array(representation)
            for representation in self.make_raw_representation()
        ]

    def mode_match(self, split_modes):
        assert split_modes.size == self.after_rootsplits_index
        self.q_params[:, :] = 0.0
        log_modes = np.log(np.clip(split_modes, 1e-6, None))
        biclipped_log_modes = np.log(np.clip(split_modes, 1e-6, 1 - 1e-6))
        # Issue #118: do we want to have the lognormal variance be in log space? Or do
        # we want to clip it?
        # Here we set the PSP parameters to zero except for the first block of
        # parameters, which are the splits.
        split_q_params = self.q_params[: self.after_rootsplits_index, :]
        split_q_params[:, 1] = -0.1 * biclipped_log_modes
        split_q_params[:, 0] = np.square(split_q_params[:, 1]) + log_modes

    def _make_lognormal_params(self, branch_representation):
        """Sum parameters across the branch representation so that we get a
        full lognormal parametrization."""
        branch_count = branch_representation.shape[1]
        lognormal_params = np.zeros((branch_count, 2))
        for psp_idx in range(3):
            lognormal_params[:, :] += self.q_params[
                branch_representation[psp_idx, :], :
            ]
        return lognormal_params

    def sample(self, px_branch_representation):
        assert len(px_branch_representation) > 0
        branch_representation_shape = px_branch_representation[0].shape
        branch_count = branch_representation_shape[1]
        px_sample = np.empty((len(px_branch_representation), branch_count))
        for particle_idx, branch_representation in enumerate(px_branch_representation):
            assert branch_representation_shape == branch_representation.shape
            # Set up a lognormal parameter matrix for this specific branch.
            lognormal_params = self._make_lognormal_params(branch_representation)
            px_sample[particle_idx, :] = np.random.lognormal(
                lognormal_params[:, 0], lognormal_params[:, 1]
            )
        return px_sample

    def sample_all(self, particle_count):
        # This is a placeholder, but it's hard to know what should go here for the PSP
        # case. This function is designed to compare things to the splits-based
        # distribution from an MCMC run. In principle we could boil the MCMC run down to
        # a summary of branch lengths per PSP, but there would be a bunch of coding on
        # that side too.
        return np.zeros((self.after_rootsplits_index, 1))

    def log_prob(self, theta_sample, px_branch_representation):
        total = 0.0
        for particle_idx, branch_representation in enumerate(px_branch_representation):
            lognormal_params = self._make_lognormal_params(branch_representation)
            total += LogNormalModel.general_log_prob(
                theta_sample[particle_idx, :],
                lognormal_params[:, 0],
                lognormal_params[:, 1],
            )
        return total

    def sample_and_gradients(self, px_branch_representation):
        particle_count = len(px_branch_representation)
        branch_representation_shape = px_branch_representation[0].shape
        sample = np.empty((particle_count, branch_representation_shape[1]))
        dg_dpsi = np.zeros((particle_count, self.scalar_model.variable_count, 2))
        dlog_qg_dpsi = np.zeros((particle_count, self.scalar_model.variable_count, 2))
        # eq:dlogqgdPsi
        dlog_qg_dpsi[:, :, 0] = -1.0
        for particle_idx, branch_representation in enumerate(px_branch_representation):
            assert branch_representation_shape == branch_representation.shape
            lognormal_params = self._make_lognormal_params(branch_representation)
            mu = lognormal_params[:, 0]
            sigma = lognormal_params[:, 1]
            sample[particle_idx, :] = np.random.lognormal(mu, sigma)
            # eq:gLogNorm
            epsilon = (np.log(sample[particle_idx, :]) - mu) / sigma
            # Loop over the rows of the branch representation, which are the root split
            # and the two PSPs. As described in the tex, we just have to set all of
            # their derivatives equal to the lognormal derivatives with respect to the
            # corresponding split.
            for which_variables in branch_representation:
                # eq:dgdPsi
                dg_dpsi[particle_idx, which_variables, 0] = sample[particle_idx, :]
                dg_dpsi[particle_idx, which_variables, 1] = (
                    sample[particle_idx, :] * epsilon
                )
                # eq:dlogqgdPsi
                dlog_qg_dpsi[particle_idx, which_variables, 1] = -epsilon - 1.0 / sigma
        return (sample, dg_dpsi, dlog_qg_dpsi)

    def scalar_grad(
        self,
        theta_sample,
        branch_gradients,
        px_branch_representation,
        dg_dpsi,
        dlog_qg_dpsi,
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
        for particle_idx, branch_representation in enumerate(px_branch_representation):
            for which_variables in branch_representation:
                for branch_idx, variable_idx in enumerate(which_variables):
                    grad[variable_idx, :] += (
                        # eq:dLdPsi
                        dlogp_dtheta[particle_idx, branch_idx]
                        * dg_dpsi[particle_idx, variable_idx, :]
                        - dlog_qg_dpsi[particle_idx, variable_idx, :]
                    )
        # We want our sentinel to stay zero (see docstring for this class).
        grad[-1, :] = 0.0
        return grad


def of_name(branch_model_name, scalar_model_name, inst):
    choices = {"split": SplitModel, "psp": PSPModel}
    if branch_model_name in choices:
        choice = choices[branch_model_name]
    else:
        raise Exception(f"BranchModel {branch_model_name} not known.")
    return choice(scalar_model_name, inst)
