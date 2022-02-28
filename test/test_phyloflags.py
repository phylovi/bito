"""Some basic testing and demo code for the bito module.

If you want to see the results of the print statements, use `pytest -s`.
"""

import json
import pprint
import pytest
import numpy as np
import bito
import bito.beagle_flags as beagle_flags
import bito.phylo_flags as flags
import bito.phylo_model_mapkeys as model_keys
import bito.phylo_gradient_mapkeys as gradient_keys
import sys

# DEMO


def gradients_with_flags_demo():
    inst = create_instance()
    initialize_model_parameters(inst)

    # Request and calculate gradients for RATIOS_ROOT_HEIGHT,
    # SUBSTITUTION_MODEL_FREQUENCIES, SUBSTITUTION_MODEL_RATES
    bito_grad = inst.phylo_gradients(
        # explicit flags: For flags that don't have associated values, can just be a list.
        [
            flags.RATIOS_ROOT_HEIGHT,
            flags.SUBSTITUTION_MODEL,
        ],
        # run_with_default_flags:
        # (1) If set to true, all fields of phylo_model_block_map are populated unless overriden by explicit flag.
        # (2) If set to false, no fields of phylo_model_block_map are populated unless overriden by explicit flag.
        False,
    )[0]

    print("bito_grad_keys: ", bito_grad.gradient.keys())
    gtr_rates_grad = np.array(
        bito_grad.gradient[gradient_keys.SUBSTITUTION_MODEL_RATES]
    )
    gtr_freqs_grad = np.array(
        bito_grad.gradient[gradient_keys.SUBSTITUTION_MODEL_FREQUENCIES]
    )
    ratios_root_height = np.array(bito_grad.gradient[gradient_keys.RATIOS_ROOT_HEIGHT])

    print("GTR rates gradient: \n", gtr_rates_grad)
    print("GTR freqs gradient: \n", gtr_freqs_grad)
    print("root height gradient: \n", ratios_root_height)

    # above works if only boolean flags are used, otherwise:
    bito_grad = inst.phylo_gradients(
        # explicit flags: For SET flags that require value, use ordered tuples. Non-SET flags just take boolean.
        [
            (flags.SET_GRADIENT_DELTA, 5.0),
        ],
        # run_with_default_flags
        True,
    )[0]

    gtr_rates_grad = np.array(
        bito_grad.gradient[gradient_keys.SUBSTITUTION_MODEL_RATES]
    )
    gtr_freqs_grad = np.array(
        bito_grad.gradient[gradient_keys.SUBSTITUTION_MODEL_FREQUENCIES]
    )
    ratios_root_height = np.array(bito_grad.gradient[gradient_keys.RATIOS_ROOT_HEIGHT])
    print("GTR rates gradient: \n", gtr_rates_grad)
    print("GTR freqs gradient: \n", gtr_freqs_grad)

    # we can also use the internal phylo_flags.
    inst.init_phylo_flags()
    inst.set_phylo_defaults(True)
    inst.set_phylo_flag(flags.SET_GRADIENT_DELTA, 5.0)
    # these can be used with or without passing flags as arguments.
    bito_grad_2 = inst.phylo_gradients()[0]
    inst.clear_phylo_flags()

    gtr_rates_grad = np.array(
        bito_grad_2.gradient[gradient_keys.SUBSTITUTION_MODEL_RATES]
    )
    gtr_freqs_grad = np.array(
        bito_grad_2.gradient[gradient_keys.SUBSTITUTION_MODEL_FREQUENCIES]
    )
    ratios_root_height = np.array(
        bito_grad_2.gradient[gradient_keys.RATIOS_ROOT_HEIGHT]
    )

    # This should trigger an exception because CLOCK_MODEL_RATES was not flagged to be computed.
    try:
        clock_grad = np.array(bito_grad.gradient[gradient_keys.CLOCK_MODEL_RATES])
        print("CHECK_THROW Failed: Key error not caught.")
    except:
        print("CHECK_THROW Successful: Error successfully caught.")
        print("ERROR: ", sys.exc_info()[0], "occurred.")

    # above works if only boolean flags are used, otherwise:
    bito_grad = inst.phylo_gradients(
        # explicit flags: For SET flags that require value, use ordered tuples. Non-SET flags just take boolean.
        [(flags.SUBSTITUTION_MODEL, False)],
        # run_with_default_flags
        True,
    )[0]
    return


# TESTS


unflagged_keys = [gradient_keys.BRANCH_LENGTHS]

include_flags_to_keys = {
    flags.SITE_MODEL: [gradient_keys.SITE_MODEL],
    flags.CLOCK_MODEL: [gradient_keys.CLOCK_MODEL],
    flags.SUBSTITUTION_MODEL: [
        gradient_keys.SUBSTITUTION_MODEL,
        gradient_keys.SUBSTITUTION_MODEL_RATES,
        gradient_keys.SUBSTITUTION_MODEL_FREQUENCIES,
    ],
    flags.RATIOS_ROOT_HEIGHT: [gradient_keys.RATIOS_ROOT_HEIGHT],
}
for flag in include_flags_to_keys:
    for unflagged_key in unflagged_keys:
        include_flags_to_keys[flag].append(unflagged_key)

exclude_flags = [
    flags.INCLUDE_LOG_DET_JACOBIAN_LIKELIHOOD,
    flags.INCLUDE_LOG_DET_JACOBIAN_GRADIENT,
]

setvalue_flags = [flags.SET_GRADIENT_DELTA]


def create_instance():
    inst = bito.rooted_instance("cheese")
    inst.read_newick_file("data/fluA.tree")
    inst.read_fasta_file("data/fluA.fa")
    inst.parse_dates_from_taxon_names(True)
    spec = bito.PhyloModelSpecification(
        substitution="GTR", site="weibull+4", clock="strict"
    )
    inst.prepare_for_phylo_likelihood(spec, 1, [beagle_flags.VECTOR_SSE], False)
    return inst


def initialize_model_parameters(inst):
    # initialize model parameters (with enums)
    phylo_model_param_block_map = inst.get_phylo_model_param_block_map()
    phylo_keys = phylo_model_param_block_map.keys()
    phylo_model_param_block_map[model_keys.SUBSTITUTION_MODEL_RATES][:] = np.repeat(
        1 / 6, 6
    )
    phylo_model_param_block_map[model_keys.SUBSTITUTION_MODEL_FREQUENCIES][
        :
    ] = np.repeat(1 / 4, 4)
    phylo_model_param_block_map[model_keys.SITE_MODEL][:] = np.array([0.5])
    phylo_model_param_block_map[model_keys.CLOCK_MODEL_RATES][:] = np.array([0.001])
    return phylo_model_param_block_map


def create_golden():
    golden_inst = create_instance()
    golden_model = initialize_model_parameters(golden_inst)
    golden_gradients = golden_inst.phylo_gradients()
    golden_likelihoods = golden_inst.log_likelihoods()
    return golden_inst, golden_model, golden_gradients, golden_likelihoods


def compare_gradient_returns(gradients, golden_gradients, test_gradient_keys):
    # test that proper keys are populated
    actual_keys = set(gradients[0].gradient.keys())
    expected_keys = set(test_gradient_keys)
    if actual_keys != expected_keys:
        print("ERROR: Output keys do not match expected keys.")
        print("Expected keys: ", expected_keys)
        print("Actual keys: ", actual_keys)
        return False
    # test that data is correct
    for key in actual_keys:
        actual_data = np.array(gradients[0].gradient[key])
        expected_data = np.array(golden_gradients[0].gradient[key])
        if np.abs(actual_data - expected_data).max() > 0.001:
            print("ERROR: Output data does not match expected data.")
            print("Expected data: ", expected_data)
            print("Actual data: ", actual_data)
            return False
    return True


def test_gradient_include_flags():
    test_passed = True
    _, _, golden_gradients, _ = create_golden()
    inst = create_instance()
    initialize_model_parameters(inst)
    inst.init_phylo_flags()
    inst.set_phylo_defaults(False)
    for flag in include_flags_to_keys:
        inst.set_phylo_flag(flag, True)
        gradients = inst.phylo_gradients([flag], False)
        single_test_passed = compare_gradient_returns(
            gradients, golden_gradients, include_flags_to_keys[flag]
        )
        if single_test_passed == False:
            test_passed = False
        inst.clear_phylo_flags()
    return test_passed


def test_gradient_exclude_flags():
    test_passed = True
    _, _, golden_gradients, golden_likelihoods = create_golden()
    inst = create_instance()
    initialize_model_parameters(inst)

    # check that liklihoods exclude log determinant
    likelihoods = inst.log_likelihoods(
        [(flags.INCLUDE_LOG_DET_JACOBIAN_LIKELIHOOD, False)], True
    )

    logdet = inst.log_det_jacobian_of_height_transform()
    likelihood_with_logdet = np.array(golden_likelihoods)
    likelihood_without_logdet = np.array(likelihoods)
    likelihood_without_logdet_plus_logdet = likelihood_without_logdet + logdet
    max_diff = np.abs(likelihood_with_logdet - likelihood_without_logdet).max()
    if max_diff < 0.001:
        print(
            "ERROR: likelihood_with_logdet == likelihood_without_logdet. max_diff: ",
            max_diff,
        )
        print("EXPECTED_DATA: ", likelihood_with_logdet)
        print("ACTUAL_DATA: ", likelihood_without_logdet)
        test_passed = False
    max_diff = np.abs(
        likelihood_with_logdet - likelihood_without_logdet_plus_logdet
    ).max()
    if max_diff > 0.001:
        print(
            "ERROR: likelihood_with_logdet != likelihood_without_logdet_plus_logdet. max_diff: ",
            max_diff,
        )
        print("EXPECTED_DATA: ", likelihood_with_logdet)
        print("ACTUAL_DATA: ", likelihood_without_logdet + logdet)
        test_passed = False

    # check that gradients exclude log determinant
    gradients = inst.phylo_gradients(
        [(flags.INCLUDE_LOG_DET_JACOBIAN_GRADIENT, False)], True
    )
    key = flags.RATIOS_ROOT_HEIGHT

    grad_without_logdet = np.array(gradients[0].gradient[key])
    grad_with_logdet = np.array(golden_gradients[0].gradient[key])
    logdet = bito.gradient_log_det_jacobian_of_height_transform(
        inst.tree_collection.trees[0]
    )
    grad_without_logdet_plus_logdet = grad_with_logdet + logdet
    max_diff = np.abs(grad_with_logdet - grad_without_logdet).max()
    if max_diff < 0.01:
        print(
            "ERROR: gradients_with_determinant == gradients_without_determinant. max_diff: ",
            max_diff,
        )
        test_passed = False
    max_diff = np.abs(grad_with_logdet - (grad_without_logdet + logdet)).max()
    if max_diff > 0.01:
        print(
            "ERROR: gradients_with_logdet != gradients_without_logdet_plus_logdet. max_diff: ",
            max_diff,
        )
        test_passed = False
    if test_passed == False:
        print("KEY:", key)
        print("GRAD_WITH_LOGDET:\n", grad_with_logdet)
        print("GRAD_WITHOUT_LOGDET:\n", grad_without_logdet)
        print("LOGDET:\n", logdet)
        print("GRAD_WITHOUT_LOGDET_PLUS_LOGDET:\n", grad_without_logdet_plus_logdet)
    return test_passed


def test_gradient_setvalue_flags():
    test_passed = False
    _, _, golden_gradients, _ = create_golden()
    inst = create_instance()
    initialize_model_parameters(inst)
    gradients = inst.phylo_gradients([(flags.SET_GRADIENT_DELTA, 5.0)], True)
    # check that gradients were modified by new value
    for key in gradients[0].gradient.keys():
        actual_data = np.array(gradients[0].gradient[key])
        expected_data = np.array(golden_gradients[0].gradient[key])
        if np.abs(actual_data - expected_data).max() > 0.001:
            test_passed = True
    return test_passed


def test_gradient_pass_flags():
    test_passed = True
    used_flags = [gradient_keys.SUBSTITUTION_MODEL]
    use_defaults = False
    # pass via internal flags
    inst_1 = create_instance()
    initialize_model_parameters(inst_1)
    inst_1.init_phylo_flags()
    inst_1.set_phylo_defaults(use_defaults)
    for flag in used_flags:
        inst_1.set_phylo_flag(flag)
    gradients_1 = inst_1.phylo_gradients()
    # pass via arguments
    inst_2 = create_instance()
    initialize_model_parameters(inst_2)
    gradients_2 = inst_2.phylo_gradients(used_flags, use_defaults)
    # compare results
    # test that proper keys are populated
    internal_keys = set(gradients_1[0].gradient.keys())
    external_keys = set(gradients_2[0].gradient.keys())
    if internal_keys != external_keys:
        print("ERROR: Internal passed keys do not match external passed keys.")
        print("internal_keys: ", internal_keys)
        print("external_keys: ", external_keys)
        return False
    # test that data is correct
    for key in internal_keys:
        data_1 = np.array(gradients_1[0].gradient[key])
        data_2 = np.array(gradients_2[0].gradient[key])
        if np.abs(data_1 - data_2).max() > 0.001:
            print("ERROR: Output data does not match expected data.")
            print("data_1: ", data_1)
            print("data_2: ", data_2)
            return False
    return test_passed


# run tests if called directly
if __name__ == "__main__":
    print("# --- TEST --- #")
    test_1 = test_gradient_include_flags()
    print("# Include Flag Test: ", test_1)
    test_2 = test_gradient_exclude_flags()
    print("# Exclude Flag Test: ", test_2)
    test_3 = test_gradient_setvalue_flags()
    print("# SetValue Flag Test: ", test_3)
    test_4 = test_gradient_pass_flags()
    print("# Internal/External Pass Flags Test: ", test_4)
    test_passed = test_1 and test_2 and test_3 and test_4
    print("# Test Results: ", test_passed)
    print("# --- DEMO --- #")
    gradients_with_flags_demo()
    print("# --- COMPLETE --- #")
