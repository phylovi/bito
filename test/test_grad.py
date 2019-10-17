import libsbn
import numpy as np

def test_grad():
    inst = libsbn.instance('charlie')
    inst.read_newick_file('data/hello.nwk')
    inst.read_fasta_file('data/hello.fasta')
    inst.make_beagle_instances(1)

    evec, ivec, eval, freqs, q_differential = [np.array(x, copy=False) for x in [inst.evec, inst.ivec, inst.eval, inst.freqs, inst.q_differential]]

    evec[:] = np.array([1.0, 1.9607843160629272, 0.0, 0.47058823704719543, 1.0, -2.040816307067871, 0.4693877696990967, 0.0, 1.0, 1.9607843160629272, 0.0, -0.529411792755127, 1.0, -2.040816307067871, -0.5306122303009033, 0.0])
    ivec[:] = np.array([0.27, 0.26, 0.24, 0.23, 0.1323, -0.1326, 0.1176, -0.1173, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0, 0.0])
    eval[:] = np.array([0.0, -1.0018032789230347, -1.4926868677139282, -1.5127228498458862])
    freqs[:] = np.array([0.27, 0.26, 0.24, 0.23])

    inst.set_beagle_subst_models()
    np.testing.assert_allclose(np.array(inst.log_likelihoods()), -84.5624981081271)

    q_differential[:] = np.array([0.0021637824504385783, -0.06502588257063284, 0.12038499624036944, -0.057522896120175214, -0.06752687732963392, 0.012181812693820238, -0.06002389095967459, 0.11536895559548827, 0.13543312424568438, -0.06502588257063284, -0.012884345554876386, -0.057522896120175214, -0.06752687732963392, 0.1304170765614889, -0.06002389095967459, -0.0028663082721803734])
    np.testing.assert_allclose(np.array(inst.parameter_gradients()), -0.019635463018150334)

if __name__ == '__main__':
    test_grad()
    
