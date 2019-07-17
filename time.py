import numpy as np
import sbn
import timeit

v1_c = sbn.make_vector()
v2_c = sbn.make_vector()
v3_c = sbn.make_vector()
a1_c = np.array(v1_c, copy=False)
a2_c = np.array(v2_c, copy=False)
a3_c = np.array(v3_c, copy=False)

print("pybind11 version: ", timeit.timeit("a3_c = a1_c + a2_c", globals=globals()))

a1 = np.random.uniform(len(a1_c))
a2 = np.random.uniform(len(a1_c))
a3 = np.random.uniform(len(a1_c))

print("numpy version: ", timeit.timeit("a3 = a1 + a2", globals=globals()))

