# File to generate radiative power values
# fmt: off
"""
@author: phineas
"""

TestPos = 0

#External modules
import numpy as np
import math as m

#Internal modules
from cea import *
import InputValues as IV
from RayMarch import getRay
from concurrent.futures import ProcessPoolExecutor
from ComputeRays import CompRay

n_cores = 18  # how many processes you want

inputs = np.array([(0, 0, 0)])
M = np.linspace(-90, 90, 181)
N = np.linspace(0, 90, 91)
M[0] = -89.99
M[-1] = 89.99
N[-1] = 89.99
#print(M)

i = 0
while i < len(M):
    j = 0
    while j < len(N):
        inputs = np.append(inputs, np.array([(TestPos, M[i]*(m.pi/180), N[j]*(m.pi/180))]), axis = 0)
        j+=1
    i+=1
inputs = np.delete(inputs, 0, 0)

def call_CompRay(args):
    return CompRay(*args)

if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        results = list(executor.map(call_CompRay, inputs))
    results = np.array(results, dtype=np.float64)
    name = str(TestPos) + "PowerSter"
    np.save(name, results)
'''
if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        futures = [executor.submit(call_CompRay, x) for x in inputs]
        finished = 0
        total = len(futures)
        for future in as_completed(futures):
            result = future.result()
            finished += 1
            print(f"{finished}/{total} finished")
    results = np.array(results, dtype=np.float64)
    name = str(TestPos) + "PowerSter"
    np.save(name, results)
'''