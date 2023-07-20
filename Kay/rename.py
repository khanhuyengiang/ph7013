# Un-comment and run this code in your notebook:
# (Note: No need to import qutip, change path file if needed)
# with open("rename.py") as f:
#    exec(f.read())

from qutip import *

annihilation_op = destroy
creation_op = create
coherent_density_matrix = coherent_dm
expectation_value = expect
fock_density_matrix = fock_dm
identity = qeye
thermal_density_matrix = thermal_dm
sigma_x = sigmax
sigma_z = sigmaz
sigma_y = sigmay

X = sigmax
Y = sigmay
Z = sigmaz
