import numpy as np
from qutip_qip.operations import Gate
from qutip.qip.operations.gates import *
import itertools
from qutip import Qobj

__all__ = ['gate_set_generator', 'matrix_list', 'add_inverse_gates']

def gates_set_generator(x):
    """ Generate a set of RX and RY gates with argument value as a multiple of pi.
    """
    return [
    Gate("RX", 0, arg_value= x * np.pi), # X Pulse
    Gate("RY", 0, arg_value= x * np.pi), # Y Pulse
    Gate("RX", 0, arg_value= x * np.pi / 2), # X Half Pulse
    Gate("RY", 0, arg_value= x * np.pi / 2), # Y Half Pulse
    Gate("RX", 0, arg_value=- x * np.pi / 2), # X Minus Half Pulse
    Gate("RY", 0, arg_value=- x * np.pi / 2), # Y Minus Half Pulse
    ]

matrix_list = [
    rx(np.pi), # X
    ry(np.pi), # Y,
    rx(np.pi / 2), # X/2
    ry(np.pi / 2), # Y/2
    rx(-np.pi / 2), # -X/2
    ry(-np.pi / 2), # -Y/2
    rx(0)
]
        
def add_inverse_gates(clifford, init_state, matrix_list = matrix_list, gates_set = gates_set_generator(1), circuit = None):
    """Add 0,1 or 2 gates that inverse the given gate(s)

        Args:
            clifford (qutip.Qobj): A qutip Qobj result from multiplying 
            the sequence of Clifford apply to the qubit.
            init_state (qutip.Qobj): The initial state.
            matrix_list (list of qutip.Qobj): A list of quantum gates in matrix form.
            gates_set (list of qutip_qip.circuit.Gate): A list of qutip Gate objects 
            that corresponds to the quantum gates given by matrix_list.
            circuit (qutip_qip.circuit.QubitCircuit): The circuit that the gates applied on.
        """
    def _inverse_search(state_before_inverse, init_state, matrix_list = matrix_list):
        """Find 2 gates from a list of gates that when apply to a given state give the initial state.

        Args:
            state_before_inverse (qutip.Qobj): The quantum state that need to be inversed.
            init_state (qutip.Qobj): The initial state.
            matrix_list (list of qutip.Qobj): A list of quantum gates in matrix form.

        Returns:
            index_list (tuple): The index (indicies) in the matrix list 
            of the two gates that inverse the state.
        """
        # Index of 2 gates that inverse the clifford sequence
        index_list = [i for i in itertools.product(range(len(matrix_list)), range(len(matrix_list)))]
        # Product of said 2 gates
        inverse_list = [i[0]*i[1] for i in itertools.product(matrix_list, matrix_list)]
        
        for i in range(len(inverse_list)):
            # Find final state after applying the 2 inverse gates and "normalize"
            final_state = inverse_list[i] * state_before_inverse
            if np.round(final_state[0][0][0],2) == 0:
                final_state = final_state/final_state[1][0][0]
            else:
                final_state = final_state/final_state[0][0][0]
            # Compare to ground state
            if np.allclose(final_state,init_state):
                return index_list[i]
    
    x = _inverse_search(clifford*init_state, init_state, matrix_list)
    if x == None:
        raise RuntimeError("Could not find an inverse Clifford")
    elif circuit == None:
        raise TypeError('NoneType object')
    elif not x[0] == 6 and not x[1] == 6:
        circuit.add_gate(gates_set[x[1]])
        circuit.add_gate(gates_set[x[0]])
    elif x[0] == 6 and not x[1] == 6:
        circuit.add_gate(gates_set[x[1]])
    elif x[1] == 6 and not x[0] == 6:
        circuit.add_gate(gates_set[x[0]])