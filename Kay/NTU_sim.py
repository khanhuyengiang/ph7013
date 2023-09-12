"""NTU_processor and NTU_compiler class, for simulating the qubit. 
Function NTU_single_sim for simulating a single simulation, 
and function NTU_sim_test_run for simulating (and average over) several simulation.
"""

__all__ = ['NTU_processor','NTU_compiler','NTU_single_sim', 'NTU_sim_test_run']

# Qutip
from qutip import (sigmax, sigmay, tensor, basis)
from qutip.metrics import fidelity
from qutip_qip.circuit import QubitCircuit
from qutip_qip.compiler import GateCompiler, Instruction
from qutip_qip.device import ModelProcessor
from qutip.qip.noise import RandomNoise
from qutip.qip.operations.gates import *

import numpy as np
import functools # for reduce
from scipy.signal import argrelextrema 
# Import function to generate a gate set and add inverse gates 
from inverse_search import gates_set_generator, matrix_list, add_inverse_gates


class NTU_processor(ModelProcessor):
    """
    Custom processor built using ModelProcessor as the base class.
    This custom processor will inherit all the methods of the base class
    such as setting up of the T1 and T2 decoherence rates in the simulations.

    Args:
        num_qubits (int): Number of qubits in the processor.
        t1, t2 (float or list): The T1 and T2 decoherence rates for the
    """

    def __init__(self, num_qubits, t1=None, t2=None):
        # call the parent class initializer
        super(NTU_processor, self).__init__(num_qubits, t1=t1, t2=t2)  
        # The control pulse is discrete or continous.
        self.pulse_mode = "discrete"
        # The dimension of each controllable quantum system
        self.model.dims = [2] * num_qubits
        self.num_qubits = num_qubits
        self.set_up_ops()  # set up the available Hamiltonians
        self.native_gates = ["RX", "RY"]

    def set_up_ops(self):
        """
        Sets up the control operators.
        """
        # sigmax pulse on m-th qubit with the corresponding pulse
        for m in range(self.num_qubits):
            self.add_control(1/2*sigmax(), m, label="sx" + str(m))
        # sy
        for m in range(self.num_qubits):
            self.add_control(1/2*sigmay(), m, label="sy" + str(m))
    
    # To remove errors arise from RandomNoise not being recognized
    # as a type Noise, may cause trouble in the future?
    def add_noise(self, noise):
        """
        Add a noise object to the processor

        Parameters
        ----------
        noise: :class:`.Noise`
            The noise object defined outside the processor
        """
        self.noise.append(noise)


class NTU_compiler(GateCompiler):
    """
    Custom compiler for generating pulses from gates using
    the base class GateCompiler.

    Args:
        num_qubits (int): The number of qubits in the processor
        params (dict): A dictionary of parameters for gate pulses
                       such as the pulse amplitude.
    """

    def __init__(self, num_qubits, params):
        super().__init__(num_qubits, params=params)
        self.params = params
        self.gate_compiler = {
            "RX": self.pulse_discretization_compiler,
            "RY": self.pulse_discretization_compiler,
        }

    def generate_pulse(self, gate, tlist, coeff, phase=0.0):
        """Generates the pulses.

        Args:
            gate (qutip_qip.circuit.Gate): A qutip Gate object.
            tlist (array): A list of times for the evolution.
            coeff (array): An array of coefficients for the gate pulses
            phase (float): The value of the phase for the gate.

        Returns:
            Instruction (qutip_qip.compiler.instruction.Instruction):
            An instruction to implement a gate containing the control pulses.
        """
        pulse_info = [
            # (control label, coeff)
            ("sx" + str(gate.targets[0]), phase[0] * coeff),
            ("sy" + str(gate.targets[0]), phase[1] * coeff),
        ]
        return [Instruction(gate, tlist=tlist, pulse_info=pulse_info)]

    def pulse_discretization_compiler(self, gate, args):
        """Compiles single-qubit gates to pulses.

        Args:
            gate (qutip_qip.circuit.Gate): A qutip Gate object.

        Returns:
            Instruction (qutip_qip.compiler.instruction.Instruction):
            An instruction to implement a gate containing the control pulses.
        """
        # gate.arg_value is the rotation angle
        if gate.name == "RX":
            phiNaught = 0
        elif gate.name == "RY":
            phiNaught = np.pi/2
        
        VNaught = self.params["VNaught"]
        VStd = self.params["VStd"]
        phaseStd = self.params["phaseStd"]
        omega = self.params["omega"]
        aNaught = self.params["aNaught"]
        detuningStd = self.params["detuningStd"]

        V = VNaught + np.random.normal(scale=VStd)
        phi = phiNaught + np.random.normal(scale=phaseStd)
        I = np.cos(phi)
        Q = np.sin(phi)

        _step_list = np.linspace(0,1,11) # a time step list so that the noise work
        coupling_time_series = np.abs(gate.arg_value) / (V*omega) * _step_list
        s = aNaught - (1 - aNaught) * np.cos(2 * np.pi * _step_list[:-1])
        #FPGA_voltage = V * omega * np.sign(gate.arg_value) * s
        FPGA_voltage = V * omega * np.sign(gate.arg_value) * (_step_list[:-1]*0 + 1)
        #dwt = np.random.normal(scale=0.1) * coupling_time_series[:-1]
        dwt = np.random.normal(scale=detuningStd) * coupling_time_series[:-1]
        phase = [- I * np.cos(dwt) + Q * np.sin(dwt),
                 I * np.sin (dwt) - Q *np.cos(dwt)]
        if gate.name == "RX":
            return self.generate_pulse(gate, tlist = coupling_time_series, coeff = FPGA_voltage, phase=phase)
        elif gate.name == "RY":
            return self.generate_pulse(gate, tlist = coupling_time_series, coeff = FPGA_voltage, phase=phase)


def NTU_single_sim(num_qubits, num_gates, gates_set, param_dict,
                   t1 = None, t2 = None,
                   add_FPGA_noise = True,):
    """
    A single simulation, with num_gates representing the number of rotations.
    
    Args:
        num_gates (int): The number of random gates to add in the simulation.
        t1, t2 (float): Decoherence time of the qubits.
        add_FPGA_noise (bool): Whether to add in gaussian FPGA noise to the simulation.
        param_dict (dictionary): Dictionary of the following parameters:
                                {"VNaught": VNaught, "VStd": VStd, "phaseStd":phaseStd,
                                "omega": omega, "aNaught": aNaught, "detuningStd": detuningStd,
                                "pulse_amplitude": 20e6}

    Returns:
        final_fidelity (float):
            Fidelity of the result state (obtained from 
            mesolve solver method) and the initial state.
    """

    # Finding the normalizing coefficient for the gate (ie finding x so that cos(x*2np.pi) = 1)

    # The actual simulation
    myprocessor = NTU_processor(num_qubits, t1 = t1, t2 = t2)
    myprocessor.native_gates = None  # Remove the native gates
    mycompiler = NTU_compiler(num_qubits, param_dict)

    # Ground state for n qubits
    init_state = functools.reduce(lambda a, b: tensor(a,b), [basis(2, 0)] * num_qubits)

    # Define a random circuit.
    circuit = QubitCircuit(num_qubits)
    clifford = rx(0)
    for ind in np.random.randint(0, 6, num_gates):
        circuit.add_gate(gates_set[ind])
        clifford = matrix_list[ind] * clifford

    # Finding inverse Clifford for the random sequence of gate
    add_inverse_gates(clifford, init_state, circuit = circuit, gates_set = gates_set)

    # Simulate the circuit.
    myprocessor.load_circuit(circuit, compiler=mycompiler)
    
    # FPGA gaussian noise
    if add_FPGA_noise == True:
        FPGA_noise = RandomNoise(dt=1e-9, indices = [0,1], rand_gen=np.random.normal, loc=0.00, scale=0.3)
        myprocessor.add_noise(FPGA_noise)
    
    # Compute results of the run using a solver of choice
    result = myprocessor.run_state(init_state, solver="mesolve")
    # Measured fidelity at the end
    final_fidelity = fidelity(result.states[0],result.states[-1])
    return final_fidelity


def NTU_sim_test_run(num_qubits: int, num_gates_list: list, num_samples: int, param_dict: dict,
                     t1 = None, t2 = None, add_FPGA_noise = True,
                    ):
    """
    Find the Clifford gate set correspond to the Hamiltonian, 
    then run a sample test run of several simulation, with 
    num_gates_list representing the list of number of rotations that looped through.
    
    Args:
        num_qubits (int): The number of qubits.
        num_gates_list (list): The number of random gates to add in the simulation.
        num_samples (int):
            Number of times to run the simulation and take average
            for a particular number of gates.
        t1, t2 (float): Decoherence time of the qubits.
        add_FPGA_noise (bool): Whether to add in gaussian FPGA noise to the simulation.
        param_dict (dictionary): 
            Dictionary of the following parameters:
            {"VNaught": VNaught, "VStd": VStd, "phaseStd":phaseStd,
            "omega": omega, "aNaught": aNaught, "detuningStd": detuningStd,
            "pulse_amplitude": 20e6}

    Returns:
        final_fidelity (float):
            Fidelity of the result state (obtained from 
            mesolve solver method) and the initial state.
        fidelity_average (list): List of average fidelity correspond to num_gates_list.
        fidelity_error (list): List of errors of the corresponding average fidelity.
    """

    fidelity_list = []
    index_list = []

    for x in np.linspace(0.01,6,100):
        myprocessor = NTU_processor(num_qubits)
        myprocessor.native_gates = None  # Remove the native gates
        mycompiler = NTU_compiler(num_qubits,param_dict)

        # Ground state for n qubits
        init_state = functools.reduce(lambda a, b: tensor(a,b), [basis(2, 0)] * num_qubits)

        # Define a random circuit.
        circuit = QubitCircuit(num_qubits)
        circuit.add_gate("RY", targets=0, arg_value=x*np.pi)

        # Simulate the circuit.
        myprocessor.load_circuit(circuit, compiler=mycompiler)
        
        # Compute results of the run using a solver of choice
        result = myprocessor.run_state(init_state, solver="mesolve")
        # Measured fidelity at the end
        fidelity_list.append(fidelity(result.states[0],result.states[-1]))
        index_list.append(x)
        
    maximum_array = argrelextrema(np.asarray(fidelity_list), np.greater,order = 10)
    first_max = maximum_array[0][0]
    pulse_coeff = index_list[first_max]/2

    if np.round(pulse_coeff,2) == 1:
        gates_set = gates_set_generator(1)
    else:
        gates_set = gates_set_generator(pulse_coeff)
    
    fidelity_average = []
    fidelity_error = []
    for num_gates in num_gates_list:
        fidelity_list = [NTU_single_sim(
            num_qubits, num_gates, gates_set = gates_set,
            t1 = t1, t2 = t2, add_FPGA_noise = add_FPGA_noise,
            param_dict = param_dict,
            ) for i in range(num_samples)]
        fidelity_average.append(np.mean(fidelity_list))
        fidelity_error.append(np.std(fidelity_list) / np.sqrt(num_samples))
    
    return fidelity_average, fidelity_error
