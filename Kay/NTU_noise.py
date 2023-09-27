# Qutip
from qutip.qip.noise import ControlAmpNoise
from qutip.qip.operations.gates import *
import numpy as np


class FPGA_noise(ControlAmpNoise):
    """
    Random noise in the amplitude of the control pulse. The arguments for
    the random generator need to be given as key word arguments.

    Parameters
    ----------
    dt: float, optional
        The time interval between two random amplitude. The coefficients
        of the noise are the same within this time range.
    rand_gen: numpy.random, optional
        A random generator in numpy.random, it has to take a ``size``
        parameter as the size of random numbers in the output array.
    indices: list of int, optional
        The indices of target pulse in the list of pulses.
    **kwargs:
        Key word arguments for the random number generator.

    Attributes
    ----------
    dt: float, optional
        The time interval between two random amplitude. The coefficients
        of the noise are the same within this time range.
    rand_gen: numpy.random, optional
        A random generator in numpy.random, it has to take a ``size``
        parameter.
    indices: list of int
        The indices of target pulse in the list of pulses.
    **kwargs:
        Key word arguments for the random number generator.

    Examples
    --------
    >>> gaussnoise = FPGA_noise( \
            dt=0.1, rand_gen=np.random.normal, loc=mean, scale=std) \
            # doctest: +SKIP
    """
    def __init__(self, dt, rand_gen, amplitude, indices=None, **kwargs):
        super(FPGA_noise, self).__init__(coeff=None, tlist=None)
        self.rand_gen = rand_gen
        self.kwargs = kwargs
        if "size" in kwargs:
            raise ValueError("size is preditermined inside the noise object.")
        self.dt = dt
        self.indices = indices
        self.amplitude = amplitude
    
    def get_noisy_dynamics(
            self, dims=None, pulses=None, systematic_noise=None):
        if pulses is None:
            pulses = []
        if self.indices is None:
            indices = range(len(pulses))
        else:
            indices = self.indices
        t_max = -np.inf
        t_min = np.inf
        for pulse in pulses:
            t_max = max(max(pulse.tlist), t_max)
            t_min = min(min(pulse.tlist), t_min)
        # create new tlist and random coeff
        num_rand = int(np.floor((t_max - t_min) / self.dt)) + 1
        tlist = (np.arange(0, self.dt*num_rand, self.dt)[:num_rand] + t_min)
        # [:num_rand] for round of error like 0.2*6=1.2000000000002

        for i in indices:
            pulse = pulses[i]
            coeff = self.amplitude * self.rand_gen(**self.kwargs, size=num_rand)
            pulses[i].add_coherent_noise(
                pulse.qobj, pulse.targets, tlist, coeff)
        return pulses, systematic_noise
