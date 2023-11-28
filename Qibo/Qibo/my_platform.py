# my_platform.py

import pathlib

from qibolab.channels import Channel, ChannelMap
from qibolab.instruments.rfsoc import RFSoC
from qibolab.instruments.rohde_schwarz import SGS100A as LocalOscillator
from qibolab.platform import Platform
from qibolab.serialize import load_qubits, load_runcard, load_settings

NAME = "my_platform"  # name of the platform
ADDRESS = "192.168.0.1"  # ip adress of the controller
PORT = 6000  # port of the controller

# path to runcard file with calibration parameter
RUNCARD = pathlib.Path.cwd() / "my_platform.yml"


def create(runcard_path=RUNCARD):
    # Instantiate controller instruments
    controller = RFSoC(NAME, ADDRESS, PORT)

    # Create channel objects and port assignment
    channels = ChannelMap()
    channels |= Channel("readout", port=controller[1])
    channels |= Channel("feedback", port=controller[0])
    channels |= Channel("drive", port=controller[0])

    # create qubit objects
    runcard = load_runcard(runcard_path)
    qubits, pairs = load_qubits(runcard)
    # assign channels to qubits
    qubits[0].readout = channels["L3-22_ro"]
    qubits[0].feedback = channels["L1-2-RO"]
    qubits[0].drive = channels["L3-22_qd"]

    instruments = {controller.name: controller}
    settings = load_settings(runcard)
    return Platform(NAME, qubits, pairs, instruments, settings, resonator_type="3D")