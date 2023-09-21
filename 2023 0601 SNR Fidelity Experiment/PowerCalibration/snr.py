import numpy as np
import matplotlib.pyplot as plt
import csv


def dBm2Watts(dBm):
  milliwatts = 10**(dBm/10)
  return milliwatts

gain = 15.5
sig_Power =  -14.77 + gain

sig_Power_w = dBm2Watts(sig_Power)

with open("20230406_dataEven_bw_100.csv", "r") as file:
  reader = csv.reader(file, delimiter=',')
  raw = list(reader)
  noise = np.array(raw[1:]).astype('float')
  noise = np.transpose(noise)


dac = noise[0]
mean = noise[1]
std = noise[2]

snr = dBm2Watts(sig_Power)/dBm2Watts(mean)


dac_fid = [0,431,862,1293,1724,2155,2586,3017,3448,3879,4311,4742,5173,5604,6035,6466,6897,7328,7759,8191
]

with open("fidelity.csv", "r") as file:
  reader = csv.reader(file, delimiter=',')
  raw = list(reader)
  data = np.array(raw[1:]).astype('float')
  data = np.transpose(data)
  fid = data[8]
  fid_err = data[9]

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(dac, snr, 'g')
ax2.plot(dac_fid, fid, 'b')

ax1.set_xlabel('DAC value')
ax1.set_ylabel('SNR', color='g')
ax1.ticklabel_format(style='plain')
ax1.set_yscale("log")
ax1.set_ylim(1e5, 1e7)
#ax2.set_ylim([10e5, 10e7])
ax2.set_ylabel('Fidelity', color='b')
ax1.grid()



plt.show()

