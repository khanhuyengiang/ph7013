reset

set term win enhanced font 'Arial, 20'

set samples 1000

#set key autotitle columnhead
set key top right
set key font ",15"
set key noreverse 
set key width -2
set datafile sep ','
set autoscale

set xlabel "DAC Value"
set ylabel "SNR \n (Signal Peak Power / Average Noise Power)" offset -2,0
set lmargin 13

set y2label "Fidelity (%)"

#set format y "%0.2t x 10^{%L}"  
set format y "10^{%L}" 


set logscale y


set ytics nomirror

set grid


gain = 15.5
output = -14.77
total = gain + output
sig_Power =  10**(total/10)/1e3 # convert from mW to W

set y2tics

##### Set initial range
set yr[1e5:5e7] 
set y2r[98:100]
set xr[-1000:9000]




f10 = "./PowerCalibration/20230406_dataEven_bw_10.csv"
f50 = "./PowerCalibration/20230406_dataEven_bw_50.csv"
f100 = "./PowerCalibration/20230406_dataEven_bw_100.csv"


p f10 using 1:(sig_Power/(10**($2/10)/1e3)) w p axis x1y1 t "SNR @ Noise BW: 10 MHz" lc 2
rep f50 using 1:(sig_Power/(10**($2/10)/1e3)) w p axis x1y1 t "SNR @ Noise BW: 50 MHz" lc 3
rep f100 using 1:(sig_Power/(10**($2/10)/1e3)) w p axis x1y1 t "SNR @ Noise BW: 100 MHz" lc 4

rep 'Fit Results.csv' using 1:9:10 w yerrorbar axis x1y2 t "Fidelity" lc 1
rep 'Fit Results.csv' using 1:9:10 w l axis x1y2 lc 1 lw 1 t ""


########## Fittings


a = 1.45
b = -5.2
c = -1.5
d = 8900
d(x) =  a*1e5*exp(c*1e-5*(x-d)-(b*1e-8)*(x-d)**2)  

set xr [0:4000]
fit d(x) f10 using 1:(sig_Power/(10**($2/10)/1e3)) via b

rep d(x) w l lw 2 t ''# 'SNR Eye Guide' lc 6
set label 3 at   0,  1.00453e+06 "SNR Guide \U+2190 " textcolor lt 6

target_fid = 99.99

set xr [0:2500]
fid1(x) = m1*x + c1
m1 = -0.000123553
c1 = 99.8598  
fit fid1(x) 'Fit Results.csv' using 1:9:10 yerror via m1,c1
rep fid1(x) w l axis x1y2 t "" #"Partial Fidelity Fit" 
set label 1 at  3398.44,  1.21606e+07 "{\U+2192 Partial Range Fidelity fit }" textcolor linetype 7
dac1 = (target_fid -c1)/m1
snr1 = a*1e5*exp(c*1e-5*(dac1-d)-(b*1e-8)*(dac1-d)**2)  
noise1 = sig_Power / snr1 # Watts
noise1_dBm = 10 * log10(noise1) + 30


set xr[-1000:9000]
fid2(x) = m2*x + c2
m2 = -0.000204255
c2 = 99.9509
fit fid2(x) 'Fit Results.csv' using 1:9:10 yerror via m2,c2
rep fid2(x) w l axis x1y2 t "" #"Full Fidelity Fit"
set label 2 at    3740.66,  4.49126e+06  "{\U+2192 Full Range Fidelity Fit}"
dac2 = (target_fid -c2)/m2
snr2 =  a*1e5*exp(c*1e-5*(dac2-d)-(b*1e-8)*(dac2-d)**2)  
noise2 = sig_Power / snr2
noise2_dBm = 10 * log10(noise2) + 30

print ""
print "########################"
print "Target Fidelity =", target_fid
print "Fit based on Red Line Fidelity:"
print "Ideal DAC for partial fit =", dac1
print "Ideal SNR for respective DAC = ", snr1/1e6, "x 10^6"
print "Required Noise Level =", noise1_dBm, "dBm"

print ""
print "Fit based on Brown Line Fidelity:"
print "Ideal DAC for full fit =", dac2
print "Ideal SNR for respective DAC = ", snr2/1e6, "x 10^6"
print "Required Noise Level =", noise2_dBm, "dBm"
print "########################"



rep

set terminal jpeg enhanced font Arial 20 size 1000,900
set encoding utf8

set output "plot.jpg"

rep

set output

set terminal win

rep

