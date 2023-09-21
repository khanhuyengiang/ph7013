
########## Save Fit Error Variables ############
array A_err[Total]
array b_err[Total]
array c_err[Total]

output(N) = sprintf("A_err[%d] = A_%d__err", N,N)


do for [i=start:Total]{
   eval(output(i))
}




########## Calculate Fidelity ###########

array fidelity[Total]
array fidelity_err[Total]

saveErr(N) = sprintf("p_std = c_%d__err", N)

do for [i=start:Total]{
   cliff_err = (1-c[i])*(1-0.5)
   gate_err = cliff_err / 1.875
   fidelity[i] = (1 - gate_err) * 100000 / 1000

   eval(saveErr(i))
   cliff_err_std = p_std * (1 - 1/2)
   gate_err_std = cliff_err_std / 1.875
   fidelity_err[i] = (gate_err_std * 100000) / 1000
}

########## Save Error ###########


set print "Fit Results.csv"
print "i ',' A ',' A_err ',' b ',' b_err ',' c ',' c_err ',' chi2 ',' Fidelity ',' Fidelity_err ','"

output(N) = sprintf("print %d, ',', A[%d], ',', A_%d__err, ',', b[%d], ',', b_%d__err, ',', c[%d], ',', c_%d__err, ',', chi2[%d], ',', fidelity[%d], ',', fidelity_err[%d]", N,N,N,N,N,N,N,N,N,N)


do for [i=start:Total]{
   eval(output(i))
}

set print

print "chi2"
do for [i=start:Total]{
#   print(chi2[i])
}


set terminal jpeg enhanced font Times 28 size 1600,1200
set output "plot.jpg"

rep

set term windows