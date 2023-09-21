set samples 1000

set term windows

set datafile sep ','

file = 'data.csv'

set key autotitle columnhead
set fit errorvariables  

start = 1
Full = 20
Total = Full

set lt cycle Total



########### Plot Datapoints ###########

set pointsize 2

plot for [i=start:Total] file using 1:2*i:2*i+1 with yerrorbar 


########### Create parameters and initialize them ###########

array A[Full] 
array b[Full] 
array c[Full] 
array chi2[Full]

initVar(N) = sprintf("A[%d] = %f; b[%d] = %f; c[%d] = %f", N, 0.4, N, -0.4, N, 1)
#initVar(N) = sprintf("A[%d] = %f; b[%d] = %f; c[%d] = %f", N, 1, N, 2, N, 0.38)


do for[i=start:Total]{
   eval(initVar(i))
}


load "2 var.plt"



########### The equation string ###########
#equ(N) = sprintf("y%d(x) = A[%d]*exp(-x/b[%d])+c[%d]",N,N,N,N)
equ(N) = sprintf("y%d(x) = A[%d] - b[%d]*c[%d]**x",N,N,N,N)


do for[i=start:Total]{
#   print(equ(i))
   eval(equ(i))
}


########### The fitting string ###########
fitEq(N) = sprintf("fit y%d(x) file using 1:2*%d:2*%d+1 yerror via A[%d], b[%d], c[%d]", N,N,N,N,N,N)
#fitEq(N) = sprintf("fit y%d(x) file using 1:2*%d:2*%d+1 yerror via b[%d], c[%d]", N,N,N,N,N)


do for [i=start:Total]{
#   eval(fitEq(i))
   chi2[i] = (FIT_STDFIT*FIT_STDFIT)
}



########### The plotting string ###########
#draw1(N) = sprintf("p y%d(x) w l t 'Fit:%d'", N, N-1)
#eval(draw1(1))
#start = start+1

draw(N) = sprintf("rep y%d(x) w l t 'Fit:%d'", N, N-1)

do for [i=start:Total]{
#   print(draw(i))
   eval(draw(i))
}

#set xr [0:1000];rep
set xr [0:*]
set yr [0.3:1]
rep

#load 'save.plt'
