set xrange [1:50]
set term eps
set output "velProfiles_x15_basic_noreg_forcing.eps"
set grid

f = 0.00005
H = 50
cs2 = 0.3333
tau = 0.6666
kvisc = cs2*(tau - 0.5)
Umax = (f*H*H)/(8*kvisc)

set object circle at first 25,Umax radius char 0.5 fillstyle empty border lc rgb '#aa1100' lw 2

plot 'velprofile_fshift.22500' using 2:3 with linespoints title 'Force Shift', 'velprofile_federico.22500' using 2:3 title 'Federico', 

f(x) = a*x*x + b*x + c
fit f(x) 'velprofile_fshift.22500' using 2:3 via a,b,c


