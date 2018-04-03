set xrange [1:50]
set term eps
set output "velProfiles_x15_basic_reg_forcing.eps"
set grid
plot 'velprofile_fshift.12000' using 2:3 with linespoints title 'FShift_noreg', 'velprofile_with_fshift_reg.10250' using 2:3 title 'FShift_withreg'

f(x) = a*x*x + b*x + c
fit f(x) 'velprofile_fshift.12000' using 2:3 via a,b,c

g(x) = d*x*x + e*x + f 
fit f(x) 'velprofile_with_fshift_reg.10250' using 2:3 via d,e,f
