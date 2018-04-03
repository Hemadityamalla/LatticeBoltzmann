set xrange [1:50]
set grid
set datafile separator " "
set term eps
set output "velProfiles_x15_basic_noreg_forcing.eps"
plot 'velProfiles_x15_noReg.dat' index 0 using 2:3 with points title 'NoReg', \
'velProfiles_x15_noReg.dat' index 1 using 2:3 dashtype 2 with lines title 'Reg'
