L = 32.0
u0 = 4.3416e-5

set grid
set xrange [0:1.0]
set yrange [0.0:2.0]

set term eps
set output "VelocityProf_ForcingVariation.eps"

set lmargin 10
set ylabel 'u / U_0' offset 1,0,0 rotate by 0
set xlabel 'y/L'
set xtics 0,0.1,1
#set ytics 0,0.25,2
set title 'Velocity Profiles at x=10 for Kn = 0.2, U_0 = 4.3416e-5, L = 32'

plot 'basic.dat' using ($2/L):($3/u0) with linespoints lt 1 ps 0.5 lw 2 title 'Basic Forcing',\
	'fshift.dat' using ($2/L):($3/u0) with linespoints lt 2 ps 0.5 lw 2 title 'Force Shift',\
	'guo.dat' using ($2/L):($3/u0) with linespoints lt 3 ps 0.5 lw 2 title 'Guo Forcing'

