set grid
set xrange [0:0.3]
set yrange [-0.2:1]

cs2 = 0.333
Uo = 1.0
Kn = 0.2
L = 128
tau = (Kn*L)/(2*cs2)
kvisc = cs2*(tau - 0.5)
t0 = (L*L)/kvisc
set term eps
set output "nondimshearprof.eps"

plot 'velAmpreg.dat' using ($1/t0):($2/Uo) with linespoints, 'velAmp.dat' using ($1/t0):($2/Uo), 'paperdata.dat' using 1:2
