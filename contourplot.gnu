set parametric 
set contour base
set view 0,0,1
unset surface
set cntrparam levels 5
splot 'vel.500' using 1:2:3 with line
