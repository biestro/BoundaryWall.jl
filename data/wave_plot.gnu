# set terminal pdfcairo size 12cm,8cm enhanced font "Consolas,8"
# set terminal epslatex # font "Consolas,8"
set terminal pngcairo  font "Consolas,8"
set termopt enhanced

set output "wave_plot.png"

set size ratio 1

wave  = 'wave_parabolic_k=1.4366627996854378.txt'
angle = 'angle_parabolic_k=1.4366627996854378.txt'
tmat  = 'tmat_parabolic_k=1.4366627996854378.txt'
billiard  = 'billiard.txt'


# set palette negative

set multiplot layout 1,2 title "Density and Angle" font ",12"
set view map 


set ylabel rotate by 0


# splot angle matrix with pm3d
# splot wave using ($1):($2):($3) with pm3d notitle
# plot tmat matrix with image
set title "Density"
set xlabel "x"
set ylabel "y"
stats wave using 3
# load '~/colormaps/scientific/grayC.pal'
# set palette negative
set palette defined ( 0 0 0 0, 1 1 1 1 )
set palette negative

set cbtics ('1' STATS_max, '0' 0)

splot wave using 1:2:3 with pm3d notitle,\
      billiard using 1:2:(0) with lines notitle

set title "Angle"
set xlabel "x"
set ylabel "y"
load '~/colormaps/scientific/vikO.pal'
stats angle using 3
set pm3d interpolate 0,0
set cbtics ('Pi' STATS_max, '0' 0, '-Pi' STATS_min)
# set palette model RGB
# set palette model HSV defined ( 0 0 1 1, 1 1 1 1 )

splot angle using 1:2:3 with pm3d notitle,\
      billiard using 1:2:(0) with lines notitle