#!/usr/bin/gnuplot

set object 1 polygon from 0,.5 to \
    0.25,0.75 to \
    .5,.5 to \
    0,.5 \
    fs transparent solid .25 fc rgb '#0060ad'
set object 2 polygon from 0,.5 to \
    0.25,0.75 to \
    0,1 to 0,.5 \
    fs transparent solid .25 fc rgb '#0060ad' 

set xlabel "x" offset 2
set ylabel "y" offset 2
set xrange [-.5:1.0]
set yrange [-.25:1.25]
set xtics 0.5,.5,.5 axis
set ytics 0.5,.5,1 axis
set style line 11 lc rgb 'black' lt 1
set style line 12 lc rgb '#808080' lt 0 lw 1
set tics nomirror
set xzeroaxis ls 11
set yzeroaxis ls 11
set arrow from 1,0 to 1.05,0 size screen 0.025,15,60 filled ls 11
set arrow from 0,1.25 to 0,1.3 size screen 0.025,15,60 filled ls 11
set noborder
set border 0 back ls 11

set label '$E_1$' center at .30,.8
set label '$E_2$' center at .05,.45
set label '$E_3$' center at .05,1.05
set label '$E_4$' center at .55,.55
set label '$T_1$' center at .0975,.75
set label '$T_2$' center at .25,.6125


# set label '$\mathcal{S}_1$' center at -.05,.75
# set label '$\mathcal{S}_1$' center at .25,.44
# set label '$\mathcal{S}_2$' center at .175,.9
set terminal tikz
set output 'triangles.tikz'
set style line 1 lc rgb 'red' pt 7 lt 1 lw 2 ps 0
set style line 2 lc rgb 'green' pt 7 lt 1 lw 2 ps 0
# set style line 3 fs transparent solid .75 lc rgb 'blue' pt 9 ps 7
set style line 3 lc rgb 'blue' pt 9 ps 7 
set nokey

plot '-' with linespoints lc rgb '#0060ad' pt 7 lt 1 lw 2 ps 1.5, \
     '-' with linespoints lc rgb '#0060ad' pt 7 lt 1 lw 2 ps 1.5, \
     '-' with linespoints lc rgb '#0060ad' pt 7 lt 1 lw 2 ps 1.5, \
     '-' with linespoints lc rgb '#0060ad' pt 7 lt 1 lw 2 ps 1.5
0 0.5
0 1
e
0 0.5
.5 .5
e
0 1
.5 .5
e
0 0.5
.25 .75
e