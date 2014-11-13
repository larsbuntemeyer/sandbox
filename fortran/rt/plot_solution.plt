set terminal postscript eps enhanced color solid 18
set output "solution.eps"
#set xlabel "Radius {/Symbol x} [r_{BE}]"
#set ylabel "Dichte D [{/Symbol r_{0}}]"
set log x
set key
set xrange [0.01:100000]
set yrange [0.000001:1.2]
#set style line 2 lt 3 lw 2 pt 3 ps 0.5
#set style function linespoints
#set ytics 0.01,0.1,1.0
plot "plot.out" using 1:2 title "Plane Parallel Atmosphere" lt 1 pt 1
