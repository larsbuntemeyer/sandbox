set termoption dashed
set style line 1 lt 1 lw 2 lc 0
set mxtics 5
set mytics 5
m_sol = 1.989e33 
set xrange [-0.55:0.55]
unset key

set style line 1 lt 1 lw 1.0 lc 3 
set style line 2 lt 2 lw 1.0

set size 1,1
set terminal postscript eps enhanced 18 size 4.0,3.0
set lmargin 8.0
set rmargin 2.5
set tmargin 3.5
set bmargin 3.5

set xlabel "x [cm]"
source_label = "log(S) [erg s^{-1} cm^{-2} sr^{-1}]"
error_label = "log(abs(S-S_d)/S_d)"

nx = 64
set key font ",15"
fi = gprintf("%04.0f",nx)
set output "diffusion_1D_".fi.".eps"
set title "nx = ".gprintf("%02.0f",nx)
set yrange [-2:6]
set format y "%.0f"
set ylabel source_label 
plot "diffusion_1D_".fi.".dat" using 1:(log10($4)) title "1D solution", \
     "diffusion_1D_".fi.".dat" using 1:(log10($3)) title "diffusion solution" with lines ls 1, \
     "diffusion_1D_".fi.".dat" using 1:(log10($2)) title "" with lines ls 2
set yrange [-5:1]
set output "diffusion_1D_err_".fi.".eps"
set ylabel error_label
unset key
plot "diffusion_1D_".fi.".dat" using 1:(log10(abs($4-$3)/$3)) with lines

nx = 32
set key font ",15"
fi = gprintf("%04.0f",nx)
set output "diffusion_1D_".fi.".eps"
set title "nx = ".gprintf("%02.0f",nx)
set yrange [-2:6]
set format y "%.0f"
set ylabel source_label 
plot "diffusion_1D_".fi.".dat" using 1:(log10($4)) title "1D solution", \
     "diffusion_1D_".fi.".dat" using 1:(log10($3)) title "diffusion solution" with lines ls 1, \
     "diffusion_1D_".fi.".dat" using 1:(log10($2)) title "" with lines ls 2
set yrange [-5:1]
set output "diffusion_1D_err_".fi.".eps"
set ylabel error_label
unset key
plot "diffusion_1D_".fi.".dat" using 1:(log10(abs($4-$3)/$3)) with lines

nx = 16
set key font ",15"
fi = gprintf("%04.0f",nx)
set output "diffusion_1D_".fi.".eps"
set title "nx = ".gprintf("%02.0f",nx)
set yrange [-2:6]
set format y "%.0f"
set ylabel source_label 
plot "diffusion_1D_".fi.".dat" using 1:(log10($4)) title "1D solution", \
     "diffusion_1D_".fi.".dat" using 1:(log10($3)) title "diffusion solution" with lines ls 1, \
     "diffusion_1D_".fi.".dat" using 1:(log10($2)) title "" with lines ls 2
set yrange [-5:1]
set output "diffusion_1D_err_".fi.".eps"
set ylabel error_label
unset key
plot "diffusion_1D_".fi.".dat" using 1:(log10(abs($4-$3)/$3)) with lines

