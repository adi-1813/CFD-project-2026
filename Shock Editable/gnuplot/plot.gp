# ---------- Terminal & Output ----------
set terminal pngcairo size 1200,900 enhanced font "Arial,14"
set output "shocktube.png"

# ---------- General Settings ----------
set grid
set border lw 1.5
set tics font ",12"

set multiplot layout 3,1 title "Shock Tube Solution" font ",16"

# ---------- Density ----------
set title "Density"
set xlabel "x"
set ylabel "Density"
set yrange [0:1.2]
plot "shocktube.dat" using 1:2 with lines lw 2 lc rgb "blue" title "Density"

# ---------- Velocity ----------
set title "Velocity"
set xlabel "x"
set ylabel "Velocity"
unset yrange
plot "shocktube.dat" using 1:3 with lines lw 2 lc rgb "red" title "Velocity"

# ---------- Temperature ----------
set title "Temperature"
set xlabel "x"
set ylabel "Temperature"
set yrange [0.8:1.1]
plot "shocktube.dat" using 1:4 with lines lw 2 lc rgb "green" title "Temperature"

unset multiplot
