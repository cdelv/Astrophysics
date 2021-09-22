set termoption enhanced
set grid
set style line 1 lc rgb "blue" lt 1 lw 2 pt 7 ps 0.5

set title "CM velocity magnitud vs. Time
set xlabel 'Time'
set ylabel 'UNITS'
set terminal jpeg enhanced
set output 'results/CMvel.jpg'
set datafile separator ","
plot 'data/CMVelocity.dat' u 1:5 w l ls 1 t 'Velocity magnitud'