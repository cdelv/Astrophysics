set termoption enhanced
set grid
set style line 1 lc rgb "blue" lt 1 lw 2 pt 7 ps 0.5
set style line 2 lc rgb "red" lt 1 lw 2 pt 7 ps 0.5

set title 'Energy time derivative'
set xlabel 'Time'
set ylabel 'UNITS'
set terminal jpeg enhanced
set output 'results/Ederivative.jpg'
set datafile separator ","
plot 'data/Energy_derivative.dat' u 1:2 w l ls 1 t 'Total system energy derivative', 'data/Energy_derivative.dat' u 1:3 w l ls 2 t 'Stars energy derivative'