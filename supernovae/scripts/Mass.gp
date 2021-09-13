set termoption enhanced
set grid
set style line 1 lc rgb "blue" lt 1 lw 2 pt 7 ps 0.5
set style line 2 lc rgb "red" lt 1 lw 2 pt 7 ps 0.5
set style line 3 lc rgb "green" lt 1 lw 2 pt 7 ps 0.5
set style line 4 lc rgb "orange" lt 1 lw 2 pt 7 ps 0.5

set title "Mass vs. Time
set xlabel 'time'
set ylabel 'mass'
set terminal jpeg enhanced
set output 'results/mass.jpg'
set datafile separator ","
plot 'data/Mass.dat' u 1:2 w l ls 1 t 'Star 1', 'data/Mass.dat' u 1:3 w l ls 2 t 'Star 2'