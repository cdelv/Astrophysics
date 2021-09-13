set termoption enhanced
set grid
set style line 1 lc rgb "blue" lt 1 lw 2 pt 7 ps 0.5
set style line 2 lc rgb "red" lt 1 lw 2 pt 7 ps 0.5
set style line 3 lc rgb "green" lt 1 lw 2 pt 7 ps 0.5
set style line 4 lc rgb "orange" lt 1 lw 2 pt 7 ps 0.5

set title "Baricenter trayectory
set xlabel 'x'
set ylabel 'y'
set zlabel 'z'
set terminal jpeg enhanced
set output 'results/baricenter.jpg'
set datafile separator ","
set parametric
set view 40, 250

splot 'data/Baricenter.dat' u 2:3:4 w l ls 1 t'All the system', 'data/Baricenter.dat' u 5:6:7 w l ls 2 t'Only stars'