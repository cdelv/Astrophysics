set termoption enhanced
set grid
set style line 1 lc rgb "blue" lt 1 lw 2 pt 7 ps 0.5
set style line 2 lc rgb "red" lt 1 lw 2 pt 7 ps 0.5
set style line 3 lc rgb "green" lt 1 lw 2 pt 7 ps 0.5
set style line 4 lc rgb "orange" lt 1 lw 2 pt 7 ps 0.5

set title "System energy vs. Time
set xlabel 'Time'
set ylabel 'Energy'
set terminal jpeg enhanced
set output 'results/tenergy.jpg'
set datafile separator ","
plot 'data/TotalEnergy.dat' u 1:2 w l ls 1 t 'K energy', 'data/TotalEnergy.dat' u 1:3 w l ls 2 t 'P energy', 'data/TotalEnergy.dat' u 1:4 w l ls 3 t 'Total E' 