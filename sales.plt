# example of plotting salesman's trajectory
set title "Gnuplot of our Salesman's Cities"
unset key
unset border
set xlabel "longitude"
set ylabel "latitude"
set style line 1 lc rgb 'gray' lt 1 lw 1

# Replace cities2.dat with your solved ordering, or add an
# additional overlay to show the solution
# use range [-180:180][-90:90] to plot the whole globe
#
# plot world map and various city locations as an 'x'
plot [-170:-60][10:70] "world.dat" with filledcu ls 1, \
"cities.dat" using 1:2 with linespoints lw 2 lc 7


#make a png file of this plot
set term push
set term png
set output "sales_NothAm.png"
replot
set output
set term pop


