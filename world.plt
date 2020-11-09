# note: due to limitations in gnuplot, the lines that are meant to
# "wrap around the back" of the map will not display correctly

# example of plotting salesman's trajectory
#set title "Gnuplot of our Salesman's Cities"
unset key
unset border

#set xlabel "longitude"
#set ylabel "latitude"
set style line 1 lc rgb 'gray' lt 1 lw 1
set tics scale 0
unset ytics
unset xtics

set arrow from 0,-86 to 0,-82 nohead lc -1
set arrow from -60,-86 to -60,-82 nohead lc -1
set arrow from 60,-86 to 60,-82 nohead lc -1


# Replace cities4.dat with your solved ordering, or add an
# additional overlay to show the solution
# use range [-180:180][-85:85] to plot the whole globe
#
# plot world map and various city locations as an 'x'
plot [-180:180][-85:85] "world_50m.txt" with filledcu ls 1, \
abs(x)>170?0:1/0 lc 0 dt 3, abs(x)>170?40:1/0 lc 0 dt 3, \
abs(x)>170?-40:1/0 lc 0 dt 3, \
"cities4.dat" using 1:2 with linespoints lw 2 lc 7




#make a png file of this plot
set term push
set term pdf
set output "world.pdf"
replot
set output
set term pop


