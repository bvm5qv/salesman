Starter code and data for traveling salesman problem


Files in this directory:

* datareader.cpp : example code to read in the data files (use Makefile)
* datareader.py  : example code to read in the data files
* cities23.dat : list of coordinates for 23 cities in North America
* cities150.dat : 150 cities in North America
* cities1k.dat : 1207 cities in North America
* cities2k.dat : 2063 cities around the world
* routeplot.py : code to plot the globe and salesman's path<br>
usage:<br>
python routeplot.py cities.dat [cities2.dat] -r [="NA"],"World"'<br>
NA = North America, World = Mercator projection of the whole earth
* earth.C : (just for fun) plotting the globe in ROOT


INSTRUCTIONS and COMMENTS:
The following command line arguments are taken to run sales.cpp: the data file of cities for the path about to be minimized, the initial temperature, the minimum temperature (under which the program halts), the cooling constant (between 0 and 1), and optionally the path of a data file to output the results of the computation (i.e. the order in which the cities are arranged using longitude and latitude for the minimum path between them).
For example, if we were finding the minimum path about the 1k cities, starting with an initial temperature of 5000km, chose a minimum temperature of 1.0E-5km, and decided upon a cooling constant of 0.995, we would run the file with the command:
/.sales cities1k.dat 5000 1.0E-5 0.995. 
On the other hand, if we wished to do the same for the 150 cities, but record our data in the file cities150OUT.dat, we would use the command: ./sales cities150.dat 5000 1.0E-5 0.995 cities150OUT.dat.

During the beginning of the annealing process, when the temperature is still hot, we wish for our program to be able to make large changes in energy - both up and down - to prevent falling into local minimums. As such, if dE is a typical change in energy during the initial "melt" of our system, the probability of the change which caused dE being accepted (if dE > 0) is roughly p = exp(-dE/T0). We want p to be large initially (say around 0.7), so we have that T0 = -dE/ln(p). 
The continental US is about 4300km wide, so when resitrcted to all of North America, T0 = 5000km is a good starting temperature following our previous logic. Similarly, the circumference of the Earth is about 40000km, so T0 = 25000km is a good initial temperature when dealing with cities across the entire Earth.

The minimum temperature should be chosen so that the probability of any typical increase in energy being accepted is small. As such, p = exp(-dE/Tmin) => Tmin = -dE/ln(p) for small p. It is clear that a value of about Tmin = 1.0E-5km works in any event.

A cooling constant of alpha = 0.995 was chosen because it cools the system slowly and allowes for many solutions to be tested (while still allowing for the program to resolve itself in a reasonable amount of time).

SOLUTIONS:
cities150 317299km 48443.4km 2.02101s
cities1k  732178km 98457.4km 16.8762s
cities2k  1.01876E+07km 295628km 34.6193s





