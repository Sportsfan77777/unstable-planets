
import numpy as np
from matplotlib import pyplot as plot
import math

import pickle

def geo_mean(array):
	gm = 1
	power = 1.0 / len(array)
	for a in array:
		gm *= a

	gm = gm ** power
	return gm

def log_mean(array):
	log_array = []
	for a in array:
		loga = math.log(a) / math.log(10.0)
		log_array.append(loga)

	log_avg = np.mean(log_array)
	final_avg = 10.0 ** log_avg
	return final_avg


x = pickle.load(open("sm_axes.p", "rb"))

ejectionTable = pickle.load(open("table_of_ejections.p", "rb"))

y = np.zeros(len(x))
y2 = np.zeros(len(x))
y3 = np.zeros(len(x))
for i in xrange(len(x)):
	row = [1000 * q for q in ejectionTable[i,:]]

	avg = np.mean(row)
	y[i] = avg

	geo = geo_mean(row)
	y2[i] = geo

	logm = log_mean(row)
	y3[i] = logm


last = -9
x = x[:last]
y = y[:last]
y2 = y2[:last]
y3 = y3[:last]

######################### PLOT BELOW ##############################


plot.title("Arithmetic Mean of Ejection Times")
plot.xlabel("semimajor axis (in scaled AU)")
plot.ylabel("AM of Ejection Time (in years)")


plot.plot(x, y)
plot.plot(x, y, 'ro')
plot.yscale('log')

plot.savefig("plot_ejectionTimes_AM.png")
plot.show()


plot.title("Geometric Mean of Ejection Times")
plot.ylabel("GM of Ejection Time (in years)")


plot.plot(x, y2)
plot.plot(x, y2, 'ro')
plot.yscale('log')

plot.savefig("plot_ejectionTimes_GM.png")
plot.show()


plot.title("Logarithmic Mean of Ejection Times")
plot.ylabel("LM of Ejection Time (in years)")


plot.plot(x, y3)
plot.plot(x, y3, 'ro')
plot.yscale('log')

plot.savefig("plot_ejectionTimes_LM.png")
plot.show()




