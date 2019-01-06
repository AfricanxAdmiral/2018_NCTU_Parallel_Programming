import matplotlib.pyplot as plt
import numpy as np
import sys

file_serial = open('time_serial.txt', 'r')
file_pthread = open('time_pthread.txt', 'r')
file_openmp = open('time_openmp.txt', 'r')

xs = []
files = [file_serial, file_pthread, file_openmp]
times = []

for file in files:
    x = []
    time = []
    for line in file.read().splitlines():
        line = line.split(' ')
        if line[0]:
            x.append(float(line[0]))
            time.append(float(line[1]))
    times.append(time)
    xs.append(x)


plt.figure()
color = ['red', 'green', 'blue']
label = ['Serial', 'Pthread', 'Openmp']
for i in range(3):
    plt.plot(np.array(xs[i]), np.array(times[i]), label=label[i], color=color[i])
plt.xlabel('Data size')
plt.ylabel('Time')
plt.legend(loc='upper left')
plt.show()
