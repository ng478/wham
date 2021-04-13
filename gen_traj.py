import numpy as np
import matplotlib.pyplot as plt

#makes the data
y1 = np.random.normal(-2, 2, 1000)
y2 = np.random.normal(2, 2, 1000)
colors = ['b','g']

outF = open("frm0_10-run.traj", "w")
outF.write('# step         dist                ')
outF.write("\n")
counter = 0
for d in y1:
    c = str(counter)
    e = str(d)
    #print(c+'    '+e)
    outF.write(c+'    '+e)
    outF.write("\n")
    counter += 1
outF.close()

outF = open("frm1_10-run.traj", "w")
outF.write('# step         dist                ')
outF.write("\n")
counter = 0
for d in y2:
    c = str(counter)
    e = str(d)
    #print(c+'    '+e)
    outF.write(c+'    '+e)
    outF.write("\n")
    counter += 1
outF.close()
#plots the histogram
fig, ax1 = plt.subplots()
ax1.hist([y1,y2],color=colors)
ax1.set_xlim(-10,10)
ax1.set_ylabel("Count")
plt.tight_layout()
plt.show()