import numpy as np
import math as mt
import matplotlib.pyplot as plt


size=1000;maxTime=1400
ez=np.zeros((size))
hy=np.zeros((size))
ez_snap=np.zeros((maxTime))
hy_snap=np.zeros((maxTime))
imp0=377.0
snap_moment=size/2-1

source_width = 5.0
delay = 10*source_width
#for snap_moment in range(5,10,1):
i=0
#while (count < 9)
for qtime in np.arange(0,maxTime+1,1):

    print(qtime)
    ez[size/2-1] += 0.5*mt.exp(-(qtime-delay) ** 2 / (2.0 * source_width**2));

    for mm in range(0,size-1):
        #  /* update magnetic field */
        hy[mm]=hy[mm]+(ez[mm+1]-ez[mm])/imp0
    for mm in range(1,size):
         #/* update electric field */
        ez[mm] = ez[mm] + (hy[mm] - hy[mm - 1]) * imp0
    #ez_snap[i]=ez[snap_moment]
    #hy_snap[i]=hy[snap_moment]
   # if((max(ez)>0)&(max(ez)<1e-15)&(qtime>50)): break
    i=i+1
    #if (ez[size/2-1]<-0.99 and ez[size/2-1]>-1) and (hy[size/2-1]<0.1 and hy[size/2-1]>0): break

fig = plt.figure()
#plt.plot(range(0,maxTime,1),hy_snap,'-')
plt.subplot(2, 1, 1)
plt.plot(range(0,size,1),ez,'-')
plt.plot(range(0,size,1),hy*imp0,'-')
plt.ylabel('Ez-field, V/m')
plt.xlabel('Spatial')
#plt.axis([0, 6, 0, 20])

plt.subplot(2, 1, 2)
plt.plot(range(0,size,1),ez,'-')
plt.ylabel('Ez-field, V/m')
plt.xlabel('Spatial')
plt.show()

fig.savefig('/Users/nikitapavlov/Library/Mobile Documents/com~apple~CloudDocs/University/Kostya/step0_inizavanish.pdf')
 #plt.pause(0.5)

