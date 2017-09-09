
import numpy as np
import matplotlib.pyplot as plt
import time
#import random poisson process data generators
from generators import *

#import bayesian blocks library
import BB

#set seed
np.random.seed(1)

#create some random data
t1 = linearPoisson(t0=0, start=10, end=20, length=5)
t2 = linearPoisson(t0=0, start=20, end=25, length=12)+t1[-1]
t3 = linearPoisson(t0=0, start=25, end=10, length=500)+t2[-1]
t = np.hstack((t1,t2,t3))

t2 = linearPoisson(t0=0, start=0, end=50, length=10)
plt.scatter(t2[1:],1/(t2[1:]-t2[:-1]),color="green")
plt.ylim(0,300)
plt.xlabel("time")
plt.ylabel("cell intensity")
plt.plot([0,10],[0,50],color="k",lw=6)
plt.title("Estimated Intensity for each cell")

for i in range(len(t2)):
    plt.plot([t2[i],t2[i]],[-1,1],lw=2,color="k")
#plt.scatter(t2[1:],np.repeat(0,len(t2)-1),marker=2,linewidths=2)
plt.xlabel("time")
plt.ylim(-0.1,0.1)
plt.xlim((0,10))
#run the bayesian blocks algorithm using linear and constant blocks(takes a couple of minutes to run)
q=BB.BayesianBlocks(t,c=4,type="linear",verbose=True,force_intercept=False)

start=time.time()
r=BB.BayesianBlocks(t,c=1 ,type="constant",verbose=True,PELT=True,steps=False)
print(str(time.time()-start)+" seconds")
####### PLOTTING #########


#plot points
plt.hist(t,color="gray",bins=40)

#plot the linear blocks
for i in range(len(q.blocks)):
    linear,=plt.plot([t[q.left[i]],t[q.right[i]]],[q.leftintensities[i],q.rightintensities[i]],color="red",lw=3)

#plot the constant blocks
for i in range(len(r.blocks)):
    constant,=plt.plot([t[r.left[i]],t[r.right[i]]],[r.intensities[i],r.intensities[i]],color="green",lw=3)


#plot the true intensities
plt.plot([t[0],5],[10,20],color="black",lw=3,linestyle="dashed")
plt.plot([5,17],[20,25],color="black",lw=3,linestyle="dashed")
true,=plt.plot([17,37],[25,50],color="black",lw=3,linestyle="dashed")


plt.ylim(0,70)
#plt.xlim(0,6)
plt.legend([constant,linear,true],["BB Estimated Constant Blocks","BB Estimated Linear Blocks","True Simulation Rates"],loc=2)
plt.title("Constant vs Linear Block Segmentation")
plt.xlabel("Time")
plt.ylabel("Intensity")

