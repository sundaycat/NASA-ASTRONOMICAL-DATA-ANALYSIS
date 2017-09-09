from sklearn.linear_model import LinearRegression
from scipy.optimize import minimize
import numpy as np
import sys
import os


############################ IGNORE FOR NOW ################################################
#Dr. Scargle's standard constant block fitness function
def constantBlockLikelihood(N, M, c):
    with np.errstate(all='ignore'):  # when N or M are zero, ignore the divide by zero warning and return 0
        return np.nan_to_num(N * np.log(np.array(N) / M)-c )


def linearBlockLikelihood( x,t,c=0):
    if len(t)==2:
        return -c

    slope=x[0]
    intercept=x[1]
    lam = (slope * ((t[:-1] + t[1:]) / 2) + intercept) * (np.diff(t))
    with np.errstate(all='ignore'):  # when N or M are zero, ignore the divide by zero warning and return 0
        LL=np.sum(np.log(lam) - lam)-c
        if LL ==np.nan:
            return -np.inf
        else:
            return LL



def linearLeastSquares( x,t,c=0):
    if len(t)==2:
        return -c
    slope=x[0]
    intercept=x[1]
    return -np.log(np.sum(np.power(1/(t[:-1] - t[1:])-((t[:-1]-t[0])*x[0]+x[1]),2)))-c



def NelderMead(t,intercept=None):
    if len(t)==2:
        return(0,1/(t[1]-t[0]))


    if intercept is None:
        ans = minimize(linearBlockLikelihood, x0=[1, 1], args=(t), method="Nelder-Mead",
                    options={'xtol': 1e-4})
    else:
        ans = minimize(linearBlockLikelihood, x0=1, args=(intercept,t), method="Nelder-Mead",
                       options={'xtol': 1e-4})

    return ((ans.x[1], ans.x[0]))




def unbinnedPoisson(t,intercept=None):
    if len(t)==2:
        return(0,1/(t[1]-t[0]))

    if intercept is None:
        b=1
    else:
        b = intercept
    slope = 1

    s = (t[1:] + t[:-1]) / 2
    for iterations in range(100):
        lastslope = slope
        lastb = b
        slope = slope * np.sum(s / (slope * s + b)) / (t[-1] * t[-1] / 2)
        if intercept is None:
            b = b * np.sum(1 / (slope * s + b)) / (t[-1])
        with np.errstate(all='ignore'):
            if np.abs(1 - (lastslope / (slope+0.01))) < 0.03 and np.abs(1 - (lastb / (b+0.01))) < 0.03:
                break
    return (slope,b)

def gradientAscent(t,intercept=None):

    a=100
    b=100
    s=(t[1:]+t[:-1])/2
    d=(t[1:]-t[:-1])
    lrb=100
    lra=100
    lastgb=1
    lastga=1
    for i in range(100):
        lasta=a
        lastb=b
        newgb=np.sum(1 / (a * s + b) - d)
        newga=np.sum(s / (a * s + b) - s * d)
        if lastgb>newgb:
            lrb*=1.2
        else:
            lrb/=3
        if lastga>newga:
            lra *= 1.2
        else:
            lra /= 3
        b += lrb*newgb
        a+=lra*newga
        lastgb = newgb
        lastga= newga
        if np.abs(lastb-b)<0.1 and np.abs(lasta-a)<0.1:
            break
    return (a,b)


def binnedPoisson(t,intercept=None):
    if len(t)==2:
        return(0,1/(t[1]-t[0]))

    if len(t)<11:
        return (1,1)

    intercept_change=1
    slope_change=1
    bins = int(np.floor(len(t) / 10))
    binwidth = (t[-1] - t[0]) / bins
    y = np.zeros(bins)
    for i in range(bins):
        y[i] = np.sum(np.logical_and(t > (i * binwidth), t < ((i + 1) * binwidth)))

    y /= binwidth
    x = np.arange(bins) * binwidth

    #y = np.reshape(y, (len(y), 1))
    y=np.ravel(y)
    x=np.ravel(x)
    x = np.reshape(x, (len(x), 1))

    m = LinearRegression()
    m.fit(x, y)
    iter=0
    while(np.abs(intercept_change)>0.01 or np.abs(slope_change)>0.01):

        last_intercept = m.intercept_
        last_slope = m.coef_
        weights = np.abs((1 / (m.intercept_ + x * m.coef_)))# / np.sum(1 / (m.intercept_ + x * m.coef_))
        weights[weights ==np.inf] = np.nanmin(np.isfinite(weights))
        weights=np.ravel(weights)
        m.fit(x, y, sample_weight=weights)
        intercept_change=1-last_intercept/(m.intercept_+0.01)
        slope_change=1-last_slope/(m.coef_+0.01)
        iter+=1
        if iter>100:
            break

    return (m.coef_[0],m.intercept_)
########################################################################################

#object consisting of the parameters of a single block
class Block(object):
    def __init__(self,type,**kwargs):
        self.left = kwargs["left"]
        self.right = kwargs["right"]
        self.LL = kwargs["LL"]
        self.type=type

        if type=="linear":
            self.leftintensity=kwargs["leftint"]
            self.rightintensity=kwargs["rightint"]

        elif type=="constant":
            self.intensity=kwargs["intensity"]

#object consisting of the blocks at a given step along with any other step information
class Step(object):
    def __init__(self,blocks, bestLL,pruned,LLs):
        self.blocks=blocks
        self.bestLL=bestLL
        self.pruned=pruned
        self.LLs=LLs


##### Bayesian Block Model Object for Poisson Process Data #########
class BayesianBlocks(object):
    def _bestLinear(self, i, c, method):

        parameters = np.zeros((i,2))
        LLs = np.ones(i) * -np.inf
        for p in (set(range(i)) - set(self.step[-1].pruned)):
        #for p in range(i):  # only compute if we have 3 or more points
            t = self.points[p:i + 1] - self.points[p]
            # print(t)
            # if self.force_intercept and self.step[p-1].blocks[-1].rightintensity:
            #    (slope, b) = method(t,self.step[p-1].blocks[-1].rightintensity)
            # else:
            if self.force_intercept==True and len(self.step[-1].blocks)>0:
                (slope, b) = method(t=t,intercept=self.step[-1].blocks[-1].rightintensity)
            else:
                (slope, b) = method(t=t)

            parameters[p,:]=[slope, b]
            LLs[p]=linearBlockLikelihood([slope, b], t, c)
            LLs[p] = linearLeastSquares([slope, b], t, c)

        totalLLs = np.array(LLs) + np.array([self.step[a].bestLL for a in range(i)])
        best_partition = np.argmax(totalLLs)

        bestslope = parameters[best_partition,0]
        bestintercept = parameters[best_partition,1]

        bestBlock = Block(type="linear", left=best_partition, right=i, leftint=bestintercept,
                          rightint=bestintercept + bestslope * (self.points[i] - self.points[best_partition]),
                          LL=LLs[best_partition])
        return (bestBlock, totalLLs)


    def _bestConstant(self,i,c):
        M = self.points[1:i+1] - self.points[:i]
        N = np.ones(len(M))
        blocksum = np.cumsum(M[::-1])[::-1]
        pointsum = np.cumsum(N[::-1])[::-1]
        #initially set all log likelihood values to -inf (used for pruning)


        if self.PELT==True:
            LLs = np.ones(i) * -np.inf
            points_to_compute=list((set(range(i)) - set(self.step[-1].pruned)))
            LLs[points_to_compute]=constantBlockLikelihood(pointsum[points_to_compute], blocksum[points_to_compute], c)
        else:
            LLs = constantBlockLikelihood(pointsum, blocksum, c)

        #    LLs[k] = constantBlockLikelihood(pointsum[k], blocksum[k], c)
        #get the summed log likelihoods
        totalLLs=LLs+np.array([self.step[a].bestLL for a in range(i)])

        #find the partition providing the largest log likelihood sum
        best_partition = np.argmax(totalLLs)
        #for this partition, calculate the intensity of the new block
        best_block_intensity = pointsum[best_partition] / blocksum[best_partition]

        #create a new "Block" object containing the new blocks attributes
        bestBlock=Block(type="constant",left=best_partition,right=i,intensity=best_block_intensity,LL=LLs[best_partition])

        #return the best block as well as the summed log likelihoods
        return (bestBlock,totalLLs)



    #this code below runs when you initialize a "BayesianBlock" object
    def __init__(self,points,c=0.5,type="constant",verbose=False,force_intercept=False,PELT=True,steps=True):
        self.points=np.append(0,points)
        self.PELT=PELT
        # if verbose is turned off(False) then turn off any printing output
        if verbose==False:
            sys.stdout = open(os.devnull, 'w')

        #create the first step with an empty block(zero points has log likelihood 0)
        self.step=[Step(blocks=[],bestLL=0,pruned=[],LLs=[])]
        self.force_intercept=force_intercept

        #start iterating through all points
        for i in np.arange(1,len(points)-1):
            #print("point: "+str(i) + " of "+str(len(points)-2)+ " pruned: "+ str(len(self.step[-1].pruned)/i))

            #find the best constant block given the first "i" points


            if type=="linear":
                (bestblock,totalLLs)=self._bestLinear(i,c=c,method=binnedPoisson)
            elif type=="constant":
                (bestblock, totalLLs) = self._bestConstant(i, c)

            #if the best block is the whole partition, append a "step" object with the single block
            if bestblock.left==0:
                self.step.append(Step([bestblock],bestblock.LL,pruned=np.where((totalLLs+c)<bestblock.LL)[0],LLs=totalLLs))

            #if a new block is found, append it to the best combination behind the new block
            else:
                newblocks=self.step[bestblock.left-1].blocks+[bestblock]
                newbestLL=self.step[bestblock.left-1].bestLL + bestblock.LL
                self.step.append(Step(newblocks,newbestLL,pruned=np.where((totalLLs+c)<newbestLL)[0],LLs=totalLLs))


        #save the last step's result
        self.blocks=self.step[-1].blocks

        #turn the printing back on
        sys.stdout = sys.__stdout__



    #convenient object properties
    @property
    def bestLLs(self):
        return np.array([self.step[a].bestLL for a in range(len(self.step))])

    @property
    def types(self):
        return np.array([self.blocks[a].type for a in range(len(self.blocks))])

    @property
    def intensities(self):
        return np.array([self.blocks[a].intensity for a in range(len(self.blocks))])

    @property
    def rightintensities(self):
        return np.array([self.blocks[a].rightintensity for a in range(len(self.blocks))])

    @property
    def leftintensities(self):
        return np.array([self.blocks[a].leftintensity for a in range(len(self.blocks))])

    @property
    def left(self):
        return np.array([self.blocks[a].left for a in range(len(self.blocks))])

    @property
    def right(self):
        return np.array([self.blocks[a].right for a in range(len(self.blocks))])



