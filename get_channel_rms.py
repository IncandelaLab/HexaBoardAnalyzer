#!/usr/bin/env python2
import glob, sys, os, array, math
import numpy as np
from scipy.stats import norm
from scipy import signal
import string
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
plt.style.use('bmh')

def get_channel_rms(fname, showplot=True): #returns dist of rms for all channels, high and low gain
#need rms of all channels -samples and events. 64 channels per event
#rms is sqrt variance or standard deviation 
#read through data and store in 2 256 row arrays for high and low gain
#every 64 data lines, read in 13 values, append to same row as 64 counter     
    opfl = open(fname,'r')
    cont = opfl.read()
    opfl.close()
    
    lines = cont.splitlines()[8:]
    print(len(lines))
    lines = [l for l in lines if not l.startswith('E')]

    opfl = open('data.txt','w')
    opfl.write('\n'.join(lines))
    opfl.close()

    data = np.loadtxt("data.txt")
    data = data[:,:-2]
    filerows = data.shape[0]
    filecols = data.shape[1]
    events = int(filerows / 512)
    cols = events*13 
    chans = 512
    print events, cols
    print "data",data
    print data.shape
    OrgArr = np.zeros((chans,cols))#rows are channels, columns events and samples. 512 rows, 64 HG and LG for 4 chips 
 
    #any time ch hits 511 spl needs to add 13
    ch=0
    spl_last = 0
    for r in range(filerows):
        spl=0+spl_last
        for c in range(filecols):
            #print data[r][c]
            #print r, ch
            #print ch, spl, r, c
            OrgArr[ch][spl]=data[r][c]
            spl+=1
        if ch < chans-1:
            ch+=1
        else: 
            ch=0#reset channel
            spl_last+=13
    print "OrgArr", OrgArr
    print spl, ch
    print OrgArr.shape

    #NOTE: MAKE SURE YOUR RMS IS RIGHT. it doesn't match sebastian's plots- both use np.std() for rms
    #MAKE SURE YOUR ORGARR IS RIGHT
    #get the rms of every row. every 64 alternate high gain/low gain 
    print "first channel ", OrgArr[0][:]
    print "first channel rms ", OrgArr[0][:].std()
    allRMS = np.array(OrgArr).std(1) #array of std of every channel
    print "allRMS shape", allRMS.shape
    print "first channel check ", allRMS[0]
    hgRMS=[]; lgRMS=[]
    i=0
    for val in allRMS:
        if i<64: #wait but what is the order... hg lg or lg hg? i dont know, good to ask
            lgRMS.append(val)
        if i>=64:
            hgRMS.append(val)
        i+=1
        if i==128:
            i=0
    print "len(RMS)", len(hgRMS), len(lgRMS)
    #print 'lg \n', lgRMS[:26], '\n hg \n', hgRMS[:26], '\n all \n', allRMS[:26] 

    #plot high and low gain distributions and fit with a gaussian
    #should always cut at rms of 10 and only fit up to 10
    lgcut=[];hgcut=[]
    for val in lgRMS:
        if val < 10.:
            lgcut.append(val)
    for val in hgRMS:
        if val < 10.:
            hgcut.append(val)

    mul, sigl = norm.fit(lgcut)
    muh, sigh = norm.fit(hgcut)

    fig, ax = plt.subplots(1,2)
    nl, binsl, pl = ax[0].hist(lgRMS, 50, normed=1, facecolor='cornflowerblue',alpha=0.7)
    nh, binsh, ph = ax[1].hist(hgRMS, 50, normed=1, facecolor='indianred', alpha=0.7)

    yl = mlab.normpdf( binsl, mul, sigl)
    fitl = ax[0].plot(binsl, yl, '--', linewidth=3)
    yh = mlab.normpdf( binsh, muh, sigh)
    fith = ax[1].plot(binsh, yh, '--', linewidth=3)
    
    # ax[0].xlabel('lgRMS')
    # ax[1].xlabel('hgRMS')
    print mul+sigl, muh+sigh
    ax[0].set_xlabel('LGrms', fontsize=20)
    ax[1].set_xlabel('HGrms', fontsize=20)
    plt.show()
    
    #next would be getting the sigmas and returning the various thresholds, plus the over 10 thr
    #the index of the allRMS array is the channel number 
    #can also maybe find the sigma for each rms value and plot noisy channels that way, continuous x axis for thr
    LGthr=[];HGthr=[] #these are the values of rms for the thr
    LGsig=[];HGsig=[]#these are the sigma values
    i=1
    while i<=5:
        f=float(i)
        #print i,f
        LGsig.append(f); HGsig.append(f)
        LGthr.append(mul+f*sigl)
        HGthr.append(muh+f*sigh)
        f=float(i)+0.5
        #print i,f
        LGsig.append(f); HGsig.append(f)
        LGthr.append(mul+f*sigl)
        HGthr.append(muh+f*sigh) 
        i+=1
    #LGthr.append(10.); HGthr.append(10.)
    #LGthr[9]=10.; HGthr[9]=10.
    #what sigma gives you 10?
    #what number of channels are noisy? this is NOT identifying the channels, it is just giving the overall number of noisy ones at each std
    LGnc=[]; HGnc=[]
    for thr in LGthr:
        ch_count=0
        for rms in lgRMS:
            if rms >= thr:
                ch_count+=1
        LGnc.append(ch_count)
    for thr in HGthr:
        ch_count=0
        for rms in hgRMS:
            if rms >= thr:
                ch_count+=1
        HGnc.append(ch_count)
    print HGnc
    print LGnc

    #plot # noisy channels as a function of SIGMA (thr)
    fig, ax = plt.subplots(1,2)
    # n, b, p = ax[0].hist(LGnc, 10, facecolor='cornflowerblue')
    # n, b, p = ax[1].hist(HGnc, 10,facecolor='indianred')
    width = .5 
    ax[0].bar(LGsig, LGnc, width, color='cornflowerblue', alpha=0.7)
    ax[1].bar(HGsig, HGnc, width, color='indianred', alpha=0.7)
    ax[0].set_xlabel('LGsigma', fontsize=20)
    ax[1].set_xlabel('HGsigma', fontsize=20)
    ax[0].set_ylabel('noisy channel count', fontsize=20)
    ax[1].set_ylabel('noisy channel count', fontsize=20)
    ax[0].xaxis.set_ticks(np.arange(1, 6, 0.5))
    ax[1].xaxis.set_ticks(np.arange(1, 6, 0.5))
    plt.show()
    print len(lgRMS)

#######################################################################
if __name__ == "__main__":

    if len(sys.argv) > 1:
        myfile = sys.argv[1]
        print '# Input files are', myfile
    else:
        print "No input files given!"
        exit(0)
    get_channel_rms(myfile)



