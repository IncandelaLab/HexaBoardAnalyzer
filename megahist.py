#!/usr/bin/env python2
import glob, sys, os, array, math
import numpy as np
from scipy.stats import norm
from scipy import signal
import string
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
plt.style.use('bmh')

def megahist(fname): #returns dist of rms for all channels, high and low gain
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
    OrgArr = np.zeros((chans,cols))#rows are channels, columns events and samples. 512 rows, 64 HG and LG for 4 chips 
 
    #any time ch hits 511 spl needs to add 13
    ch=0
    spl_last = 0
    for r in range(filerows):
        spl=0+spl_last
        for c in range(filecols):
            OrgArr[ch][spl]=data[r][c]
            spl+=1
        if ch < chans-1:
            ch+=1
        else: 
            ch=0#reset channel
            spl_last+=13

    hg_unirms=5.535; lg_unirms=5.858
    #resort OrgArr
    ci=0
    for chan in range(OrgArr.shape[0]):
        si=0
        for value in range(OrgArr.shape[1]):
            #print ci, si
            #print OrgArr[ci][si]
            if (si+1)%13==0:#take mean!
                #print "13 factor ", si+1
                mean_list=[]; j=0
                while j<13:
                    mean_list.append(OrgArr[ci][si-j]) #spl, spl-1...sp-12
                    #print OrgArr[ci][si-j], mean_list
                    j+=1
                mean = np.mean(mean_list)
                #print "mean ", mean, len(mean_list)

		j=0
                while j<13:
                    OrgArr[ci][si-j]-=mean #spl, spl-1...sp-12
                    #print OrgArr[ci][si-j], mean_list
                    if (ci<64) or (ci>=128 and ci<192) or (ci>=256 and ci<320) or (ci>=384 and ci<448):#lg
                        OrgArr[ci][si-j]=OrgArr[ci][si-j]/lg_unirms
                    else: #hg
                        OrgArr[ci][si-j]=OrgArr[ci][si-j]/hg_unirms
                    #print OrgArr[ci][si-j]
                    j+=1
            si+=1
        ci+=1
    #print OrgArr.shape
    #print OrgArr[:5][:15]

    #loop through and sort 1 more time!
    LG_Arr=np.zeros((256,OrgArr.shape[1])); HG_Arr=np.zeros((256,OrgArr.shape[1]))
    ci=0; li=0; hi=0
    #print "ci", ci              
    for chan in range(OrgArr.shape[0]):
        si=0
        for value in range(OrgArr.shape[1]):
            if (ci<64) or (ci>=128 and ci<192) or (ci>=256 and ci<320) or (ci>=384 and ci<448):#lg
                #print "LG", li, ci, si
                LG_Arr[li][si]=OrgArr[ci][si]
                lgbool=True
            elif (ci>=64 and ci<128) or (ci>=192 and ci<256) or (ci>=320 and ci<384) or (ci>=448): #hg
                #print "HG", hi, ci, si
                HG_Arr[hi][si]=OrgArr[ci][si]
                lgbool=False
            si+=1
        ci+=1
        if lgbool:
            li+=1
        else:
            hi+=1
    # print "lg", LG_Arr
    # print "hg", HG_Arr
    # print "org",OrgArr

    return LG_Arr, HG_Arr
#######################################################################
if __name__ == "__main__":
    #LGMegaArr=np.empty((256,482560000)); HGMegaArr=np.empty((256,482560000))
    LGMegaArr=[];HGMegaArr=[]
    if len(sys.argv) > 1:
        print "Number of files: ", len(sys.argv)
        print '# Input files are',sys.argv
        i=0
        for file in sys.argv:
            LGarr=[];HGarr=[]
            #LGbig=[];HGbig=[]
            if i>0: #there is a problem with the data for C11
                print i
                LGarr, HGarr = megahist(file)
                if i==1:
                    LGMegaArr=LGarr; HGMegaArr=HGarr;
                    print LGarr
                else:
                    LGMegaArr=np.concatenate((LGMegaArr, LGarr), axis=1)
                    HGMegaArr=np.concatenate((HGMegaArr, HGarr), axis=1)
                    print "shape", LGMegaArr.shape, HGMegaArr.shape
            i+=1
            #print LGMegaArr, HGMegaArr
    else:
        print "No input files given!"
        exit(0)

    LGMegaArr=np.ravel(LGMegaArr)
    HGMegaArr=np.ravel(HGMegaArr)
    print "shape", LGMegaArr.shape
    print LGMegaArr, HGMegaArr
    fig, ax = plt.subplots(1,2)
    nl, binsl, pl = ax[0].hist(LGMegaArr, 100, normed=0, facecolor='cornflowerblue',alpha=0.7)
    nh, binsh, ph = ax[1].hist(HGMegaArr, 100, normed=0, facecolor='indianred', alpha=0.7)

    LG_med=np.median(LGMegaArr)
    HG_med=np.median(HGMegaArr)
    ax[0].axvline(LG_med, color='k', linestyle='dashed', linewidth=2)
    ax[1].axvline(HG_med, color='k', linestyle='dashed', linewidth=2)
    
    print "lgmed", LG_med , "hgmed", HG_med
    print "number entries", len(LGMegaArr), len(HGMegaArr)
    ax[0].set_xlabel('LG noise for all channels, all pcbs', fontsize=20)
    ax[1].set_xlabel('HG noise for all channels, all pcbs', fontsize=20)
    plt.show()

 
