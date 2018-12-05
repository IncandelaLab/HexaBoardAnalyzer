#!/usr/bin/env python2
import glob, sys, os, array, math
import numpy as np
from scipy.stats import norm
from scipy import signal
import string
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
plt.style.use('bmh')

def chip_median_rms(fname, showplot=False, medplot=False): #returns dist of rms for all channels, high and low gain
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
    print "chan 58 chip 4 check", OrgArr[384+58][:]
    #chk=OrgArr[0][:]
    chk=OrgArr[384+58][:]
    #NOTE: MAKE SURE YOUR RMS IS RIGHT. it doesn't match sebastian's plots- both use np.std() for rms
    allRMS = np.array(OrgArr).std(1) #array of std of every channel
    hgRMS=[]; lgRMS=[]
    i=0
    for val in allRMS:
        if i<64: 
            lgRMS.append(val)
        if i>=64:
            hgRMS.append(val)
        i+=1
        if i==128:
            i=0
    #print "len(RMS)", len(hgRMS), len(lgRMS)

    ch1_hgRMS=[]; ch1_lgRMS=[]; ch2_hgRMS=[]; ch2_lgRMS=[]; ch3_hgRMS=[]; ch3_lgRMS=[]; ch4_hgRMS=[]; ch4_lgRMS=[]; 
    i=0
    for val in lgRMS:
        if i<64:
            ch1_lgRMS.append(val)
        if i>=64 and i<128:
            ch2_lgRMS.append(val)
        if i>=128 and i<192:
            ch3_lgRMS.append(val)
        if i>=192 and i<256:
            ch4_lgRMS.append(val)
        i+=1
    i=0
    for val in hgRMS:
        if i<64:
            ch1_hgRMS.append(val)
        if i>=64 and i<128:
            ch2_hgRMS.append(val)
        if i>=128 and i<192:
            ch3_hgRMS.append(val)
        if i>=192 and i<256:
            ch4_hgRMS.append(val)
        i+=1
    
    #get the medians
    LG_med=np.median(lgRMS)
    HG_med=np.median(hgRMS)
    ch1_lgmed=np.median(ch1_lgRMS); ch2_lgmed=np.median(ch2_lgRMS); ch3_lgmed=np.median(ch3_lgRMS); ch4_lgmed=np.median(ch4_lgRMS); 
    ch1_hgmed=np.median(ch1_hgRMS); ch2_hgmed=np.median(ch2_hgRMS); ch3_hgmed=np.median(ch3_hgRMS); ch4_hgmed=np.median(ch4_hgRMS); 
    # print "LG: ch1 ", ch1_lgmed,"ch2 ", ch2_lgmed,"ch3 ", ch3_lgmed,"ch4 ", ch4_lgmed
    # print "HG: ch1 ", ch1_hgmed,"ch2 ", ch2_hgmed,"ch3 ", ch3_hgmed,"ch4 ", ch4_hgmed
    lgm = [ch1_lgmed, ch2_lgmed, ch3_lgmed, ch4_lgmed]; hgm = [ch1_hgmed, ch2_hgmed, ch3_hgmed, ch4_hgmed]
    print "LG: ", lgm; print "HG: ", hgm

    nl, binsl, pl = plt.hist(chk, 50, normed=0, facecolor='magenta',alpha=0.7)
    plt.show()
    #plot the distribution of rmses for all channels
    if showplot:
        fig, ax = plt.subplots(1,2)
        nl, binsl, pl = ax[0].hist(lgRMS, 50, normed=0, facecolor='cornflowerblue',alpha=0.7)
        nh, binsh, ph = ax[1].hist(hgRMS, 50, normed=0, facecolor='indianred', alpha=0.7)
        ax[0].axvline(LG_med, color='k', linestyle='dashed', linewidth=2)
        ax[1].axvline(HG_med, color='k', linestyle='dashed', linewidth=2)
        print "lgmed", LG_med , "hgmed", HG_med
        ax[0].set_xlabel('LGrms', fontsize=20)
        ax[1].set_xlabel('HGrms', fontsize=20)
        #plt.show()
        plt.close()
    #find threshold, defined as 3x,4x, or 5x the median rms for this pcb/chip, and find number of channels greater than that threshold
    #ch1_lgthr=[]; ch2_lgthr=[]; ch3_lgthr=[]; ch4_lgthr=[]; ch1_hgthr=[]; ch2_hgthr=[]; ch3_hgthr=[]; ch4_hgthr=[]; 
    facs=[2.0,3.0,4.0,5.0]
    fxmedlg=np.zeros((len(facs),4)); fxmedhg=np.zeros((len(facs),4)) #will end up with 4 channel values for 4 chips
    #print "shape check", fxmedhg.shape
    for l0,l1,l2,l3,h0,h1,h2,h3 in zip(ch1_lgRMS, ch2_lgRMS, ch3_lgRMS, ch4_lgRMS, ch1_hgRMS, ch2_hgRMS, ch3_hgRMS, ch4_hgRMS):
        j=0
        for f in facs:
            #print f,j
            if l0>=f*lgm[0]: fxmedlg[j][0]+=1; print l0, f*lgm[0], "chan", ch1_lgRMS.index(l0)
            if l1>=f*lgm[1]: fxmedlg[j][1]+=1; print l1, f*lgm[1], "chan", ch2_lgRMS.index(l1)
            if l2>=f*lgm[2]: fxmedlg[j][2]+=1; print l2, f*lgm[2], "chan", ch3_lgRMS.index(l2)
            if l3>=f*lgm[3]: fxmedlg[j][3]+=1; print l3, f*lgm[3], "chan", ch4_lgRMS.index(l3)

            if h0>=f*hgm[0]: fxmedhg[j][0]+=1; print h0, f*hgm[0]
            if h1>=f*hgm[1]: fxmedhg[j][1]+=1; print h1, f*hgm[1]
            if h2>=f*hgm[2]: fxmedhg[j][2]+=1; print h2, f*hgm[2]
            if h3>=f*hgm[3]: fxmedhg[j][3]+=1; print h3, f*hgm[3]
            j+=1
    print "channel counts: LG ",fxmedlg, "HG", fxmedhg
    print "shape check", fxmedhg.shape
    if medplot: #return median rms
        return ch1_lgmed, ch2_lgmed, ch3_lgmed, ch4_lgmed, ch1_hgmed, ch2_hgmed, ch3_hgmed, ch4_hgmed
    else: #return # chans over median rms for each chip for hg and lg  
        return fxmedlg, fxmedhg
#######################################################################
if __name__ == "__main__":

    if len(sys.argv) > 1:
        myfile = sys.argv[1]
        print '# Input files are', myfile
    else:
        print "No input files given!"
        exit(0)
    chip_median_rms(myfile)
