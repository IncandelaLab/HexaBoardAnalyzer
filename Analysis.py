import glob, sys, os, array, math
import numpy as np
from scipy.stats import norm
from scipy import signal
import string
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
plt.style.use('seaborn-deep')
#plt.style.use('bmh')
import ROOT as rt
from readHexConvData import *
from plotNoise import *
#matplotlib.rcParams.update({'font.size': 22})

def read_root_file(filename):
    #reads the root tree
    #returns a nested list {[sample1(ch2 events),(ch1 events)],[sample2(ch1 events),(ch2 events)]...} 
    #rootfile stores channel by channel a lg and hg array of samples 
    LG_arr=[]; HG_arr=[]
    for sample in range(13):
        lgsam=[];hgsam=[]
        treedata=readTree(filename, chip='all',sca=sample)#sca is sample#
        chans = treedata.keys()
        gains = treedata[chans[0]].keys()
        all_chan_data = { chan:{gain:[] for gain in gains} for chan in chans}
        LG_ev_arr=[]; HG_ev_arr=[]
        for chan in chans:
            chandata = treedata[chan]
            #print chandata
            for gain,values in chandata.items():
                #print gain,values
                valarr=values
                if gain == 'lg' and len(values)!=0: LG_ev_arr.append(valarr)#tolist()
                if gain == 'hg' and len(values)!=0: HG_ev_arr.append(valarr)
            lgsam.append(LG_ev_arr); hgsam.append(HG_ev_arr)
        LG_arr.append(lgsam); HG_arr.append(hgsam)
    print "len arrays",len(LG_arr), len(HG_arr)
    return LG_arr, HG_arr

def get_chan_list(array, ped=False):
#returns list of lists of samples for all events for better data organization
#subtracts mean adc count (pedestal) from all events in 1 sample if ped=true and returns a list of the means per channel as well if ped=True
    print "array length (should be 13)", len(array)
    if ped:
        print "returning array of averages over each event"
        mean_list=[[]for Null in range(256)]
    chan_list=[[]for Null in range(256)]#list of 256 lists
    for num in range(256):#only get one channel at a time!
        for i, sam in enumerate(array):#loop over 13 samples
            #print "len sam",len(sam)
            for c,chan in enumerate(sam):
                if c==num: #only work with one channel at a time
                #print "events in chan", len(chan)
                    for e,event in enumerate(chan): #loop over 256 channels
                        if e==num:
                        #print "len check",len(event)
                            if not ped:
                                chan_list[num].extend(event) 
                            else: #return mean averaged events
                                meane=np.mean(event)
                                mean_evs=[]
                                for evi in event:
                                    mean_evs.append(evi-meane)
                                    #print evi
                                mean_list[num].append(meane)
                                chan_list[num].extend(mean_evs)
    print  "len chan_list[0]", len(chan_list[0])
    if not ped:
        return chan_list
    else:
        return chan_list, mean_list

def get_channel_rms(array, plot=False, returnarr=False, color='cornflowerblue'):
#plots histogram of the rms for each channel (all events and samples) in a pcb if plot is true and returns the median rms over all channels in that board
    chan_rms=[]
    chan_list=get_chan_list(array)
    for chan in chan_list:
        chan_rms.append(np.std(chan))#rms of all samples and events in that channel
    print "len(chan_rms)", len(chan_rms)
    medrms=np.median(chan_rms)#median rms for that BOARD not channel
    print "median channel rms", medrms 
    if plot:
        n, bins, p = plt.hist(chan_rms, 50, normed=0, facecolor=color,alpha=0.7)
        plt.axvline(medrms, color='k', linestyle='dashed', linewidth=2)
        plt.xlabel('RMS (std of all events and samples) for all channels', fontsize=20)
        plt.show()
    if returnarr:
        return medrms, chan_rms
    else:
        return medrms

def get_channel_ped(array,plotmean=False,plotped=False, color='cornflowerblue'):
#returns pedestal plot and median pedestal for a board in hg or low gain
#the mean is the mean over the event, for each sample
    chan_ped=[]; chan_mean=[]
    chan_list, mean_list=get_chan_list(array, ped=True)
    print "mean list len (13?)", len(mean_list[0])
    for avgofev in mean_list:
        chan_mean.extend(avgofev)
    samp_mean=np.mean(chan_mean)
    print "number of entries", len(chan_mean)
    if plotmean:
        n, bins, p = plt.hist(chan_mean, 50, normed=0, facecolor=color,alpha=0.7)
        plt.axvline(samp_mean, color='k', linestyle='dashed', linewidth=2)
        plt.xlabel('Mean over events for all samples and all channels', fontsize=20)
        plt.show()
    #plot the pedestal subtracted samples:
    for chan in chan_list:
        for val in chan:
            chan_ped.append(val)
    if plotped:
        n, bins, p = plt.hist(chan_ped, 50, normed=0, facecolor=color,alpha=0.7)#indianred
        plt.axvline(np.mean(chan_ped), color='k', linestyle='dashed', linewidth=2)
        plt.xlabel('Mean subtracted events for all samples and all channels', fontsize=20)
        plt.show()
    print "chan_ped entries", len(chan_ped)
    return samp_mean, chan_ped

def plot_pcb_med_rms(medrmslist, color='cornflowerblue'):
#plots the median rms given an array of median board rmses. finds the overall median rms for all boards
    unimed=np.median(medrmslist)
    n, bins, p = plt.hist(medrmslist, 50, normed=0, facecolor=color,alpha=0.7)
    plt.axvline(unimed, color='k', linestyle='dashed', linewidth=2)
    plt.xlabel('med rms of all channels in each pcb', fontsize=20)
    plt.show()
    return unimed

def plot_noise_dist(bigarray, unimed, color='cornflowerblue'):#still have missing events! 160 with rollnoise
#plots the noise distribution for all boards channels events and samples for pedestal subracted samples divided by universal rms
    for elem in bigarray:
        elem = elem/unimed
    #mu, sig = norm.fit(bigarray)

    n, bins, p = plt.hist(bigarray, 50, normed=0, facecolor=color,alpha=0.7)
    plt.axvline(np.mean(bigarray), color='k', linestyle='dashed', linewidth=2)
    plt.xlabel('noise values for all pcbs in terms sigmas of universal medians', fontsize=20)
    #y = mlab.normpdf(bins, mu, sig)
    #fit = plt.plot(bins, y, '--', linewidth=3)
    plt.show()
    
def num_noisy_chans(array, thr=3.0):#thr=[3.0, 4.0, 5.0]):#should make it take an array of thresh you want to test and return subplots
#thr is factor to multiply the median
#if channel rms is over some factor* board's median rms, that channel is noisy
    med, rms_arr = get_channel_rms(array,returnarr=True)
    print "thresholds are", thr
    #NCarr=[]
    #for t in thr:
    numnoisychans=0
    for r, rms in enumerate(rms_arr):
        if rms>thr*med:
            numnoisychans+=1
            print "found noisy channel at channel", r, "thr",thr,"times the median rms"
    return numnoisychans#, thr#NCarr

def plot_noisy_chans(ncarray, thr=3.0, color='cornflowerblue'):
#plots number of noisy channels for all boards
#should be subplots in the future
    nbins=len(ncarray)
    n, bins, p = plt.hist(ncarray, nbins+1, normed=0, facecolor=color,alpha=0.7)
    #plt.axvline(np.mean(bigarray), color='k', linestyle='dashed', linewidth=2)
    #plt.xlabel('number of noisy channels for all boards', fontsize=20)
    plt.title(thr)
    plt.show()

#######################################################################
if __name__ == "__main__":
    #LGMegaArr=np.empty((256,482560000)); HGMegaArr=np.empty((256,482560000))
    LGmedrms=[]; HGmedrms=[]
    LGnoise=[]; HGnoise=[]
    LGtotnc3=[]; HGtotnc3=[];LGtotnc4=[]; HGtotnc4=[];LGtotnc5=[]; HGtotnc5=[];#number of total noisy channels in each board

    if len(sys.argv) > 1:
        print "Number of files: ", len(sys.argv)
        print '# Input files are',sys.argv
        i=0
        for file in sys.argv:
            #LGarr=[];HGarr=[]           
            if i>0: 
                print "File # ", i
                #get_npz_noisevals(file)
                lgdata,hgdata=read_root_file(file)
                lgmedrms= get_channel_rms(lgdata,plot=True); hgmedrms=get_channel_rms(hgdata,plot=True,color='indianred')
                LGmedrms.append(lgmedrms); HGmedrms.append(hgmedrms);
                lgmean,lgnoise=get_channel_ped(lgdata, plotped=True); hgmean,hgnoise=get_channel_ped(hgdata, plotped=True,color='indianred');
                #LGnoise.extend(lgnoise); HGnoise.extend(hgnoise)
                LGtotnc3.append(num_noisy_chans(lgdata,thr=3.0)); LGtotnc4.append(num_noisy_chans(lgdata,thr=4.0)); LGtotnc5.append(num_noisy_chans(lgdata,thr=5.0))
                HGtotnc3.append(num_noisy_chans(hgdata,thr=3.0)); HGtotnc4.append(num_noisy_chans(hgdata,thr=4.0)); HGtotnc5.append(num_noisy_chans(hgdata,thr=5.0))
            i+=1
    else:
        print "No input files given!"
        exit(0)
    print "med len (30)", len(LGmedrms), len(HGmedrms)
    print "noise len (30*~100*13)", len(LGnoise),len(HGnoise)
    #plot general pcb stuff:
    #lg_unimed=plot_pcb_med_rms(LGmedrms); hg_unimed=plot_pcb_med_rms(HGmedrms,color='indianred');
    #plot_noise_dist(LGnoise,lg_unimed); plot_noise_dist(HGnoise,hg_unimed,color='indianred')
    plot_noisy_chans(LGtotnc3, thr=3.0); plot_noisy_chans(HGtotnc3, thr=3.0, color='indianred')
    plot_noisy_chans(LGtotnc4, thr=4.0); plot_noisy_chans(HGtotnc4, thr=4.0, color='indianred')
    plot_noisy_chans(LGtotnc5, thr=5.0); plot_noisy_chans(HGtotnc5, thr=5.0, color='indianred')
