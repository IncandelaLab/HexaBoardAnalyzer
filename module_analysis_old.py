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
import matplotlib
matplotlib.rc("figure", facecolor="white")
matplotlib.rcParams.update({'font.size': 22})
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
#chan_list has all samples and events in that channel
#mean_list has mean of all events in that sample in that channel
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

def get_channel_rms(array, plot=False, returnarr=False, gain='LG',Dir='./'):
#plots histogram of the rms for each channel (all events and samples) in a pcb if plot is true and returns the median rms over all channels in that board
    if gain=='HG'or gain=='hg': color='indianred'
    else: color='cornflowerblue'
    chan_rms=[]
    chan_list=get_chan_list(array)
    for chan in chan_list:
        chan_rms.append(np.std(chan))#rms of all samples and events in that channel
    print "len(chan_rms)", len(chan_rms)
    medrms=np.median(chan_rms)#median rms for that BOARD not channel
    print "median channel rms", medrms 
    if plot:
        plt.figure(figsize=(20,10))
        plt.hist(chan_rms, 50, normed=0, facecolor=color,alpha=0.7)
        plt.axvline(medrms, color='k', linestyle='dashed', linewidth=2)
        plt.xlabel('RMS (std of all events and samples) for all channels', fontsize=20)
        #plt.show()
        plt.savefig(Dir+gain+'board_rms_dist'+savetype)
    if returnarr:
        return medrms, chan_rms
    else:
        return medrms

def get_adc_dist(array, gain='LG',Dir='./'):
#returns distribution of the adc counts for all channels in a board
    if gain=='HG'or gain=='hg': color='indianred'
    else: color='cornflowerblue'
    adc_list=[]
    get_chans=get_chan_list(array)
    for chan in chan_list:
        adc_list.extend(chan)#literally just an array of all the adc counts in that board
    plt.figure(figsize=(20,10))
    n, bins, p = plt.hist(adc_list, 50, normed=0, facecolor=color,alpha=0.7)#indianred
    plt.axvline(np.mean(adc_list), color='k', linestyle='dashed', linewidth=2)
    plt.xlabel('ADC counts for all channels', fontsize=20)
    #plt.show()
    plt.savefig(Dir+gain+'board_adc_dist'+savetype)

def get_channel_ped(array,plotmean=False,plotped=False, gain='LG',Dir='./'):
#returns pedestal plot and median pedestal for a board in hg or low gain
#the mean is the mean over the event, for each sample
    if gain=='HG'or gain=='hg': color='indianred'
    else: color='cornflowerblue'
    chan_ped=[]; chan_mean=[]
    chan_list, mean_list=get_chan_list(array, ped=True)
    print "mean list len (13?)", len(mean_list[0])
    for avgofev in mean_list:
        chan_mean.extend(avgofev)
    samp_mean=np.mean(chan_mean)
    print "number of entries", len(chan_mean)
    if plotmean:
        plt.figure(figsize=(20,10))
        n, bins, p = plt.hist(chan_mean, 50, normed=0, facecolor=color,alpha=0.7)
        plt.axvline(samp_mean, color='k', linestyle='dashed', linewidth=2)
        plt.xlabel('Mean over events for all samples and all channels', fontsize=20)
        plt.savefig(Dir+gain+'plot_mean'+savetype)
        #plt.show()
    #plot the pedestal subtracted samples:
    for chan in chan_list:
        for val in chan:
            chan_ped.append(val)
    if plotped:
        plt.figure(figsize=(20,10))
        n, bins, p = plt.hist(chan_ped, 50, normed=0, facecolor=color,alpha=0.7)#indianred
        plt.axvline(np.mean(chan_ped), color='k', linestyle='dashed', linewidth=2)
        plt.xlabel('Mean subtracted events for all samples and all channels', fontsize=20)
        #plt.show()
        plt.savefig(Dir+gain+'meansub_ped'+savetype)
    print "chan_ped entries", len(chan_ped)
    return samp_mean, chan_ped

def plot_pcb_med_rms(medrmslist, gain='LG',Dir='./'): 
#plots the median rms given an array of median board rmses. finds the overall median rms for all boards
    if gain=='HG'or gain=='hg': color='indianred'
    else: color='cornflowerblue'
    unimed=np.median(medrmslist)
    plt.figure(figsize=(20,10))
    n, bins, p = plt.hist(medrmslist, 50, normed=0, facecolor=color,alpha=0.7)
    plt.axvline(unimed, color='k', linestyle='dashed', linewidth=2)
    plt.xlabel('med rms of all channels in each pcb', fontsize=20)
    #plt.show()
    plt.savefig(Dir+gain+'plot_pcb_med_rms'+savetype)   # save the figure to file
    # plt.close(fig) 
    return unimed

def board_rms_compare(before_array, after_array, showdist=True, gain='LG',Dir='./'):
#takes 2 data arrays and plots the difference in rms per channel for each channel in those arrays
#do subplots for active and inactive channels
    if gain=='HG'or gain=='hg': color='indianred'
    else: color='cornflowerblue'
    diff_rms=[]; sub_rms=[]
    chans=range(0, 256)
    bef_med,bef_rms = get_channel_rms(before_array, returnarr=True)
    aft_med,aft_rms = get_channel_rms(after_array, returnarr=True)
    for bef,aft in zip(bef_rms,aft_rms):
        diff_rms.append(aft/bef)
        sub_rms.append(aft-bef)
    print "len arr (256)", len(diff_rms)
    evn_diff_med=np.median(diff_rms[::2])#aft_med/bef_med
    odd_diff_med=np.median(diff_rms[1::2])
    width=.35
    print "diff arr", bef_rms[0:5],aft_rms[0:5], diff_rms[0:5]
    fig,ax=plt.subplots(2,1,figsize=(30,20))
    ax[0].bar(chans[::2], diff_rms[::2], width, color=color, alpha=0.7,label='RMS change')#r'$\Delta$'+'RMS'
    ax[0].axhline(evn_diff_med, color='k', linestyle='dashed', linewidth=2,label='med RMS change')
    ax[0].set_xlabel('Channel number')
    ax[0].set_ylabel('(new channel rms)/(old channel rms)')
    ax[0].set_title('Factor of RMS change (new/old) per Channel (Active Channels)')
    ax[0].legend(loc='lower right')
    ax[1].bar(chans[1::2], diff_rms[1::2], width, color=color, alpha=0.7,label=r'RMS change')
    ax[1].axhline(odd_diff_med, color='k', linestyle='dashed', linewidth=2,label='med RMS change')
    ax[1].set_xlabel('Channel number')
    ax[1].set_ylabel('(new channel rms)/(old channel rms)')
    ax[1].set_title('Factor of RMS change (new/old) per Channel (Inactive Channels)')
    ax[1].legend(loc='lower right')
    #plt.legend(bbox_to_anchor=(0, 1), loc='upper left', ncol=1)
    ##plt.show()
    plt.savefig(Dir+gain+'rms_chan_comp'+savetype)
    if showdist:
        plt.figure(figsize=(20,10))
        n, bins, p = plt.hist(sub_rms[::2], 80, normed=0, facecolor=color,alpha=0.7,label='active channels')
        plt.hist(sub_rms[1::2], bins, normed=0, facecolor='silver',alpha=0.4,label='inactive channels')
        plt.xlabel('(new channel rms)-(old channel rms)')
        plt.title('RMS Difference (new-old) Distribution')
        plt.legend(loc='upper right')
        #plt.show()
        plt.savefig(Dir+gain+'rms_diff_dist'+savetype)

def board_mean_adc_compare(before_array, after_array, showdist=True, gain='LG',Dir='./'):
#takes 2 data arrays and plots the difference in rms per channel for each channel in those arrays
    if gain=='HG'or gain=='hg': color='indianred'
    else: color='cornflowerblue'
    diff_mean=[]; sub_mean=[]
    chans=range(0, 256)
    bef_list= get_chan_list(before_array)
    aft_list= get_chan_list(after_array)
    for aftchan, befchan in zip(aft_list,bef_list):
        bef_chan_mean=np.mean(befchan)#take mean of all events and samples in that channel
        aft_chan_mean=np.mean(aftchan)
        diff_mean.append(aft_chan_mean/bef_chan_mean)
        sub_mean.append(aft_chan_mean-bef_chan_mean)
    print "len arr (256)", len(diff_mean)
    evn_diff_totmean=np.mean(diff_mean[::2])#mean of means
    odd_diff_totmean=np.mean(diff_mean[1::2])#mean of means
    width=.35
    fig,ax=plt.subplots(2,1,figsize=(30,20))
    ax[0].bar(chans[::2], diff_mean[::2], width, color=color, alpha=0.7,label='mean adc change')
    ax[0].axhline(evn_diff_totmean, color='k', linestyle='dashed', linewidth=2,label='overall mean adc change')
    ax[0].set_xlabel('Channel number')
    ax[0].set_ylabel('(new channel mean adc)/(old channel mean adc)')
    ax[0].set_title('Factor of Mean ADC change (new/old) per Channel (Active Channels)')
    ax[0].legend(loc='lower right')
    ax[1].bar(chans[1::2], diff_mean[1::2], width, color=color, alpha=0.7,label='mean adc change')
    ax[1].axhline(odd_diff_totmean, color='k', linestyle='dashed', linewidth=2,label='overall mean adc change')
    ax[1].set_xlabel('Channel number')
    ax[1].set_ylabel('(new channel mean adc)/(old channel mean adc)')
    ax[1].set_title('Factor of Mean ADC change (new/old) per Channel (Inactive Channels)')
    ax[1].legend(loc='lower right')
    ##plt.show()
    plt.savefig(Dir+gain+'adc_chan_comp'+savetype)
    if showdist:
        plt.figure(figsize=(20,10))
        n, bins, p = plt.hist(sub_mean[::2], 80, normed=0, facecolor=color,alpha=0.7,label='active channels')
        plt.hist(sub_mean[1::2], bins, normed=0, facecolor='silver',alpha=0.4,label='inactive channels')
        plt.xlabel('(new channel mean adc)-(old channel mean adc)')
        plt.title('Mean ADC Difference (new-old) Distribution')
        plt.legend(loc='upper right')
        plt.savefig(Dir+gain+'adc_diff_dist'+savetype)

def plot_noise_dist(bigarray, unimed, gain='LG',Dir='./'):#still have missing events! 160 with rollnoise
#plots the noise distribution for all boards channels events and samples for pedestal subracted samples divided by universal rms
    if gain=='HG'or gain=='hg': color='indianred'
    else: color='cornflowerblue'
    for elem in bigarray:
        elem = elem/unimed
    n, bins, p = plt.hist(bigarray, 50, normed=0, facecolor=color,alpha=0.7)
    plt.axvline(np.mean(bigarray), color='k', linestyle='dashed', linewidth=2)
    plt.xlabel('noise values for all pcbs in terms sigmas of universal medians', fontsize=20)
    ##plt.show()
    plt.savefig(Dir+gain+'allpcb_noise_dist'+savetype)

def num_noisy_chans(array, thr=[3.0,4.0,5.0], gain='LG', Dir='./', pltboardnc=False):
#gets med-rms based noisy channel info for board and plots as # per thr and prints what channels were noisy
#thr is factor to multiply the median
#if channel rms is over some factor* board's median rms, that channel is noisy
    NCarr=[]; chan_label=[];
    if gain=='HG'or gain=='hg': color='indianred'
    else: color='cornflowerblue'
    med, rms_arr = get_channel_rms(array,returnarr=True)
    print "thresholds are", thr
    for i,t in enumerate(thr):
        thrlist=[]
        numnoisychans=0
        chantag=''
        for r, rms in enumerate(rms_arr):
            if rms>t*med:
                numnoisychans+=1
                chantag+=(str(r)+',')
                print "found noisy channel at channel", r, "thr",thr,"times the median rms"
        NCarr.append(numnoisychans)# should be [#nc thr1, #nc thr2,...]
        if chantag=='': chantag='none'
        if chantag.endswith(','): chantag=chantag[:-1]#strip comma at end
        chan_label.append(chantag)
    print 'chan_label',chan_label
    if pltboardnc:
        plt.figure(figsize=(20,10))
        mybar=plt.bar(thr, NCarr, .5, color=color, alpha=0.7)#r'$\Delta$'+'RMS'
        plt.xlabel('Threshold (factor x median channel rms)')
        plt.ylabel('# Noisy Channels (with channel number labels)')
        plt.title('PCB # Noisy Channels per threshold (with noisy channel numbers)')
        #put channels that were noisy in test above the bar
        for j,rect in enumerate(mybar):
            height = rect.get_height()
            plt.text(rect.get_x() + rect.get_width()/2.0, height, '%s' % chan_label[j], ha='center', va='bottom')
        ##plt.show()
        plt.savefig(Dir+gain+'board_nc_per_thr'+savetype)
    return NCarr, thr

def plot_all_noisy_chans(ncarray, thr, gain='LG',Dir='./'):
#plots number of noisy channels for all boards as distribution, as # per pcb
#ncarray has rows per thr and cols per board
    if gain=='HG'or gain=='hg': color='indianred'
    else: color='cornflowerblue'
    #histogram of # pcbs with certain # of noisy channels per thr
    fig,ax=plt.subplots(1,len(thr),figsize=(30,20))
    i=0
    for axi,t in zip(ax,thr):
        nbins=len(ncarray[i,:])
        print 'nbins',nbins
        axi.hist(ncarray[i,:], nbins+1, normed=0, facecolor=color,alpha=0.7)
        axi.set_title('Number of Noisy Channels for all PCBS at thr '+str(t))
        axi.set_ylabel('Number of pcbs with that # of noisy channels')
        axi.set_xlabel('Number of noisy channels')
        i+=1
    #plt.show()
    plt.savefig(Dir+gain+'allpcb_noisy_chans_dist'+savetype)
    #bar graph of # nc per pcb
    fig,ax=plt.subplots(len(thr),1,figsize=(25,20))
    pcbs=range(0,len(ncarray[0,:]))#list of pcbs 1-30 or whatever
    i=0
    for axi,t in zip(ax,thr):
        axi.bar(pcbs, ncarray[i,:], 0.4, color=color,alpha=0.7)
        axi.set_title('Number of Noisy Channels per PCB at thr '+str(t))
        axi.set_ylabel('Number of noisy channels')
        axi.set_xlabel('PCB #')
        i+=1
    #plt.show()
    plt.savefig(Dir+gain+'noisy_chans_per_pcb'+savetype)

def FileReader():#dont need to put in sys.argv
#reads the files given accordingly, returns list of files and whether it is for comparing module stages or not
    compare=False
    filelist=[]
    if len(sys.argv)==1:
        print "No files given!"
        exit(0)
    elif len(sys.argv) > 1:
        if len(sys.argv)==3:
            compare=True
            print "2 files found: Treating as assembly stage comparison."
        else: print "Analyzing",len(sys.argv)-1, " files"
        filelist=sys.argv[1:]
        #find out whether you need to create the rootfiles or not
        if any(".txt" in f for f in filelist):
            print "raw data given. creating rootfiles."
            for fname in filelist:
                if ".root" not in fname:
                    #do things in analyze hexanoise, returning list with root files instead of txt files
                    rt.gROOT.SetBatch(True)
                    fname = createTree(fname)
        elif any(".txt" not in f for f in filelist):
            print "rootfiles given"
    return filelist, compare

def AssemblyStageCompare(twofile_list):
#takes 2 root files for before and after assembly and spits out adc/rms comparison plots for lg and hg
#saves in folder with BOTH run names
    f1 = twofile_list[0].replace('.root',''); f2 = twofile_list[1].replace('.root','')
    comppltdir =  f1+'_'+f2+'_Comp_'+ '_plots/'
    if not os.path.exists(comppltdir): os.makedirs(comppltdir)
    print("Output dir: " + comppltdir)

    lgdatalast,hgdatalast=read_root_file(twofile_list[0]); lgdata,hgdata=read_root_file(twofile_list[1])
    board_rms_compare(lgdatalast,lgdata,Dir=comppltdir); board_rms_compare(hgdatalast,hgdata,gain='HG',Dir=comppltdir)
    board_mean_adc_compare(lgdatalast,lgdata,Dir=comppltdir); board_mean_adc_compare(hgdatalast,hgdata,gain='HG',Dir=comppltdir)

def PCBAnalyzer(file_list, multipcb=False, withPlotNoise=True):
#plots and saves everything, including plotNoise graphs and creating rootfiles if needed, in the folder of origin (where plotNoise saves things anyway)
#for functions that use full list of pcbs: noisy channels, noise of all pcbs, etc enable multipcb
    LGmedrms=[]; HGmedrms=[]; LGnoise=[]; HGnoise=[]; LGtotNC=[]; HGtotNC=[]; LGncarr=[]; HGncarr=[]
    for i, filename in enumerate(file_list):
        print "File # ", i
        if withPlotNoise:
            runPlotNoise(filename)
        pltdir = filename.replace('.root',''); pltdir = pltdir + '_plots/'
        if not os.path.exists(pltdir): os.makedirs(pltdir)
        print("Output dir: " + pltdir)

        lgdata,hgdata = read_root_file(filename)
        lgmedrms= get_channel_rms(lgdata,plot=True,Dir=pltdir); hgmedrms=get_channel_rms(hgdata,plot=True,gain='HG',Dir=pltdir)
        lgmean,lgnoise=get_channel_ped(lgdata, plotped=True,Dir=pltdir); hgmean,hgnoise=get_channel_ped(hgdata, plotped=True,gain='HG',Dir=pltdir);
        lgncs,lgThr=num_noisy_chans(lgdata,pltboardnc=(not multipcb),Dir=pltdir); hgncs,hgThr=num_noisy_chans(lgdata,gain='HG',pltboardnc=(not multipcb),Dir=pltdir)
        if multipcb:
            LGmedrms.append(lgmedrms); HGmedrms.append(hgmedrms); LGnoise.extend(lgnoise); HGnoise.extend(hgnoise)
            if i==0:
                LGncarr=lgncs; HGncarr=hgncs
            else:
                print 'LGncarr',LGncarr
                LGncarr=np.vstack((LGncarr,lgncs)); HGncarr=np.vstack((HGncarr,hgncs))

    if multipcb:#MAKE SURE ALL PLOTS PUT THINGS IN RIGHT DIRECTORY
        LGncarr=np.transpose(LGncarr); HGncarr=np.transpose(HGncarr)
        multipltdir='full_'+str(len(file_list))+'_pcb_plots/'
        if not os.path.exists(multipltdir): os.makedirs(multipltdir)
        print "Output dir: ", multipltdir
        lg_unimed=plot_pcb_med_rms(LGmedrms,Dir=multipltdir); hg_unimed=plot_pcb_med_rms(HGmedrms,gain='HG',Dir=multipltdir);
        plot_noise_dist(LGnoise,lg_unimed,Dir=multipltdir); plot_noise_dist(HGnoise,hg_unimed,gain='HG',Dir=multipltdir)
        plot_all_noisy_chans(LGncarr, lgThr,Dir=multipltdir); plot_all_noisy_chans(HGncarr, hgThr, gain='HG',Dir=multipltdir)

#######################################################################
if __name__ == "__main__":
    #pltdir='/' #note this variable is REDEFINED not PASSED as each function does its job
    savetype='.png'
    file_list, comparebool= FileReader()
    if comparebool:
        AssemblyStageCompare(file_list)
    elif not comparebool:
        if len(file_list)==1:
            print "analyzing one pcb..."
            PCBAnalyzer(file_list, multipcb=False)
        if len(file_list)>1:
            "analyzing set of ", len(file_list),"pcbs..."
            PCBAnalyzer(file_list, multipcb=True)
