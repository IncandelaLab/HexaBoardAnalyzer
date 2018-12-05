#!/usr/bin/env python2
import glob, sys, os, array, math
import numpy as np
from scipy.stats import norm
from scipy import signal
import string
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from chip_median_rms import *
plt.style.use('bmh')


#########################################################
if __name__ == "__main__":
    lgMedArr=[]; hgMedArr=[]
    thr2mlg=[]; thr3mlg=[]; thr4mlg=[]; thr5mlg=[]
    thr2mhg=[]; thr3mhg=[]; thr4mhg=[]; thr5mhg=[]
    thr2mlg_pcb=[]; thr3mlg_pcb=[]; thr4mlg_pcb=[]; thr5mlg_pcb=[]
    thr2mhg_pcb=[]; thr3mhg_pcb=[]; thr4mhg_pcb=[]; thr5mhg_pcb=[]

    lgMedArr= [5.52391721697 , 5.74734575189 , 5.87303901589 ,  6.09715088671
, 5.72736953215 , 5.93391584053 , 5.60980550746 ,  5.98104233799
, 5.89807005104 , 6.10599646155 , 5.92877431897 ,  5.79153547118
, 5.88043870137 , 6.16185060079 , 6.06264258486 ,  5.81585781998
, 5.91606174754 , 6.22296586638 , 5.75897247076 ,  5.94959737686
, 5.76348196644 , 5.89567639852 , 6.10492007063 ,  5.97798456884
, 5.74239123776 , 5.74037113689 , 6.01821901746 ,  5.79533030903
, 6.29807291602 , 6.21966806906 , 6.1981337748 ,  5.91798049263
, 5.80367513129 , 5.89172089115 , 5.99576129904 ,  5.6577616593
, 5.92584709025 , 5.76756436528 , 5.7200183507 ,  5.63890666849
, 5.94610747014 , 5.86695084327 , 6.07343823154 ,  5.75622101486
, 5.73258562816 , 5.7037581074 , 5.83116105986 ,  5.87571461243
, 5.72237942308 , 6.08687171218 , 5.99254817366 ,  5.89666724109
, 6.31325122596 , 5.81271102151 , 6.21073201291 ,  5.68660325312
, 6.11032141381 , 6.13116355409 , 6.18663989887 ,  5.45393327157
, 5.25514655869 , 5.67652866208 , 5.55894262659 ,  5.83448150509
, 5.63354449576 , 6.18803517867 , 6.10820939371 ,  6.09850581864
, 5.992647173 , 5.76596500546 , 6.10027085779 ,  6.31037543495
, 6.09461993386 , 6.07337362505 , 5.99273648699 ,  6.31136008592
, 6.05795820922 , 5.99072134752 , 6.2757935243 ,  5.82310973854
, 5.89581920588 , 5.73901888991 , 5.87826063115 ,  5.80162956958
, 5.98946015011 , 5.55748974392 , 5.70079453304 ,  5.44989596498
, 5.70977243318 , 6.08785429662 , 5.43618638359 ,  5.84899163183
, 5.69130379658 , 5.930020298 , 5.78753604884 ,  5.60488699586
, 5.5997092866 , 5.44369255427 , 5.40137695737 ,  5.59757649239
, 5.50365468911 , 5.65126290684 , 5.83237710869 ,  5.66511040667
, 5.66339105662 , 5.66521730274 , 5.90363032207 ,  5.97834851626
, 5.64508866634 , 5.93127790574 , 5.94631549608 ,  5.51786738043
, 5.79890600072 , 5.35078346766 , 5.75935836347 ,  5.77180315182]

    hgMedArr = [5.08514566554 , 5.418615047 , 5.6809600769 ,  5.92708731299
, 5.28574175907 , 5.52295691103 , 5.50096886371 ,  5.58759943359
, 5.42514505879 , 5.89644791392 , 5.65138348741 ,  5.41124646641
, 5.29141996498 , 5.68563119608 , 5.87799295554 ,  5.40889561183
, 5.51419888877 , 5.917464152 , 5.54506836738 ,  5.59748760765
, 5.25693727532 , 5.57746817602 , 5.8420151879 ,  5.58840848446
, 5.27976197208 , 5.74829507537 , 6.13825931482 ,  5.79866799877
, 5.65716094463 , 5.97629738365 , 6.20568889788 ,  5.67859088563
, 5.14234722591 , 5.54588157353 , 5.84455891332 ,  5.3633599174
, 5.38546200585 , 5.41317303357 , 5.39358526895 ,  5.1638184392
, 5.37819720224 , 5.77375838361 , 5.91345352026 ,  5.40410661794
, 5.305210213 , 5.42313220535 , 5.66284387342 ,  5.48781868891
, 5.13268681047 , 6.00548769357 , 5.82911536427 ,  5.66892963206
, 5.71593935385 , 5.66778164555 , 6.10402296629 ,  5.37773897421
, 5.39545741142 , 5.90567878737 , 5.90101720335 ,  5.01596491235
, 4.84097124582 , 5.3530542913 , 5.49667186008 ,  5.37058405929
, 5.08975425644 , 5.79374214924 , 5.89572982621 ,  5.48959210929
, 5.5633572123 , 5.63349725643 , 5.8938817495 ,  5.84431801366
, 5.68477477397 , 5.93585548461 , 5.75692063647 ,  5.88112747818
, 5.45837825127 , 5.67249005631 , 6.12968474776 ,  5.45646235334
, 5.29639439171 , 5.44479487311 , 5.62474600165 ,  5.26524073921
, 5.5247009123 , 5.26811346702 , 5.6046048796 ,  5.1541130325
, 5.21315169857 , 5.97123922142 , 5.51695647729 ,  5.41028703896
, 5.22595236797 , 5.79886467705 , 5.85330359146 ,  5.33554856647
, 5.08428121466 , 5.44151249269 , 5.48769621228 ,  5.37012368658
, 5.18838757307 , 5.60941716108 , 5.78017002259 ,  5.19816388295
, 5.25327381013 , 5.63733411051 , 5.94591292565 ,  5.54943036727
, 5.25677370413 , 5.75285594818 , 5.79739474033 ,  5.26087167531
, 5.37123913108 , 5.34191847666 , 5.74831754928 , 5.46697161316]

    boardplot=True
    print "Number of files: ", len(sys.argv)
    if len(sys.argv) > 1:
        print '# Input files are',sys.argv
        i=0
        for file in sys.argv:
            if i!=0: #there is a problem with the data for C11
                #lg1,lg2,lg3,lg4,hg1,hg2,hg3,hg4=chip_median_rms(file,medplot=True)
                #lgMedArr.append(lg1); lgMedArr.append(lg2); lgMedArr.append(lg3); lgMedArr.append(lg4)
                #hgMedArr.append(hg1); hgMedArr.append(hg2); hgMedArr.append(hg3); hgMedArr.append(hg4)
                lgchans,hgchans=chip_median_rms(file)
                lgpcb2m=0;lgpcb3m=0;lgpcb4m=0;lgpcb5m=0;hgpcb2m=0;hgpcb3m=0;hgpcb4m=0;hgpcb5m=0
                j=0
                while j<4:
                    thr2mlg.append(lgchans[0][j]); thr3mlg.append(lgchans[1][j]); thr4mlg.append(lgchans[2][j]); thr5mlg.append(lgchans[3][j])
                    thr2mhg.append(hgchans[0][j]); thr3mhg.append(hgchans[1][j]); thr4mhg.append(hgchans[2][j]); thr5mhg.append(hgchans[3][j])
                    if boardplot:
                        lgpcb2m+=lgchans[0][j];lgpcb3m+=lgchans[1][j];lgpcb4m+=lgchans[2][j];lgpcb5m+=lgchans[3][j]
                        hgpcb2m+=hgchans[0][j];hgpcb3m+=hgchans[1][j];hgpcb4m+=hgchans[2][j];hgpcb5m+=hgchans[3][j]
                    j+=1
                if lgpcb2m>0:
                    thr2mlg_pcb.append(lgpcb2m)
                if lgpcb3m>0:
                    thr3mlg_pcb.append(lgpcb3m)
                if lgpcb4m>0:
                    thr4mlg_pcb.append(lgpcb4m)
                if lgpcb5m>0:
                    thr5mlg_pcb.append(lgpcb5m)
                if hgpcb2m>0:
                    thr2mhg_pcb.append(hgpcb2m)
                if hgpcb3m>0:
                    thr3mhg_pcb.append(hgpcb3m)
                if hgpcb4m>0:
                    thr4mhg_pcb.append(hgpcb4m)
                if hgpcb5m>0:
                    thr5mhg_pcb.append(hgpcb5m)
                print "thr2mlg", thr2mlg, lgchans[0][:]
                print "thr2mhg", thr2mhg, hgchans[0][:]
            i+=1
    else:
        print "No input files given!"
        #exit(0)
        #arrays for all but C11

    fig, ax = plt.subplots(1,2)
    nl, binsl, pl = ax[0].hist(lgMedArr, 15, normed=0, facecolor='cornflowerblue',alpha=0.7)
    nh, binsh, ph = ax[1].hist(hgMedArr, 15, normed=0, facecolor='indianred', alpha=0.7)

    LG_med=np.median(lgMedArr)
    HG_med=np.median(hgMedArr)
    ax[0].axvline(LG_med, color='k', linestyle='dashed', linewidth=2)
    ax[1].axvline(HG_med, color='k', linestyle='dashed', linewidth=2)
    
    print "lgmed", LG_med , "hgmed", HG_med
    print "number entries", len(lgMedArr), len(hgMedArr)
    print "number entries", len(lgMedArr), len(hgMedArr)
    ax[0].set_xlabel('LG median rms all chips all pcbs', fontsize=20)
    ax[1].set_xlabel('HG median rms all chips all pcbs', fontsize=20)
    plt.show()
    plt.close()

    # lgdim = lgchans.size[0]
    # hgdim = hgchans.size[0]
    # print lgdim

    #should loop and create number of arrays to match number of factors but whatever. add it later
    #for lc,hc in zip(lgchans, hgchans):
    print "#entries", thr2mlg_pcb
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8) ) = plt.subplots(4,2)
    nl, binsl, pl = ax1.hist(thr2mlg, 10, normed=0, facecolor='cornflowerblue',alpha=0.7)
    nl, binsl, pl = ax3.hist(thr3mlg, 10, normed=0, facecolor='cornflowerblue',alpha=0.7)
    nl, binsl, pl = ax5.hist(thr4mlg, 10, normed=0, facecolor='cornflowerblue',alpha=0.7)
    nl, binsl, pl = ax7.hist(thr5mlg, 10, normed=0, facecolor='cornflowerblue',alpha=0.7)
    nh, binsh, ph = ax2.hist(thr2mhg, 10, normed=0, facecolor='indianred', alpha=0.7)
    nh, binsh, ph = ax4.hist(thr3mhg, 10, normed=0, facecolor='indianred', alpha=0.7)
    nh, binsh, ph = ax6.hist(thr4mhg, 10, normed=0, facecolor='indianred', alpha=0.7)
    nh, binsh, ph = ax8.hist(thr5mhg, 10, normed=0, facecolor='indianred', alpha=0.7)
    ax1.set_xlabel('LG # noisy channels thr 2*med all chips all pcbs', fontsize=20)
    ax3.set_xlabel('LG # noisy channels thr 3*med all chips all pcbs', fontsize=20)
    ax5.set_xlabel('LG # noisy channels thr 4*med all chips all pcbs', fontsize=20)
    ax7.set_xlabel('LG # noisy channels thr 5*med all chips all pcbs', fontsize=20)
    ax2.set_xlabel('HG # noisy channels thr 2*med all chips all pcbs', fontsize=20)
    ax4.set_xlabel('HG # noisy channels thr 3*med all chips all pcbs', fontsize=20)
    ax6.set_xlabel('HG # noisy channels thr 4*med all chips all pcbs', fontsize=20)
    ax8.set_xlabel('HG # noisy channels thr 5*med all chips all pcbs', fontsize=20)
    plt.show()

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8) ) = plt.subplots(4,2)
    nl, binsl, pl = ax1.hist(thr2mlg_pcb, 10, normed=0, facecolor='cornflowerblue',alpha=0.7)
    nl, binsl, pl = ax3.hist(thr3mlg_pcb, 10, normed=0, facecolor='cornflowerblue',alpha=0.7)
    nl, binsl, pl = ax5.hist(thr4mlg_pcb, 10, normed=0, facecolor='cornflowerblue',alpha=0.7)
    nl, binsl, pl = ax7.hist(thr5mlg_pcb, 10, normed=0, facecolor='cornflowerblue',alpha=0.7)
    nh, binsh, ph = ax2.hist(thr2mhg_pcb, 10, normed=0, facecolor='indianred', alpha=0.7)
    nh, binsh, ph = ax4.hist(thr3mhg_pcb, 10, normed=0, facecolor='indianred', alpha=0.7)
    nh, binsh, ph = ax6.hist(thr4mhg_pcb, 10, normed=0, facecolor='indianred', alpha=0.7)
    nh, binsh, ph = ax8.hist(thr5mhg_pcb, 10, normed=0, facecolor='indianred', alpha=0.7)
    # hist, bins = np.histogram(gaussian_numbers, bins=12)
    # hist[np.where(hist <= freq)] = 0
    ax1.set_xlabel('LG # noisy channels thr 2*med all chips all pcbs', fontsize=20)
    ax3.set_xlabel('LG # noisy channels thr 3*med all chips all pcbs', fontsize=20)
    ax5.set_xlabel('LG # noisy channels thr 4*med all chips all pcbs', fontsize=20)
    ax7.set_xlabel('LG # noisy channels thr 5*med all chips all pcbs', fontsize=20)
    ax2.set_xlabel('HG # noisy channels thr 2*med all chips all pcbs', fontsize=20)
    ax4.set_xlabel('HG # noisy channels thr 3*med all chips all pcbs', fontsize=20)
    ax6.set_xlabel('HG # noisy channels thr 4*med all chips all pcbs', fontsize=20)
    ax8.set_xlabel('HG # noisy channels thr 5*med all chips all pcbs', fontsize=20)
    plt.show()

    #plt.close()

    
