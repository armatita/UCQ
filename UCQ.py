# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 15:06:21 2014

@author: pedro.correia
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def do_pearson_correlation(a,b):
    return np.sum((a-a.mean())*(b-b.mean()))/np.sqrt(np.sum((a-a.mean())**2)*np.sum((b-b.mean())**2))
    
def do_reflective_similarity(a,b):
    return np.sum(a*b)/np.sqrt(np.sum(a**2)*np.sum(b**2))
    
def do_quasi_correlation(a,b):
    return 2*np.sum(a*b)/(np.sum(a**2)+np.sum(b**2))
    
def do_squared_quasi_correlation(a,b):
    return (2*np.sum(a*b)/(np.sum(a**2)+np.sum(b**2)))**2

def do_spearman_correlation(a,b):
    return np.sum((a-a.mean())*(b-b.mean()))/np.sqrt(np.sum((a-a.mean())**2)*np.sum((b-b.mean())**2))
    
def do_cosine_correlation(a,b):
    return np.sum(a*b)/np.sqrt(np.sum(a**2)*np.sum(b**2))

def convolve_trace(trace,wavelet):
    reflectivity = trace[:,-1].copy()
    for i in xrange(trace.shape[0]-1):
        reflectivity[i]=(reflectivity[i+1]-reflectivity[i])/(reflectivity[i+1]+reflectivity[i])
    reflectivity[-1]=0 #grid[:,:,-2]
    synthetic = trace[:,-1].copy()
    synthetic[:] = 0
    h_size=(wavelet.shape[0]-1)/2
    for i in xrange(trace.shape[0]):
        if i-h_size<0:
            wa=h_size-i
            a=0
        else:
            wa=0
            a=i-h_size
        if i+h_size>trace.shape[0]:
            wb=h_size+i-trace.shape[0]
            b=trace.shape[0]
        else:
            wb=2*h_size+1
            b=i+h_size
        #print synthetic[a:b].shape,wavelet[wa:(2*h_size-wb)].shape
        synthetic[a:b]=synthetic[a:b]+reflectivity[i]*wavelet[wa:(2*h_size-wb)]
    return synthetic
    
def initial_stats():
    ua = np.loadtxt('logs_AI_1.prn')
    ub = np.loadtxt('logs_AI_2.prn')
    s = np.load('seismic.npy')
    ta = s[np.int_(ua[:,0]),np.int_(ua[:,1]),np.int_(ua[:,2])]
    tb = s[np.int_(ub[:,0]),np.int_(ub[:,1]),np.int_(ub[:,2])]
    w  = np.loadtxt('wavelet.txt')[:,1]
    ca = convolve_trace(ua,w)
    cb = convolve_trace(ub,w)
    print 'Old well A quasi-correlation: ',do_quasi_correlation(ca,ta)
    print 'Old well B quasi-correlation: ',do_quasi_correlation(cb,tb)
    return w,ta,tb
    
def load_las():
    apath = 'abc.las'
    bpath = 'def.las'
    
    a = np.loadtxt(apath,skiprows=44)
    b = np.loadtxt(bpath,skiprows=45)
    
    # A 1:~A DEPTH 2:TVD 3:TVDSS 4:DEVX 5:DEVY 6:TIME 7:PHIE 8:DTS 9:IS 10:DT 11:AI 12:GR 13:RHOB 
    # B 1:~A DEPTH 2:TVD 3:TVDSS 4:DEVX 5:DEVY 6:TIME 7:PHIE 8:DTS 9:DT 10:AI 11:GR 12:RHOB 13:IS
    
    awell = a[np.where(a[:,10]!=-999.25)[0],:] #,[1,3,4,10]]
    bwell = b[np.where(b[:,9]!=-999.25)[0],:] #,[1,3,4,9]]
    awell = awell[:,[3,4,5,10]]
    bwell = bwell[:,[3,4,5,9]]
    awell = awell[np.where(awell[:,2]!=-999.25),:]
    bwell = bwell[np.where(bwell[:,2]!=-999.25),:]
    awell = awell[0,:,:]
    bwell = bwell[0,:,:]
    
    ha = np.histogram(awell[:,2],bins=np.arange(awell[:,2].min(),awell[:,2].max()+4,4))
    hb = np.histogram(bwell[:,2],bins=np.arange(bwell[:,2].min(),bwell[:,2].max()+4,4))
    h1a = ha[1][213]
    h2a = ha[1][321]
    h1b = hb[1][194]
    xa = (awell[np.where((awell[:,2]>h1a) & (awell[:,2]<h2a))[0],2]-awell[np.where((awell[:,2]>h1a) & (awell[:,2]<h2a))[0],2].min())/4
    ya = awell[np.where((awell[:,2]>h1a) & (awell[:,2]<h2a))[0],3]
    xb = (bwell[np.where((bwell[:,2]>h1b))[0],2]-bwell[np.where((bwell[:,2]>h1b))[0],2].min())/4
    yb = bwell[np.where((bwell[:,2]>h1b))[0],3]
    return xa,ya,xb,yb
    
    
def load_stuff():
    apath = 'abc.las'
    bpath = 'def.las'
    
    a = np.loadtxt(apath,skiprows=44)
    b = np.loadtxt(bpath,skiprows=45)
    
    # A 1:~A DEPTH 2:TVD 3:TVDSS 4:DEVX 5:DEVY 6:TIME 7:PHIE 8:DTS 9:IS 10:DT 11:AI 12:GR 13:RHOB 
    # B 1:~A DEPTH 2:TVD 3:TVDSS 4:DEVX 5:DEVY 6:TIME 7:PHIE 8:DTS 9:DT 10:AI 11:GR 12:RHOB 13:IS
    
    awell = a[np.where(a[:,10]!=-999.25)[0],:] #,[1,3,4,10]]
    bwell = b[np.where(b[:,9]!=-999.25)[0],:] #,[1,3,4,9]]
    awell = awell[:,[3,4,5,10]]
    bwell = bwell[:,[3,4,5,9]]
    awell = awell[np.where(awell[:,2]!=-999.25),:]
    bwell = bwell[np.where(bwell[:,2]!=-999.25),:]
    awell = awell[0,:,:]
    bwell = bwell[0,:,:]
    
    ha = np.histogram(awell[:,2],bins=np.arange(awell[:,2].min(),awell[:,2].max()+4,4))
    hb = np.histogram(bwell[:,2],bins=np.arange(bwell[:,2].min(),bwell[:,2].max()+4,4))
    
    ca = 0
    for i in xrange(ha[1].shape[0]):
        ca = ca + 1
    cb = 0
    for i in xrange(hb[1].shape[0]):
        cb = cb + 1
    
    h1a = ha[1][213:321]
    h1b = hb[1][194::]
    return awell,bwell,h1a,h1b
    
def iterate(a,b,ha,hb,w,ta,tb,ites=1000):
    #na = a[:,:]
    #nb = b[:,:]
    fa = np.zeros((ha.shape[0]-1,2))
    fb = np.zeros((hb.shape[0]-1,2))
    for i in xrange(ha.shape[0]-1):
        ind = np.where((a[:,2]>=ha[i]) & (a[:,2]<ha[i+1]))
        fa[i,0] = a[ind[0],2].mean()
        p = np.random.randint(0,ind[0].shape[0]-1)
        fa[i,1] = a[ind[0][p],3]
    for i in xrange(hb.shape[0]-1):
        ind = np.where((b[:,2]>=hb[i]) & (b[:,2]<hb[i+1]))
        fb[i,0] = b[ind[0],2].mean()
        p = np.random.randint(0,ind[0].shape[0]-1)
        fb[i,1] = b[ind[0][p],3]
    na = fa[:,:]
    nb = fb[:,:]
    #print fb
    ca = convolve_trace(fa,w)
    cb = convolve_trace(fb,w)
    qa = np.array(do_quasi_correlation(ca,ta))
    qb = np.array(do_quasi_correlation(cb,tb))
    qqa = np.copy(qa)
    qqb = np.copy(qb)
    na = fa[:,:]
    nb = fb[:,:]
    print 'START well A quasi-correlation: ',qa
    print 'START well B quasi-correlation: ',qb
    #"""
    for i in xrange(ites):
        sa = np.random.randint(0,ha.shape[0]-1)
        ind_a = np.where((a[:,2]>=ha[sa]) & (a[:,2]<ha[sa+1]))
        pa = np.random.randint(0,ind_a[0].shape[0])
        appex = fa[sa,1]
        #print ind[0],pa,ind[0].shape[0]
        fa[sa,1] = a[ind_a[0][pa],3]
        ca = convolve_trace(fa,w)
        q2a = do_quasi_correlation(ca,ta)
        print 'ITE ',i,' CORRELATION FOR A IS: ',q2a
        if q2a >qa:
            qa = np.copy(q2a)
        else:
            fa[sa,1] = appex
    #"""
    for i in xrange(ites):
        sb = np.random.randint(0,hb.shape[0]-1)
        #print sb
        ind_b = np.where((b[:,2]>=hb[sb]) & (b[:,2]<hb[sb+1]))
        pb = np.random.randint(0,ind_b[0].shape[0])
        #print pb,ind_b[0].shape[0]
        appex = fb[sb,1]
        #print ind[0],pa,ind[0].shape[0]
        fb[sb,1] = b[ind_b[0][pb],3]
        cb = convolve_trace(fb,w)
        q2b = do_quasi_correlation(cb,tb)
        print 'ITE ',i,' CORRELATION FOR B IS: ',q2b
        if q2b >qb:
            qb = np.copy(q2b)
        else:
            fb[sb,1] = appex
    #"""
    #print na.shape
    initial_stats()
    print '#####################################'
    print 'START well A quasi-correlation: ',qqa
    print 'END well A quasi-correlation: ',qa
    print '#####################################'
    print 'START well B quasi-correlation: ',qqb
    print 'END well B quasi-correlation: ',qb
    xa,ya,xb,yb = load_las()
    ua = np.loadtxt('logs_AI_1.prn')
    ub = np.loadtxt('logs_AI_2.prn')
    #"""
    plt.plot(ua[:,2],ua[:,3],color='black',linewidth=3,label='old upscale')
    plt.plot(ua[:,2],fa[:,1],color='red',linewidth=3,label='new uspcale')
    plt.plot(xa,ya,color='black',linestyle='dashed',label='original',alpha=0.3)
    plt.fill_between(ua[:,2],ua[:,3],fa[:,1],color='pink',alpha=0.3)
    plt.xlim(ua[:,2].min(),ua[:,2].max())
    plt.ylim(min([ua[:,3].min(),fa[:,1].min()]),max([ua[:,3].max(),fa[:,1].max()]))
    plt.legend()
    plt.show()
    
    plt.hist([ua[:,3],fa[:,1],ya],bins=20,color=['black','red','gray'],normed=True)
    plt.show()
    #"""
    # #######################################################################    
    plt.plot(ub[:,2],ub[:,3],color='black',linewidth=3,label='old uspcale')
    plt.plot(ub[:,2],fb[:,1],color='red',linewidth=3,label='new uspcale')
    plt.plot(xb+1,yb,color='black',linestyle='dashed',label='original',alpha=0.3)
    plt.fill_between(ub[:,2],ub[:,3],fb[:,1],color='pink',alpha=0.3)
    plt.xlim(ub[:,2].min(),ub[:,2].max())
    plt.ylim(min([ub[:,3].min(),fb[:,1].min()]),max([ub[:,3].max(),fb[:,1].max()]))
    plt.legend()
    plt.show()
    
    plt.hist([ub[:,3],fb[:,1],yb],bins=20,color=['black','red','gray'],normed=True)
    plt.show()
        
w,ta,tb = initial_stats()
a,b,ha,hb = load_stuff()
iterate(a,b,ha,hb,w,ta,tb,10000)


    
