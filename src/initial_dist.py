'''
Created on Mar 27, 2016
generate ring shape intial distribution beam. uniform_ringwmom is the beam with the initial momenetum spread, it need initial ek with certain emittance. 
and uniform_ring function shows the zero emittance distribution
@author: wange
'''
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import norm
import scipy.constants.constants as conts
#import pprint
import time
#unit:cm
start_time=time.time()
c=conts.c
mu0=conts.mu_0
eps0=conts.epsilon_0
e=conts.e
e0=0.511*10**6

def uniform_ringwmom(tranE,in_rad,out_rad,partic_numbs,bunch_charge,ek,bunch_length):
    z_mom=(ek+e0/e0)*(1-(e0/(e0+ek))**2)**0.5#beta*gamma longitudinal
    bl=(1-(ek+e0/e0)**-2)**0.5*c*bunch_length*10**-9# bunch length :cm
    t_mom=(tranE+e0/e0)*(1-(e0/(e0+tranE))**2)**0.5#beta*gamma transverse
    radius=out_rad*0.1*np.sqrt(np.random.uniform(in_rad/out_rad,out_rad/out_rad,partic_numbs))#cm
    theta_pos=np.random.uniform(-math.pi,math.pi,partic_numbs)
    theta_mom=np.random.uniform(-math.pi,math.pi,partic_numbs)

    trans_mom=np.random.uniform(0,t_mom,partic_numbs)
    charge=bunch_charge*10**(-9)/partic_numbs#C/per-maro-p
    
    xcoordi=radius*np.cos([theta_pos])#cm
    x_mom=trans_mom*np.cos([theta_mom])#rad
    ycoordi=radius*np.sin([theta_pos])#cm
    y_mom=trans_mom*np.sin([theta_mom])#rad
    zcoordi=np.array([np.random.uniform(-bl,0,partic_numbs)])#cm
    z_momi=np.empty(partic_numbs); z_momi.fill(z_mom)
    z_momi=np.array([z_momi])#cm/s
    mcharge=np.empty(partic_numbs); mcharge.fill(charge)
    mcharge=np.array([mcharge])#C/pmaro
    
    initial_beam=np.concatenate((xcoordi.T,x_mom.T,ycoordi.T,y_mom.T,zcoordi.T,z_momi.T,mcharge.T),axis=1)
    #initial_beam=np.concatenate((xcoordi,x_mom,ycoordi,y_mom,zcoordi,z_voli,mcharge),axis=0)
    #initial_beam=initial_beam.T
    #print(trans_mom)
    #plt.hist(trans_mom,bins=200,normed=False)
    #plt.show()
    #print(y_mom)

    return initial_beam
    #print(initial_beam, initial_beam[:,5])

def uniform_ring(freq,in_rad,out_rad,partic_numbs,bunch_length):#uniform distribution without mom
    del_phs= freq*10**6*360*bunch_length*10**-9
    del_phsi=np.random.uniform(-del_phs/2.0,del_phs/2.0,partic_numbs)
    del_phsi=np.array([del_phsi])
    radius=out_rad*0.1*np.sqrt(np.random.uniform(in_rad/out_rad,out_rad/out_rad,partic_numbs))#cm
    theta_pos=np.random.uniform(-math.pi,math.pi,partic_numbs)
    xcoordi=radius*np.cos([theta_pos])#cm
    ycoordi=radius*np.sin([theta_pos])#cm
    x_momi=np.zeros(partic_numbs)
    x_momi=np.array([x_momi])
    y_momi=np.zeros(partic_numbs)
    y_momi=np.array([y_momi])
    del_w=np.zeros(partic_numbs)
    del_w=np.array([del_w])
    initial_beam=np.concatenate((xcoordi.T,x_momi.T,ycoordi.T,y_momi.T,del_phsi.T,del_w.T),axis=1)
    #print(initial_beam)
    return initial_beam
    
def plot_dis(des_dist,out_rad):
    out_rad=out_rad*0.1
    plt.figure(figsize=(10,10))
    plt.subplot(221)
    plt.xlim(-2*out_rad,2*out_rad)
    plt.ylim(-2*out_rad,2*out_rad)
    plt.plot(des_dist[:,0],des_dist[:,2],'b,') 
    plt.subplot(223)
    if np.mean(des_dist[:,1])!=0:
        plt.hist(np.arctan(des_dist[:,3]/des_dist[:,1]), bins=30, normed=False) 
    else:
        plt.hist(des_dist[:,1], bins=30, normed=False) 
    plt.subplot(222)
    plt.hist(des_dist[:,4],bins=20,normed=False)
    plt.subplot(224)
    #plt.xlim(0,10000)
    plt.hist((des_dist[:,3]**2+des_dist[:,1]**2)**0.5,bins=30,normed=False)
    plt.show()

def charact(des_dist,out_rad,in_rad,partic_numbs):  
    rmsx=np.std(des_dist[:,0])*10
    rmsy=np.std(des_dist[:,2])*10
    temit_x=(np.std(des_dist[:,0])**2*np.std(des_dist[:,1])**2-np.mean(des_dist[:,0]*des_dist[:,1])**2)**0.5*10**4
    temit_y=(np.std(des_dist[:,2])**2*np.std(des_dist[:,3])**2-np.mean(des_dist[:,2]*des_dist[:,3])**2)**0.5*10**4
    mean_z=np.mean(des_dist[4])
    thermemit=(temit_x*temit_y)**0.5/(out_rad/2)
    with open('distri_od_%.1fid_%.1fpn_%dte_%.2f.txt' %(out_rad,in_rad,partic_numbs,thermemit),'wb') as dat:
        np.savetxt(dat,des_dist,fmt='%.8e',delimiter=' ')
    print('Zmean=%.4f[mm/deg],Xrms= %.4f[mm], Yrms= %.4f[mm],Xemit=%.4f[mm-mrad],Yemit=%.4f[mm-mrad],thermal_emittance=%.4f[mrad]' %(mean_z,rmsx,rmsy,temit_x,temit_y,thermemit) )

if __name__ == '__main__':
    bunch_charge=-5.3#nC,-:e
    tranE=0.1 # eV
    out_rad=6 #mm
    in_rad=0 #mm
    partic_numbs=9999 #1
    bunch_length=1.0#ns
    ek=0.4#eV
    freq=112#MHz
    if tranE>=0.00001:
        des_dist=uniform_ringwmom(tranE,in_rad,out_rad,partic_numbs,bunch_charge,ek,bunch_length)
    elif tranE<0.00001:
        des_dist=uniform_ring(freq,in_rad,out_rad,partic_numbs,bunch_length)
    print('current density=%.2f[A/cm^2]' %(-bunch_charge*100/bunch_length/(math.pi*(out_rad**2-in_rad**2))))
    charact(des_dist,out_rad,in_rad,partic_numbs)
    plot_dis(des_dist,out_rad)
    #print(des_dist)
    
    