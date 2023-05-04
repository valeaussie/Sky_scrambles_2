import numpy as np
import csv
import pandas as pd
import pdb
import time
import random
import matplotlib.pyplot as plt
import scipy

def compute_match(orf_true, phi1, phi2, Omega, f, Pij):

    n = len(Pij)
    ai = []
    ai_num = []
    for i in range(n):
        a = (f**(-7/3))/(Pij[i])
        a_num = (f**(-7/3))/(Pij[i])
        a_num = a_num.squeeze()
        a = a.squeeze()
        ai_num.append(a_num)
        ai.append(a)
    ai = np.array(ai)
    ai_num = np.array(ai_num)
    aij = []
    aij_num = []
    for i in range(n):
        for j in range(i+1, n):
            a = ai[i]*ai[j]
            #a_num = a*np.cos(2*phi2[j] - 2*np.array(phi1[i]))
            a_num = a*np.cos(phi1[j] - phi1[i] + phi2[j] - phi2[i])
            a = a.squeeze()
            a_num = a_num.squeeze()
            aij.append(a)
            aij_num.append(a_num)
    aij = np.array(aij)
    aij_num = np.array(aij_num)
    aij_s = []
    aij_num_s = []
    
    for i in range(len(orf_true)):
        a = np.sum(aij[i])
        a_num = np.sum(aij_num[i])
        a = a.squeeze()
        a_num = a_num.squeeze()
        aij_s.append(a)
        aij_num_s.append(a_num)
    aij_s = np.array(aij_s)
    aij_num_s = np.array(aij_num_s)


    num = sum(orf_true*orf_true*aij_num_s)
    den = sum(orf_true*orf_true*aij_s)

    match = np.abs(num/den)
    

    return match

def HD(x):                                                                      
    return 1.5*(1./3. + (1.-x)/2.*(np.log((1.-x)/2.)-1./6.))

def compute_orf(ptheta, pphi):                                                  
    npsr = len(ptheta)                                                          
    pos = [ np.array([np.cos(phi)*np.sin(theta),                                
                      np.sin(phi)*np.sin(theta),                                
                      np.cos(theta)]) for phi, theta in zip(pphi, ptheta) ]     
                                                                                
    x = []                                                                      
    for i in range(npsr):                                                       
        for j in range(i+1,npsr):                                               
            x.append(np.dot(pos[i], pos[j]))                                    
    x = np.array(x)                                                             
    orf = HD(x)                                                                 
                                                                                
    return orf


filename = '/fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles/IPTA_ORF_true.csv'

Omega_str = '/fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles/IPTA_Omega.csv' 
f_str = '/fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles/IPTA_frequencies.csv'
Pij_str = '/fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles/IPTA_S_I.csv'



Omega = []                                                                      
with open(Omega_str, 'r') as csv_file:                                          
    csv_reader = csv.reader(csv_file)                                           
    for row in csv_reader:                                                      
        Omega.append(row)                                                       
Omega = np.array(Omega)                                                         
Omega = Omega.astype(float)

f = []                                                                          
with open(f_str, 'r') as csv_file:                                              
    csv_reader = csv.reader(csv_file)                                           
    for row in csv_reader:                                                      
        f.append(row)                                                           
f = np.array(f)                                                                 
f = f.astype(float)

Pij = []                                                                        
with open(Pij_str, 'r') as csv_file:                                            
    csv_reader = csv.reader(csv_file)                                           
    for row in csv_reader:                                                      
        Pij.append(row)                                                         


Pij = np.array(Pij)                                                             
Pij = Pij.astype(float)                                                         


orf_true = np.loadtxt(filename, delimiter = ',') 

def compute_npsrs(npairs):                                                      
    n = (1 + np.sqrt(1+8*npairs))/2                                             
    return int(n)

def get_scrambles(orf_true, Omega, f, Pij,  N=10, thresh=0.1, resume =False):    
    #npsr = compute_npsrs(len(orf_true))
    npsr = 65
    if resume == False:
        print("not resume")
        orfs = []
        orfs.append(orf_true)
        phis = []
        phis.append(np.zeros([npsr,f.size]))
        n_accep = 1
        n_iter = 1
        rejections = []
    elif resume == True:
        print("resume")
        phis = list(np.load("/fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles/results/phase_IPTA_phis.npy"))
        rejections = list(pd.read_csv("/fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles/results/phase_IPTA.csv",index_col=None).values.squeeze())
        n_iter = len(rejections)
        n_accep = np.cumsum(rejections)[-1]

    st = time.time()
    while n_accep < N:
        matchs = []
        n_iter += 1

        phi2 = np.random.uniform(0, 2*np.pi,[npsr,f.size])

        for ii in range(n_accep):
            match = compute_match(orf_true, phis[ii], phi2, Omega, f, Pij)
            matchs = np.append(matchs, abs(match))

        matchs = np.array(matchs)
        reject = any(matchs > thresh)
        ct_rej = 1
        if reject == True:
            ct_rej = 0
            et = time.time()
            if et - st > 7200:
                return (matchs, n_accep, n_iter, rejections) 
        rejections.append(ct_rej)
        if reject == False:
            phis.append(phi2)
            n_accep += 1
            np.save("/fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles/results/phase_IPTA_phis",np.array(phis))
            
            print("accepted phase scrambles = {}".format(n_accep))  
            st = time.time()
            pd.DataFrame(rejections).to_csv('/fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles/results/phase_IPTA.csv', header=None, index=None) 
            plt.plot(np.cumsum(np.array(rejections)))
            plt.title("Phase scrambling; threshold: 0.1")
            plt.xlabel("Number of iterations")
            plt.ylabel("Cumulative number of accepted scrambles")
            plt.savefig("/fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles/plots/phase_IPTA.png")
            plt.clf()




    return (matchs, n_accep, n_iter, rejections)

#Tspan = 474459299.74052906
#newf = np.arange(2.1076623443715344e-09,2.00000000e-07,1/Tspan)
newf = np.arange(f[0].min(),2.00000000e-07,f[0].min())

newPij = []

for P in Pij:
    func = scipy.interpolate.interp1d(f[0],P)
    newP = func(newf)
    newPij.append(newP)

newPij = np.array(newPij)

iters = []
#while i < 10:
m, n_a, n_i, c = get_scrambles(orf_true, Omega, newf, newPij, 10000, 0.1, resume=True)
#iters.append(n_i)
#i = i+1


print(n_i)
pd.DataFrame(c).to_csv('/fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles/results/phase_IPTA.csv', header=None, index=None) 
plt.plot(np.cumsum(np.array(c)))
plt.title("Phase scrambling; threshold: 0.1")
plt.xlabel("Number of iterations")
plt.ylabel("Cumulative number of accepted scrambles")
plt.savefig("/fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles/plots/phase_IPTA.png")
plt.clf()
