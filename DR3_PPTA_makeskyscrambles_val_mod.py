import numpy as np
import csv
import pandas as pd
import time
import pdb

def compute_match(orf1, orf2, Omega, f, Pij):
    # hack
    #Omega = 0*Omega + 1
    #f = 0*f + 1
    #Pij = 0*Pij + 1

    n = len(Pij)
    ai = []
    for i in range(n):
        a = f**(-7/3)/Pij[i]
        a = a.squeeze()
        ai.append(a)
    ai = np.array(ai) 
    aij = []
    for i in range(n):
        for j in range(i+1, n):
            a = ai[i]*ai[j]
            a = a.squeeze()
            aij.append(a)
    aij = np.array(aij)
    aij_s = []
    for i in range(len(orf1)):
        a = np.sum(aij[i])
        a = a.squeeze()
        aij_s.append(a)
    aij_s = np.array(aij_s)
    num = sum(orf1*orf2*aij_s)
    den1 = np.sqrt(sum(orf1*orf1*aij_s))
    den2 = np.sqrt(sum(orf2*orf2*aij_s))
    den = den1*den2
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

filename = '/fred/oz002/vdimarco/sky_scrambles/skies/sensitivity_curves/ppta_dr3/PPTA_ORF_true_DR3.csv'
Omega_str = '/fred/oz002/vdimarco/sky_scrambles/skies/sensitivity_curves/results/PPTA_DR3_Omega_alltimes.csv' 
f_str = '/fred/oz002/vdimarco/sky_scrambles/skies/sensitivity_curves/results/PPTA_DR3_frequencies_alltimes.csv'
Pij_str = '/fred/oz002/vdimarco/sky_scrambles/skies/sensitivity_curves/results/PPTA_DR3_S_I_alltimes.csv'

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

orf_true= []
with open(filename, 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        orf_true.append(row)

#orf_true = np.loadtxt(filename, delimiter = ',') 

orf_true = np.array(orf_true)
orf_true = orf_true.astype(float)
orf_true = orf_true.squeeze()

def compute_npsrs(npairs):                                                      
    n = (1 + np.sqrt(1+8*npairs))/2                                             
    return int(n)

def get_scrambles(orf_true, Omega, f, Pij,  N=10, thresh=0.1):
    npsr = compute_npsrs(len(orf_true))
    orfs = []
    thetas = []
    phis = []
    orfs.append(orf_true)
    n_accep = 1
    n_iter = 1
    rejections = []
    start = time.time()
    while n_accep < N:
        matchs = []
        n_iter += 1
        #print("this is iteration = {}".format(n_iter))
        #print("accepted skies = {}".format(n_accep))   
        ptheta = np.arccos(np.random.uniform(-1, 1, npsr))
        pphi = np.random.uniform(0, 2*np.pi, npsr)
        new_orf = compute_orf(ptheta, pphi)
        for ii in range(n_accep):
            match = compute_match(new_orf, orfs[ii], Omega, f, Pij)
            matchs = np.append(matchs, abs(match))
            #print(match)
        matchs = np.array(matchs)
        reject = any(matchs > thresh)
        #pdb.set_trace()
        ct_rej = 1
        if reject == True:
            ct_rej = 0
            end = time.time()
            #print("time {}".format((end-start)))
            if end - start > 7200:
                return(n_accep, n_iter, rejections, orfs, thetas, phis)
        rejections.append(ct_rej)
        if reject == False:
            orfs.append(new_orf)
            thetas.append(ptheta)
            phis.append(pphi)
            n_accep += 1
            start = time.time()
            #print("accepted")
        #else:
            #print("rejected")
    return(n_accep, n_iter, rejections, orfs, thetas, phis)


n_a, n_i, c, orfs, thetas, phis = get_scrambles(orf_true, Omega, f, Pij, 1000, 0.1)

#print(c)

pd.DataFrame(c).to_csv('/fred/oz002/vdimarco/sky_scrambles/skies/results/DR3_PPTA_n_tries_many_mod_2h_2.csv', header=None, index=None) 
pd.DataFrame(orfs).to_csv('/fred/oz002/vdimarco/sky_scrambles/skies/results/DR3_mod_orfs.csv', header=None, index=None)
pd.DataFrame(thetas).to_csv('/fred/oz002/vdimarco/sky_scrambles/skies/results/DR3_mod_thetas.csv', header=None, index=None)
pd.DataFrame(phis).to_csv('/fred/oz002/vdimarco/sky_scrambles/skies/results/DR3_mod_phis.csv', header=None, index=None)  
