import numpy as np
import scipy.linalg as sl
import matplotlib.pyplot as plt
import pandas as pd

import glob, pickle, json, csv

import hasasia.sensitivity as hsen
import hasasia.sim as hsim
import hasasia.skymap as hsky

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['figure.figsize'] = [5,3]
mpl.rcParams['text.usetex'] = True

from enterprise.pulsar import Pulsar as ePulsar

#import par ant tim files and noise dictionary for DR2 PPTA
par_tim_dir = '/fred/oz002/vdimarco/sky_scrambles/sensitivity_curves/ppta_dr2_noise_analysis/data/'
noise_dir = '/fred/oz002/vdimarco/sky_scrambles/sensitivity_curves/ppta_dr2_noise_analysis/noisefiles/'
pars = sorted(glob.glob(par_tim_dir+'*.par'))
tims = sorted(glob.glob(par_tim_dir+'*.tim'))
noise_files = sorted(glob.glob(noise_dir+'*.json'))

#this file contains the median TOAs errors
#ensure this is in the correct order (i.e. sorted by pulsar name as below)
toas_med_err_file = '/fred/oz002/vdimarco/sky_scrambles/sensitivity_curves/PPTA_med_TOAerr3.csv'
toas_med_err = []
with open(toas_med_err_file, 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        toas_med_err.append(row)
toas_med_err = np.array(toas_med_err)
toas_med_err = toas_med_err.astype(float)

#extract name of 26 pulsars
psr_list = noise_files[:]
for i in range(len(noise_files)):
    psr_list[i] = psr_list[i].replace("_noise.json", "")
    psr_list[i] = psr_list[i].replace("/fred/oz002/vdimarco/sky_scrambles/sensitivity_curves/ppta_dr2_noise_analysis/noisefiles/", "")
#n is the number of pulsars, should be 26
n = len(psr_list)
print(n)


#finalise par and tim files and noise dictionary
noise = {}
for i in range(len(noise_files)):
    with open(noise_files[i]) as json_file:
        n = json.load(json_file)
    noise.update(n)


#load pulsar into enterprise.pulsar.Pulsar class instance
ePsrs = []
for par,tim in zip(pars,tims):
    ePsr = ePulsar(par, tim,  ephem='DE421')
    ePsrs.append(ePsr)
    print('\rPSR {0} complete'.format(ePsr.name),end='',flush=True)


#this is a simplified function that constructs the correlation matrix
def make_corr(toas_med):
    med_avg = (1/n)*sum(toas_med)
    scaled = toas_med/med_avg
    corr = np.zeros((n,n))
    sigma_sqr = np.zeros(n)
    sigma_sqr = np.power(scaled,2)
    corr = np.diag(sigma_sqr)
    return corr

#def make_corr(psr):
#    N = psr.toaerrs.size
#    corr = np.zeros((N,N))
#    _, _, fl, _, bi = hsen.quantize_fast(psr.toas,psr.toaerrs,
#                                         flags=psr.flags['f'],dt=1)
#    keys = [ky for ky in noise.keys() if psr.name in ky]
#    backends = np.unique(psr.flags['f'])
#    sigma_sqr = np.zeros(N)
#    ecorrs = np.zeros_like(fl,dtype=float)
#    for be in backends:
#        mask = np.where(psr.flags['f']==be)
#        key_ef = '{0}_{1}_{2}'.format(psr.name,be,'efac')
#        key_eq = '{0}_{1}_log10_{2}'.format(psr.name,be,'equad')
#        sigma_sqr[mask] = (noise[key_ef]**2 * (psr.toaerrs[mask]**2)
#                           + (10**noise[key_eq])**2)
#        mask_ec = np.where(fl==be)
#        key_ec = '{0}_{1}_log10_{2}'.format(psr.name,be,'ecorr')
#        ecorrs[mask_ec] = np.ones_like(mask_ec) * (10**noise[key_ec])
#    j = [ecorrs[ii]**2*np.ones((len(bucket),len(bucket)))
#         for ii, bucket in enumerate(bi)]
#
#    J = sl.block_diag(*j)
#    corr = np.diag(sigma_sqr) + J
#    return corr
#


#define red noise values for 12 years data
key_logA = 'red_noise_log10_A'
key_gamma = 'red_noise_gamma'
log10_A_1 = [value for key, value in noise.items() if key_logA in key]
log10_A_2 = [10**x for x in log10_A_1]
gamma = [value for key, value in noise.items() if key_gamma in key]
name = [key for key, value in noise.items() if key_logA in key]
for i in range(len(name)):
    name[i] = name[i].replace("_red_noise_log10_A", "")
values = [list(t) for t in zip(log10_A_2, gamma)]
rn_psrs = {}
for i in name:
    for j in values:
        rn_psrs[i] = j
        values.remove(j)
        break

#retrieve timespan across set of pulsars
Tspan = hsen.get_Tspan(ePsrs)

#set the frequency array
fyr = 1/(365.25*24*3600)
freqs = np.logspace(np.log10(1/(5*Tspan)),np.log10(2e-7),600)

#instantiate hasasia.Pulsar class instances using enterprise
psrs = []
thin = 1
for ePsr in ePsrs:
    corr = make_corr(ePsr)
    plaw = hsen.red_noise_powerlaw(A=9e-16, gamma=13/3., freqs=freqs)
    if ePsr.name in rn_psrs.keys():
        Amp, gam = rn_psrs[ePsr.name]
        plaw += hsen.red_noise_powerlaw(A=Amp, gamma=gam, freqs=freqs)
/Users/vdim0001/Desktop/PPTA_med_TOAerr.csv 
    corr += hsen.corr_from_psd(freqs=freqs, psd=plaw,
                               toas=ePsr.toas[::thin])
    psr = hsen.Pulsar(toas=ePsr.toas[::thin],
                      toaerrs=ePsr.toaerrs[::thin],
                      phi=ePsr.phi,theta=ePsr.theta,
                      N=corr, designmatrix=ePsr.Mmat[::thin,:])
    psr.name = ePsr.name
    psrs.append(psr)
    del ePsr
    print('\rPSR {0} complete'.format(psr.name),end='',flush=True)

##instantiate hasasia.Pulsar class instances using enterprise
#psrs = []
#thin = 20
#for ePsr in ePsrs:
#    corr = make_corr(ePsr)[::thin,::thin]
#    plaw = hsen.red_noise_powerlaw(A=9e-16, gamma=13/3., freqs=freqs)
#    if ePsr.name in rn_psrs.keys():
#        Amp, gam = rn_psrs[ePsr.name]
#        plaw += hsen.red_noise_powerlaw(A=Amp, gamma=gam, freqs=freqs)
#/Users/vdim0001/Desktop/PPTA_med_TOAerr.csv 
#    corr += hsen.corr_from_psd(freqs=freqs, psd=plaw,
#                               toas=ePsr.toas[::thin])
#    psr = hsen.Pulsar(toas=ePsr.toas[::thin],
#                      toaerrs=ePsr.toaerrs[::thin],
#                      phi=ePsr.phi,theta=ePsr.theta,
#                      N=corr, designmatrix=ePsr.Mmat[::thin,:])
#    psr.name = ePsr.name
#    psrs.append(psr)
#    del ePsr
#    print('\rPSR {0} complete'.format(psr.name),end='',flush=True)
#
#instantiate hasasia.Spectrum
specs = []
for p in psrs:
    sp = hsen.Spectrum(p, freqs=freqs)
    _ = sp.NcalInv
    specs.append(sp)
    print('\rPSR {0} complete'.format(p.name),end='',flush=True)

#print outputs
with open("PPTA_S_I_2.csv", "w") as S_I_file:
    writer = csv.writer(S_I_file)
    for i in range(26):
        writer.writerow(specs[i].S_I)
with open("PPTA_frequencies_2.csv", "w") as freqs_file:
    writer = csv.writer(freqs_file)
    writer.writerow(specs[0].freqs)
sc = hsen.DeterSensitivityCurve(specs)
with open("PPTA_Omega_2.csv", "w") as Omega_file:
    writer = csv.writer(Omega_file)
    writer.writerow(sc.Omega_gw)
