import numpy as np
import scipy.linalg as sl
import matplotlib.pyplot as plt
import pandas as pd

import glob, pickle, json, csv, os

import hasasia.sensitivity as hsen
import hasasia.sim as hsim
import hasasia.skymap as hsky

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['figure.figsize'] = [5,3]
mpl.rcParams['text.usetex'] = True

from enterprise.pulsar import Pulsar as ePulsar

#extract name of pulsars and order them
IPTA_dir = '/fred/oz002/vdimarco/IPTA_scrambles/DR2/release/VersionB'
psr_names = [name for name in os.listdir(IPTA_dir)]
#psr_names = ['J1824-2452A']
psr_list = sorted(psr_names)

#extract par and tim files
dirs = sorted(glob.glob('/fred/oz002/vdimarco/IPTA_scrambles/DR2/release/VersionB/*/'))
tims = []
pars = []
for directory in dirs:
    tim = glob.glob(directory+'*.tim')
    #tim = glob.glob(directory+'*J1824-2452A*.tim')
    par = glob.glob(directory+'*2.par')
    #par = glob.glob(directory+'*J1824-2452A*2.par')
    tims.extend(tim)
    pars.extend(par)

#create the noise dictionary
lines_all = []
for i in range(len(psr_list)):
    with open(pars[i], 'r') as file:
        lines = file.readlines()
    lines_red = lines[-4:-2]
    lines_all.extend(lines_red)

for i, string in enumerate(lines_all):
    lines_all[i] = " ".join(string.split()[1:])
    
noise_vals = [float(x) for x in lines_all]
print(noise_vals)

noise_strings = []
for psr in psr_list:
    noise_strings.append(psr + "_red_noise_log10_A")
    noise_strings.append(psr + "_red_noise_gamma")
noise = dict(zip(noise_strings, noise_vals))


#load pulsar into enterprise.pulsar.Pulsar class instance
ePsrs = []
for par,tim in zip(pars,tims):
    ePsr = ePulsar(par, tim)
    ePsrs.append(ePsr)
    print('\rPSR {0} complete'.format(ePsr.name),end='',flush=True)



#construct the correlation matrix
def make_corr(psr):
    N = psr.toaerrs.size
    corr = np.zeros((N,N))
    sigma_sqr = np.zeros(N)
    sigma_sqr = psr.toaerrs**2
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



#define red noise values
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
freqs = np.linspace(1/(Tspan),2e-7,600)

#instantiate hasasia.Pulsar class instances using enterprise
psrs = []
thin = 1
for ePsr in ePsrs:
    corr = make_corr(ePsr)[::thin,::thin]
    plaw = hsen.red_noise_powerlaw(A=9e-16, gamma=13/3., freqs=freqs)
    if ePsr.name in rn_psrs.keys():
        Amp, gam = rn_psrs[ePsr.name]
        plaw += hsen.red_noise_powerlaw(A=Amp, gamma=gam, freqs=freqs)

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

#instantiate hasasia.Spectrum
specs = []
for p in psrs:
    sp = hsen.Spectrum(p, freqs=freqs)
    _ = sp.NcalInv
    specs.append(sp)
    print('\rPSR {0} complete'.format(p.name),end='',flush=True)


with open("/fred/oz002/vdimarco/IPTA_scrambles/results/IPTA_S_I.csv", "w") as S_I_file:
    writer = csv.writer(S_I_file)
    for i in range(len(psr_list)):
        writer.writerow(specs[i].S_I)
with open("/fred/oz002/vdimarco/IPTA_scrambles/results/IPTA_frequencies.csv", "w") as freqs_file:
    writer = csv.writer(freqs_file)
    writer.writerow(specs[0].freqs)
sc = hsen.DeterSensitivityCurve(specs)
with open("/fred/oz002/vdimarco/IPTA_scrambles/results/IPTA_Omega.csv", "w") as Omega_file:
    writer = csv.writer(Omega_file)
    writer.writerow(sc.Omega_gw)

