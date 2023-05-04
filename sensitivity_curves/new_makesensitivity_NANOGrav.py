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

#import par ant tim files and noise dictionary for the 12y NANOGrav
pardir = '/fred/oz002/vdimarco/sky_scrambles/skies/sensitivity_curves/12p5yr_stochastic_analysis/data/par/'
timdir = '/fred/oz002/vdimarco/sky_scrambles/skies/sensitivity_curves/12p5yr_stochastic_analysis/data/tim/'
noise_dir = '/fred/oz002/vdimarco/sky_scrambles/skies/sensitivity_curves/12p5yr_stochastic_analysis/data/'
pars = sorted(glob.glob(pardir+'*.par'))
tims = sorted(glob.glob(timdir+'*.tim'))
noise_files = sorted(glob.glob(noise_dir+'*.json'))

#import name of pulsars
data = pd.read_csv('/fred/oz002/vdimarco/sky_scrambles/skies/sensitivity_curves/NANOGrav_psr_list')
pulsar_names = data.iloc[:, 0].values
psr_list = pulsar_names.tolist()
psr_list.insert(0, 'B1855+09')

print(psr_list)

#finalise par and tim files and noise dictionary
def get_psrname(file,name_sep='_'):
    return file.split('/')[-1].split(name_sep)[0]

pars = [f for f in pars if get_psrname(f) in psr_list]
tims = [f for f in tims if get_psrname(f) in psr_list]
noise = {}
with open(noise_files[2]) as json_file:
    noise = json.load(json_file)

#load pulsar into enterprise.pulsar.Pulsar class instance
ePsrs = []
for par,tim in zip(pars,tims):
    ePsr = ePulsar(par, tim,  ephem='DE436')
    ePsrs.append(ePsr)
    print('\rPSR {0} complete'.format(ePsr.name),end='',flush=True)

#construct the correlation matrix
def make_corr(psr):
    N = psr.toaerrs.size
    corr = np.zeros((N,N))
    _, _, fl, _, bi = hsen.quantize_fast(psr.toas,psr.toaerrs,
                                         flags=psr.flags['f'],dt=1)
    keys = [ky for ky in noise.keys() if psr.name in ky]
    backends = np.unique(psr.flags['f'])
    sigma_sqr = np.zeros(N)
    ecorrs = np.zeros_like(fl,dtype=float)
    for be in backends:
        mask = np.where(psr.flags['f']==be)
        key_ef = '{0}_{1}_{2}'.format(psr.name,be,'efac')
        key_eq = '{0}_{1}_log10_{2}'.format(psr.name,be,'equad')
        sigma_sqr[mask] = (noise[key_ef]**2 * (psr.toaerrs[mask]**2)
                           + (10**noise[key_eq])**2)
        mask_ec = np.where(fl==be)
        key_ec = '{0}_{1}_log10_{2}'.format(psr.name,be,'ecorr')
        ecorrs[mask_ec] = np.ones_like(mask_ec) * (10**noise[key_ec])
    j = [ecorrs[ii]**2*np.ones((len(bucket),len(bucket)))
         for ii, bucket in enumerate(bi)]

    J = sl.block_diag(*j)
    corr = np.diag(sigma_sqr) + J
    return corr

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

#define red noise values for significant pulsars for 11 years data
#rn_psrs = {'B1855+09':[10**-13.7707, 3.6081],
#           'B1937+21':[10**-13.2393, 2.46521],
#           'J0030+0451':[10**-14.0649, 4.15366],
#           'J0613-0200':[10**-13.1403, 1.24571],
#           'J1012+5307':[10**-12.6833, 0.975424],
#           'J1643-1224':[10**-12.245, 1.32361],
#           'J1713+0747':[10**-14.3746, 3.06793],
#           'J1747-4036':[10**-12.2165, 1.40842],
#           'J1903+0327':[10**-12.2461, 2.16108],
#           'J1909-3744':[10**-13.9429, 2.38219],
#           'J2145-0750':[10**-12.6893, 1.32307],
#           }

#retrieve timespan across set of pulsars
Tspan = hsen.get_Tspan(ePsrs)

#set the frequency array
fyr = 1/(365.25*24*3600)
freqs = np.linspace(1/(Tspan),2e-7,600)

#instantiate hasasia.Pulsar class instances using enterprise
#note thinning by a facrtor of 10
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

print(specs[0].S_I)
print('test')
print(specs[1].S_I)

with open("NANOGrav_S_I_45_psrs.csv", "w") as S_I_file:
    writer = csv.writer(S_I_file)
    for i in range(len(psr_list)):
        writer.writerow(specs[i].S_I)
with open("NANOGrav_frequencies_45_psrs.csv", "w") as freqs_file:
    writer = csv.writer(freqs_file)
    writer.writerow(specs[0].freqs)
sc = hsen.DeterSensitivityCurve(specs)
with open("NANOGrav_Omega_45_psrs.csv", "w") as Omega_file:
    writer = csv.writer(Omega_file)
    writer.writerow(sc.Omega_gw)
with open("S_I.csv", "w") as S_I_file:
    writer = csv.writer(S_I_file)
    writer.writerow(sc.S_I)

#plot a samlple of sensitivity curves
#fig=plt.figure(figsize=[15,45])
#j = 1
#names = ['B1937+21','J0340+4130','J1024-0719',
#         'J1713+0747','J1853+1303','J1909-3744',]
#for sp,p in zip(specs,psrs):
#    if p.name in names:
#        fig.add_subplot(12,3,j)
#        a = sp.h_c[0]/2*1e-14
#        if p.name == 'J1024-0719':
#            alp = -5/2
#            a *= 8e-10
#            plt.loglog(sp.freqs[:150],a*(sp.freqs[:150])**(alp),
#                       color='C2',label=r'$f^{-5/2}$')
#        else:
#            alp = -3/2
#            plt.loglog(sp.freqs[:150],a*(sp.freqs[:150])**(alp),
#                       color='C1',label=r'$f^{-3/2}$')
#        plt.ylim(2e-15,2e-10)
#        plt.loglog(sp.freqs,sp.h_c, color='C0')
#        plt.rc('text', usetex=True)
#        plt.xlabel('Frequency [Hz]')
#        plt.ylabel('Characteristic Strain, $h_c$')
#        plt.legend(loc='upper left')
#        plt.title(p.name)
#        j+=1
#fig.tight_layout()
#plt.savefig('sensitivity_curves_NANOGrav.png', bbox_inches='tight')
#plt.close()
