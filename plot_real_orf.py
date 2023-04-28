import os
import json
import matplotlib.pyplot as plt
import numpy as np
import glob
import csv
import pandas as pd

from enterprise.signals import signal_base, selections
from enterprise.signals import gp_signals
from enterprise.pulsar import Pulsar
from enterprise_extensions import model_utils, blocks
from enterprise_extensions.frequentist import optimal_statistic as opt_stat

plt.style.use('plot_style.txt')

vec_angsep = []

vec_crosscorr1 = []
vec_crosscorr_error1 = []
vec_OSHD1 = []
vec_OSHD_error1 = []
vec_SNR_OSHD1 = []

    
# Load up the noise dictionary to get values for the white noise parameters
with open('/fred/oz002/vdimarco/sims_june_22/26Pulsars_chann_norn.json', 'r') as f:                          
    noisedict = json.load(f) 

# Load up the properly specified maximum-likelihood values for the pulsars' timing noise parameter                                                               
with open('/fred/oz002/vdimarco/sims_feb_23/26Pulsars_maxlike_new_missp_and_GW_min15_4.json', 'r') as f:                                                                        
    ml_params1 = json.load(f)


# set the data directory
datadir = '/fred/oz002/vdimarco/sims_june_22/26Pulsars_timn_only_100/output/real_61'

#set up pulsar objects 
parfiles = sorted(glob.glob(datadir + '/*par'))
timfiles = sorted(glob.glob(datadir + '/*tim'))

psrs = []
ephemeris = 'DE421'
for p, t in zip(parfiles, timfiles):
    psr = Pulsar(p, t, ephem=ephemeris)
    psrs.append(psr)

#set up model
Tspan = model_utils.get_tspan(psrs)
s = gp_signals.TimingModel()
by_band = selections.by_band
no_select = selections.no_selection
s += blocks.white_noise_block(vary=False, inc_ecorr=False, tnequad=True, select=no_select)
s += blocks.red_noise_block(prior='log-uniform', Tspan=Tspan, components=30)
s += blocks.common_red_noise_block(psd='powerlaw', prior='log-uniform', Tspan=Tspan, components=5, gamma_val=4.33333, name='gw')

# Set up the PTA object using the signal we defined above and the pulsars
pta = signal_base.PTA([s(p) for p in psrs])

# Set the white noise parameters to the values in the noise dictionary
pta.set_default_params(noisedict)

ostat = opt_stat.OptimalStatistic(psrs, pta=pta, orf='hd')

# Compute the optimal statistic
# The optimal statistic returns five quantities:
#  - xi: an array of the angular separations between the pulsar pairs (in radians)
#  - rho: an array of the cross-correlations between the pulsar pairs
#  - sig: an array of the uncertainty in the cross-correlations
#  - OS: the value of the optimal statistic
#  - OS_sig: the uncertainty in the optimal statistic
 
xi, rho, sig, OS, OS_sig = ostat.compute_os(params=ml_params1)

SNR_HD = OS/OS_sig

def get_HD_curve(zeta):
    
    coszeta = np.cos(zeta*np.pi/180.)
    xip = (1.-coszeta) / 2.
    HD = 3.*( 1./3. + xip * ( np.log(xip) -1./6.) )
    
    return HD/2


def weightedavg(rho, sig):
    weights, avg = 0., 0.
    for r,s in zip(rho,sig):
        weights += 1./(s*s)
        avg += r/(s*s)
        
    return avg/weights, np.sqrt(1./weights)

# sort the cross-correlations by the angle of separation

idx = np.argsort(xi)
angsep_sorted = xi[idx]
cross_sorted = rho[idx]
cross_error_sorted = sig[idx]

#-- bin the cross-correlations so that there are the same number of pairs per bin
npairs = 24

angsep_mean = []
i = 0
while i < len(angsep_sorted):
    angsep_mean.append(np.mean(angsep_sorted[i:npairs+i]))
    i += npairs

cross_avg = []
cross_error_avg = []
i = 0
while i < len(angsep_sorted):
    r, s = weightedavg(cross_sorted[i:npairs+i], cross_error_sorted[i:npairs+i])
    i += npairs
    print(r)    
    cross_avg.append(r)
    cross_error_avg.append(s)
    
angsep_mean_array = np.array(angsep_mean)


print("the signal to noise ratio is")
print(SNR_HD)


# plot the binned cross correlations, the HD curve, the monopole and the dipole
zeta = np.linspace(0.01,180,100)
HD = get_HD_curve(zeta+1)

(_, caps, _) = plt.errorbar(angsep_mean_array*180/np.pi, cross_avg, yerr=cross_error_avg, ls='', color='0.1', fmt='o', capsize=4, elinewidth=1.2)
plt.plot(zeta, OS*HD, ls='--', label='Hellings-Downs', color='C0', lw=1.5)

plt.xlim(0, 180)
#plt.title('Correlations for realisation 61 with SNR {}'.format(SNR_HD))
plt.ylabel(r'$\hat{A}^2 \Gamma_{ab}(\zeta)$')
plt.xlabel(r'$\zeta$ (deg)')

plt.savefig(('orf_true.pdf'),dpi=300, bbox_inches='tight')  


pd.DataFrame(rho).to_csv('/fred/oz002/vdimarco/sky_scrambles/cross_real_61_missp_4.csv', header=None, index=None)
