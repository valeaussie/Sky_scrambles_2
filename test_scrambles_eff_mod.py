import numpy as np
import makeskyscrambles_eff_mod as sky
import open_par
import timeit
import pandas as pd
import csv
import matplotlib.pyplot as plt 

t_0 = timeit.default_timer()

filename = '/fred/oz002/vdimarco/sims_feb_23/shuffles/json_files_uplim/json_61/outputs/crosscorr_random_missp_GWmin15_5.csv'
Omega_str = '/fred/oz002/vdimarco/sky_scrambles/sensitivity_curves/PPTA_Omega_thin.csv'
f_str = '/fred/oz002/vdimarco/sky_scrambles/sensitivity_curves/PPTA_frequencies_thin.csv'
Pij_str = '/fred/oz002/vdimarco/sky_scrambles/sensitivity_curves/PPTA_S_I_thin.csv'

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

scrambles = sky.get_scrambles(orf_true[1], Omega, f, Pij, 10, 1000000000000, 0.1, save=False, filename='/fred/oz002/vdimarco/sky_scrambles/final_scrambles_eff_mod_10', resume=False)

#def count_zeros(arr):
#    result = []
#    count = 0
#    for i in arr:
#        if i == 0:
#            count += 1
#        else:
#            result.append(count+1)
#            count = 0
#    return result

#efficiency = count_zeros(scrambles[4])
#print(efficiency)

#cum_accepted = np.arange(len(efficiency))
#print(cum_accepted)
#plt.scatter(efficiency, cum_accepted)
#plt.yscale('log')
#plt.savefig('test_efficiency_mod_1e4.png')

#pd.DataFrame(scrambles[4]).to_csv('/fred/oz002/vdimarco/sky_scrambles/test_scrambles_eff_3e5.csv', header=None, index=None)

#t_1 = timeit.default_timer()
#elapsed_time = round((t_1 - t_0) * 10 ** 6, 3)
#print(f"Elapsed time: {elapsed_time} µs")

