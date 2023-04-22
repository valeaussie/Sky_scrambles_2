import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib import rc
import os
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin/'
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.style.use('plot_style.txt')


def SciNotation(num,sig):
    x = '%.1e'  %num  #<-- Instead of 2, input sig here
    #formatted_x = '%.1e' % x
    x = x.split('e')
    if num == 0:
        return r"$0$"
    if (x[1])[0] == "-":
        return r"${}\times 10^{{{}}}$".format(x[0],x[1].lstrip('0'))
    else:
        return r"${}\times 10^{{{}}}$".format(x[0],x[1][1:].lstrip('0'))

with open('/fred/oz002/vdimarco/sky_scrambles/skies/results/final_PPTA_n_tries_many_mod_2h.csv', 'r') as file:
    csvreader = csv.reader(file)
    n_t = []    
    for row in csvreader:
        n_t.append(row)

n_tries = [int(item[0]) for item in n_t]


def count_zeros(arr):
    result = []
    count = 0
    for i in arr:
        if i == 0:
            count += 1
        else:
            result.append(count+1)
            count = 0
    return result

def cumulative(lists): 
    cu_list = [] 
    length = len(lists) 
    cu_list = [sum(lists[0:x:1]) for x in range(0, length+1)] 
    return cu_list[1:]


num_of_tries = count_zeros(n_tries)
#efficiency = np.reciprocal([float(i) for i in num_of_tries])
#cum_eff = cumulative(efficiency)
cum_num_of_tries = cumulative(num_of_tries)
cum_accepted = np.arange(len(num_of_tries))


plt.step(cum_num_of_tries, cum_accepted)
plt.minorticks_on()

#plt.grid(linestyle = ':', color = 'k', alpha = 0.3)
#plt.xscale('log') 
#plt.yscale('log')
plt.xlabel('$N_{{\mathrm{{trial}}}}$', fontsize=14)
plt.ylabel('$N_{{\mathrm{{accept}}}}$', fontsize=14)
#labs = plt.xticklabels()
locs, labs = plt.xticks()

plt.xticks(locs, [SciNotation(l, 2) for l in locs])
plt.xlim(-100, None)
plt.ylim(0,None)

plt.savefig('final_PPTA_mod.pdf', bbox_inches='tight' )


