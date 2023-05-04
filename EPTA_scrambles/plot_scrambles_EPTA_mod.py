import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib import rc
#import matplotlib.ticker as tick
#from matplotlib.ticker import FuncFormatter
import os


os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin/'
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.style.use('plot_style.txt')


#-- function formtting tick lables with 10^x at every tick
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

n_tries_PPTA = [int(item[0]) for item in n_t]

with open('/fred/oz002/vdimarco/sky_scrambles/skies/results/2_final_NANOGrav_n_tries_many_mod.csv', 'r') as file:
    csvreader = csv.reader(file)
    n_t = []
    for row in csvreader:
        n_t.append(row)
        
n_tries_NANOGrav = [int(item[0]) for item in n_t]

with open('/fred/oz002/vdimarco/IPTA_scrambles/results/IPTA_n_tries_many_mod.csv', 'r') as file:
    csvreader = csv.reader(file)
    n_t = []
    for row in csvreader:
        n_t.append(row)

n_tries_IPTA = [int(item[0]) for item in n_t] 

with open('/fred/oz002/vdimarco/EPTA_scrambles/results/EPTA_n_tries_many_mod.csv', 'r') as file:
    csvreader = csv.reader(file)
    n_t = []
    for row in csvreader:
        n_t.append(row)

n_tries_EPTA = [int(item[0]) for item in n_t]

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

cum_accepted_PPTA = np.cumsum(n_tries_PPTA)
cum_accepted_NANOGrav = np.cumsum(n_tries_NANOGrav)
cum_accepted_IPTA = np.cumsum(n_tries_IPTA)
cum_accepted_EPTA = np.cumsum(n_tries_EPTA) 
cum_num_of_tries_PPTA = np.arange(0, len(n_tries_PPTA))
cum_num_of_tries_NANOGrav = np.arange(0, len(n_tries_NANOGrav))
cum_num_of_tries_IPTA = np.arange(0, len(n_tries_IPTA))
cum_num_of_tries_EPTA = np.arange(0, len(n_tries_EPTA)) 

print('PPTA')
print(len(cum_num_of_tries_PPTA))
print(cum_accepted_PPTA)
print('NG')
print(len(cum_num_of_tries_NANOGrav))
print(cum_accepted_NANOGrav)
print('IPTA')
print(len(cum_num_of_tries_IPTA))
print(cum_accepted_IPTA)
print('EPTA')
print(len(cum_num_of_tries_EPTA))
print(cum_accepted_EPTA)

#-- old code
#num_of_tries = count_zeros(n_tries)
#efficiency = np.reciprocal([float(i) for i in num_of_tries])
#cum_eff = cumulative(efficiency)

#cum_num_of_tries_old = cumulative(num_of_tries)
#cum_accepted_old = np.arange(len(num_of_tries))
#cum_accepted_old = np.sum(n_tries)

#-- old code
#plt.step(cum_num_of_tries, cum_accepted)
#plt.minorticks_on()

#plt.xlabel('$N_{{\mathrm{{proposed}}}}$', fontsize=18)
#plt.ylabel('$N_{{\mathrm{{accepted}}}}$', fontsize=18)
#labs = plt.xticklabels()
#locs, labs = plt.xticks()

#plt.xticks(locs, [SciNotation(l, 2) for l in locs])
#plt.xlim(-100, None)
#plt.ylim(0,None)


#-- new code
fig, ax = plt.subplots()
plt.step(cum_num_of_tries_PPTA, cum_accepted_PPTA, color = 'darkviolet', label='PPTA')
plt.step(cum_num_of_tries_NANOGrav, cum_accepted_NANOGrav, linestyle='--', color = 'darkorange',label='NANOGrav')
plt.step(cum_num_of_tries_IPTA, cum_accepted_IPTA, linestyle='dotted', color = 'forestgreen', label='IPTA') 
plt.step(cum_num_of_tries_EPTA, cum_accepted_EPTA, linestyle='dashdot', color = 'blue', label='EPTA')
#xticks = np.linspace(0, 40000, 9)
#xticklabels = ['{:.1f}'.format(i/10000) for i in xticks]
#ax.set_xticks(xticks)
#ax.set_xticklabels(xticklabels) 


#locs, labs = plt.xticks()
#plt.xticks(locs, [SciNotation(l, 2) for l in locs]) 


plt.minorticks_on()
plt.xlabel('$N_{{\mathrm{{proposed}}}}$', fontsize=18)
plt.ylabel('$N_{{\mathrm{{accepted}}}}$', fontsize=18)
plt.legend(loc='upper right')
plt.xlim(-100, None)
plt.ylim(0, 35)

xlim = ax.get_xlim() 
#ax.annotate('x10$^4$', xy=(xlim[1], 0), xytext=(5, -20), textcoords='offset points', ha='right', va='top', fontsize=12)

plt.savefig('/fred/oz002/vdimarco/IPTA_scrambles/plots/IPTA_PPTA_NG_EPTA_scrambles_mod.pdf', bbox_inches='tight' )


plt.show()
