from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import csv

raj_file = '/fred/oz002/vdimarco/EPTA_scrambles/results/EPTA_psr_rajd.txt'
dec_file = '/fred/oz002/vdimarco/EPTA_scrambles/results/EPTA_psr_dec.txt'

with open(raj_file, "r") as f:
    raj = f.readlines()

rajs = [float(line.strip()) for line in raj]

with open(dec_file, "r") as f:
    dec = f.readlines()

decs = [float(line.strip()) for line in dec]

def get_angles(ra_in, dec_in):
    ra_deg = ra_in * u.deg
    dec_deg = dec_in * u.deg
    coord = SkyCoord(ra=ra_deg, dec=dec_deg)

    n = len(ra_in)
    angles = []
    for i in range(n):
        for j in range(i+1, n):
            angle = coord[i].separation(coord[j]).degree
            angles.append(angle)
            
    return angles


ang_deg = get_angles(rajs, decs)
ang_deg = np.array(ang_deg)


def get_HD_curve(zeta):

    #coszeta = np.cos(zeta)
    coszeta = np.cos(zeta*np.pi/180.)
    xip = (1.-coszeta) / 2.
    HD = 3.*( 1./3. + xip * ( np.log(xip) -1./6.) )
    
    return HD/2

HD = get_HD_curve(ang_deg)

print(len(HD))

with open('/fred/oz002/vdimarco/EPTA_scrambles/results/EPTA_ORF_true.csv', mode = "w") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(HD)

#plt.scatter(ang_deg, HD)
#plt.savefig('HD_plot_dr3.pdf', bbox_inches='tight' )
