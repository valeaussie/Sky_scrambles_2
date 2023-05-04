import numpy as np                                                             
import pandas as pd                                                            
                                                                               
import glob, pickle, json, csv, os                                             
                                                                               
#extract name of pulsars and order them                                        
IPTA_dir = '/fred/oz002/vdimarco/IPTA_scrambles/DR2/release/VersionB'          
psr_names = [name for name in os.listdir(IPTA_dir)]                            
psr_list = sorted(psr_names) 


np.savetxt("/fred/oz002/vdimarco/IPTA_scrambles/results/IPTA_psr_list.txt", psr_list, fmt='%s')
