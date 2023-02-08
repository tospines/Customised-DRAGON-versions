import numpy as np
import os
import sys
from scipy.interpolate import RegularGridInterpolator



ChannelID = 'bb'
PartID = 'antideuterons'
Interp_Mass = 100 #GeV
Mode = 'Tune'  ## 'PPPC' 'Pythia' 'Tune'

## --------
DM_path = '/home/tospines/DRAGON/dragon-3.1.0/DM/'

if PartID ==   'antideuterons':  ID = 'Ad'
elif PartID == 'antiHelium3':  ID = 'AHe'
elif PartID == 'antiprotons':  ID = 'Ap'
elif PartID == 'positrons':  ID = 'Pos'
elif PartID == 'gammas':  ID = 'Gamma'
else:    print('Choose a correct Parent particle!!')
## --------

##### Define it depending on your interpolated file!!! ####

E_Flux, Masses = [], []
for myfile in os.listdir(DM_path):
    if Mode == 'PPPC':
        FileCondition =  ('-PPPC' in myfile) and ('{}_'.format(ChannelID) in myfile) and ('{}'.format(ID) in myfile)
        Strip = 'GeV_{}-PPPC.dat'.format(ID)
    elif Mode == 'Tune':
        FileCondition =  ('_tune' in myfile) and ('{}_'.format(ChannelID) in myfile) and ('{}'.format(ID) in myfile)
        Strip = 'GeV_{}_tune.dat'.format(ID)
    else:
        FileCondition =  ('-PPPC' not in myfile) and ('_tune' not in myfile) and ('{}_'.format(ChannelID) in myfile) and ('{}'.format(ID) in myfile)
        Strip = 'GeV_{}-PPPC.dat'.format(ID)
    if (FileCondition):
        E, Flux = np.loadtxt(DM_path + myfile, skiprows=0, usecols=(0, 1), unpack=True)
        E_Flux.append(np.array([E, Flux]))
        Masses.append(float(myfile.strip('{}_'.format(ChannelID)).strip(Strip) ))

        
### Interpolation...
Mass_vec, E_Flux = zip(*sorted(zip(Masses, E_Flux)))
Energy = np.array(E_Flux[-1][0])
inter = []
for im, mM in enumerate(Mass_vec):
    Gen_Flux = np.interp(Energy, E_Flux[im][0], E_Flux[im][1])
    Gen_Flux[Energy> Energy[np.argmin(np.abs(Energy-(mM/2.))) ]] = 0. 
    inter.append(Gen_Flux)
inter = np.array(inter).reshape(len(Mass_vec), len(Energy))
my_interp = RegularGridInterpolator((Mass_vec, Energy), inter)



### Writing table...
En_vec = np.logspace(np.log10(max(Interp_Mass/2.*1e-9, Energy[0]+1e-12)), np.log10(Interp_Mass/2.),  20*9) ## 20 points per decade!
f = open(ChannelID + "_" + str(Interp_Mass) + "GeV_" + ID + "-" + Mode + "_Interp.dat", 'w')
[f.write( "{0:.10e}".format(En) + "\t" + "{0:.10e}".format(my_interp((Interp_Mass, En))) + "\n") for iE, En in enumerate(En_vec)]
f.close()
