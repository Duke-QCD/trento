#!/usr/bin/env python3

import numpy as np
import sys
import os
import subprocess as sp

config="""# TRENTo config
reduced-thickness = 0.0
constit-number = {nparton}
nucleon-width = {wnucleon}
constit-width = {wparton}
nucleon-min-dist = 0
fluctuation = {variance}

### Fragmentation region
kT-min = {kTmin} # controls the magnitude of dET/deta of frag region
shape-alpha = {alpha} # controls the shape of fragmentation region
shape-beta = {beta} # controls the shape of fragmentation region

# controls the TRENTo normalzation at eta_s=0
# dET/deta[etas=0] = mid-norm  * (sqrt/2M)^mid-power
mid-power = {midpower}
mid-norm = {norm}
flatness = {flatness}

# number of rapidity points of the 3D dET/deta/dx/dy
nsteps-etas = {rapgrid}
"""


label = 'param-set-1'
with open(label+".txt",'w') as f:
    f.write(config.format(
        nparton=3, wnucleon=.6, wparton=.4, 
        variance=.3,
        alpha=4, beta=.4, flatness=1.25
        midpower=0.37, norm=0.4, 
        rapgrid=31, kTmin=.2,
        )
    )


cmd = """trento-3 {proj} {targ} {N} --sqrts {sqrts} -c ./{label}.txt > {filename}.dat
"""

folder = 'Events/'+label+'/'
os.makedirs(folder,exist_ok=True)
for (proj, targ), sqrtslist in zip(
    [('p','p'),
     ('p','Pb'),
     ('Pb','Pb'),
     ('Xe','Xe'),
     ('Au','Au'),
     ('Cu','Cu'),
     ('U','U'),
     ('Cu','Au'),('p','Al'),('p','Au') ,('d','Au') ,('He3','Au')
    ],
    [[53,200,410,546,900,2760,5020,7000],
     [5020, 8160],
     [17.3,2760,5020],
     [5440],
     [200,130,62.4,39,27,19.6,14.5,7.7],
     [200,130,62.4,22.4],
     [193],
     [200],[200],[200],[200],[200]
     ]
    ):
    N = 400
    for sqrts in sqrtslist:
        if proj=='He3':
            Tproj='He3.hdf'
        elif proj == 'Xe':
            Tproj = 'Xe2'
        elif proj == 'U':
            Tproj = 'U3'
        else:
            Tproj=proj
 
        if targ=='He3':
            Ttarg='He3.hdf'
        if targ == 'Xe':
            Ttarg = 'Xe2'
        elif targ == 'U':
            Ttarg = 'U3' 
        else:
            Ttarg = targ
        CMD = cmd.format(proj=Tproj,
                     targ=Ttarg,N=N,sqrts=sqrts,label=label,
                     filename=folder+'/'+proj+'-'+targ+'-{}'.format(sqrts))
        print(CMD)
        sp.call(CMD, shell=True)

