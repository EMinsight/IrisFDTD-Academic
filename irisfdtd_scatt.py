"""
Copyright (C) 2019-2020 Sergio G Rodrigo <sergut@unizar.es>

**Important:** 
*IrisFDTD* is licensed under the AGPL and it is free to use. 
if you are using, or plan to use this example, specially if it is 
for research or academic purposes, please send an email with your name, 
institution and a brief description of your interest for this program. 
"""

'''
 Plot transmission and reflection obtained by the IrisFDTD program
'''
import numpy as np
import IrisFDTD_scattering as scatt

nowav=500       # number of wavelengths (nof in inputFDTD)
nout=50         # times output files are update (no_downloads in inputFDTD)
deltawav=5.0    # wavelength step in nm (dlambda in inputFDTD)
iwav=400.0      # initinal wavelength in nm (lno in inputFDTD)
filewav = "wavelength_list.dat" # in case wavelengths introduced from a file
path="./example_1/" # Don't forget / at the end    
file1="1_scatt5.dat" # Current in reflection region (z=szi)
file2="1_scatt6.dat" # Current in transmission region (z=szf)
lambdai=400.0;lambdaf=3000.0 # Spectral range shown by plot
pngname= "T-R"
videoname="T-R"

r = scatt.fdtdscatt(path,file1,filewav,nowav,nout,deltawav,iwav)
ref = r.load_spectra() 
t = scatt.fdtdscatt(path,file2,filewav,nowav,nout,deltawav,iwav)
trans = t.load_spectra() 

iter=input("Introduce output step (from 0 to "+str(nout)+") = ")
iter=int(iter)

x,y=r.get_single_spectra(ref,iter)
xx,yy=t.get_single_spectra(trans,iter)

total_curr=y+yy
y=y/total_curr
yy=yy/total_curr

x=np.vstack((xx, x))
y=np.vstack((yy,y))
nofiles=int(x.size/nowav)
print(nofiles)
scatt.plot_scatt(nofiles,x,y,path,pngname,lambdai=lambdai,lambdaf=lambdaf)


