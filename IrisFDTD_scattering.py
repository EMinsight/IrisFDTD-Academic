"""
Copyright (C) 2019-2020 Sergio G Rodrigo <sergut@unizar.es>

**Important:** 
*IrisFDTD* is licensed under the AGPL and it is free to use. 
if you are using, or plan to use this example, specially if it is 
for research or academic purposes, please send an email with your name, 
institution and a brief description of your interest for this program. 
"""


class fdtdscatt:

 def __init__(self,path,filescatt,filewav,nowav,nout,deltawav,iwav):
    self.nowav = nowav
    self.nout  = nout
    self.deltawav = deltawav
    self.iwav  = iwav
    self.filescatt= filescatt
    self.filewav = filewav
    self.path   = path
    pass
     
 def load_spectra(self):
    import os        
        
    nowav=self.nowav 
    nout= self.nout    
    filescatt=self.filescatt     
    path=self.path    
    
    f = []       
    # File containing scattering coefficients    
    file=os.path.exists(path+filescatt)
    print(path+filescatt)
    if (file):           
       for line in open(path+filescatt, 'r'):
          values = line.split()
          if(values[0] !='#'):
              f.append(values[0])     
    else: 
        print(filescatt+" file doesn't exist!. You have to run first IrisFDTD.exe.")
    
    lst=[]
    # Note: f[0:nowav] may contain wavelengths, like in trans.dat and ref.dat
    for i in range(0,(nout+1)*(nowav)-1,nowav): #range(start, stop, step)
       lst.append((f[i:i+nowav]))
   
    return lst
       
 def get_single_spectra(self,full_spectra,iter):
     
    import numpy as np
    import os
    nowav=self.nowav     
    deltawav=self.deltawav
    iwav=self.iwav    
    filewav=self.filewav
    path=self.path
     
    fwav=iwav+deltawav*nowav # final wavelength in nm (see IrisFDTD)
    
    lst = full_spectra
    wl = lst[0] 
    S = lst[iter]
    map(float, wl) # Convert the list to float numbers
    map(float, S) # Convert the list to float numbers
    wl=np.array(wl)
    S=np.array(S)     
        
    # List of wavelengths if not included in the files can be introduced either
    # from a file
    fileME=os.path.exists(path+filewav)
    print(path+filewav)
    if (fileME):           
       i=0
       for line in open(path+filewav, 'r'):
          values = line.split()
          wl[i]=values[0]
          i=i+1
    else: 
        print(filewav+" File of wavelengths doesn't exist!.")
    
    # When wavelengths are included in the files to plot (trans.dat,ref.dat...)
    x = np.linspace(0, 10, nowav)
    y = np.linspace(0, 10, nowav)    
    for i in range(0,nowav):   
       x[i]=wl[i]
       y[i]=S[i]
       
    # from iwav,fwav,deltawav
    if(deltawav != 0.0): # If deltawav (iwav,fwav also) is provided 
        x = np.linspace(iwav, fwav, nowav)    
    
    return x,y

# Can plot several scattering calculations in the same
# figure.
def plot_scatt(nfiles,x,y,path,pngname,**kwargs):
    import matplotlib.pyplot as plt    
    plt.style.use('ggplot') #('seaborn-whitegrid') 
    
    xmin=kwargs.get('lambdai')
    xmax=kwargs.get('lambdaf')
    
    labels,colors,line = [],[],[]    
    colors = ["red","green","blue"]
    line = ["-","--",".."]
    labels = ["T","R","A"]
    fig, ax = plt.subplots()    
    for i in range(0,nfiles):
        ax.plot(x[i],y[i], color=colors[i],ls=line[i], label=labels[i])
    txt = "$\lambda (nm)$"
    ax.set_xlabel(txt)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-0.01,1.01)
    #ax.set_ylabel("Scattering coefficients")
    ax.legend()
    fig.savefig(path+pngname+".png", dpi=200, facecolor="#f1f1f1")
    plt.show()
pass


def movie_scatt(path,pngname,imgo,imgf,videoname,removeimg):
    import os
    
    print("Creating animation")
    command="ffmpeg -framerate 1 -start_number "+str(imgo)+" -i "+path+pngname+"%d.png "+path+str(videoname)+".avi"
    print(command)
    videoexist=os.path.exists(path+videoname+".avi")
    print("The avi file already exist =",videoexist)
    if(videoexist): 
        os.remove(path+videoname+".avi")
        print("The previous avi file removed")
    ok=os.system(command)    
    if(ok==1):
        ok=False
    else:
        ok=True
    print("A new avi file has been created =",ok)
    if(ok==False):
        print("Check ffmpeg path and dimensions for MPEG-4 (dpi,figx and figy)")
    print("Creating animation done")    

    if(removeimg):
     for l in range(imgo,imgf):       
         fname = path+pngname+str(l)+".png"
         os.remove(fname)
     print("Files output*.png removed")
    else:
     print("Files output*.png NOT removed!!")
           
  
    pass