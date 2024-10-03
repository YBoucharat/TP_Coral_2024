import matplotlib.pyplot as plt
import numpy as np
from Dicts import Dicos
import matplotlib.cm as mcm
from tools import (FindKey, 
                   readfile)
from tools_models import (ModelParameters,  
                          shore,
                          ZarrName)
from os import path, makedirs
# import matplotlib.colors as colors
from cmcrameri import cm
# import pandas as pd
from SL_MIS import MISage
import xarray as xr
from datetime import datetime as dtime

#********************************************************************************************
#********************************************************************************************

class Profile(plt.Figure):
    """ 
    Figure to plot the cross-section of the simulation with age of construction and presence of clastic sediments 
      
    
    Parameters
    ----------
    ds: xarray dataset
        Simulation out with profile__z and depot__dS at each timestep
        
    **kwargs:
        - width, height: ints, floats, optional
            Figure dimensions
        - xmin, xmax, ymin, ymax: floats, optional
            Profile boundaries. By default, values are set by the model to include 
            any, but only, evolution of the profile
        - tmax: Float, optional
            Time of evolution from the begining of the simulation. 
            By default, tmax is present. Value must be equal to a timestep.
        - store: 'str, optional
            Options to store the final image: 'Profile' (default), 'Anim'
        - fs: float, optional
            Fontsize (14 by default)
            
    Returns
    -------
    A figure environment ready to plot stuff in it
    
    """
    
#*********************************************************************
#*********************************************************************
    
    def __init__(self, 
                 ds,   
                 fs = 14, 
                 width=10, height=7,
                 xmin=0, xmax=0, ymin=0, ymax=0, 
                 store='Profile',
                 tmax = -1
                ):
        
        # Variables initialization
        self.dt = ds.time[1].values-ds.time[0].values   # Simulation timestep
        self.dico = Dicos()
        self.fs = fs
        self.store = store
        
        # Retreiving SL keyword from SL filename
        self.SL = str(ds.SLstory__RSLin.values)        
        # Saving entire SL curve and time for the ASL plot if any
        self.time_ASL = ds.out[1:].values.copy()
        print('time_ASL', self.time_ASL)
        self.ASL = ds.sealevel__asl[1:].values.copy()
        
        # tmax for the figure
        if tmax == -1:
            self.tmax = ds.time[-1].values
            self.ds = ds
#             self.itmax = -1
            
        elif tmax < ds.out[1]:
            self.tmax = ds.out[1].values
            self.ds = ds.isel(out=slice(0, 2), x = slice(None, None)) #, time=slice(0, np.argmax(ds.time.values>=self.tmax)+1))
            
        elif tmax == 0:
            self.ds = ds.isel(time=slice(0, 1), x=slice(None, None), out=slice(0, 1))
            
        else:
            self.ioutmax = np.argmax(ds.out.values == tmax)
            self.tmax = tmax
            # Slicing dataset to tmax 
            self.ds = ds.isel(out=slice(None, self.ioutmax+1), x=slice(None, None)) #, time=slice(0, np.argmax(ds.time.values == self.tmax)+1))
        
        # Defines profile boundaries
        self.AxeBoundaries(xmin, xmax, ymin, ymax)
        # Slicing dataset to profile boundaries
        self.ds = self.ds.isel(out=slice(None, None), x=slice(self.xmin,self.xmax))
        
        # Model parameters to define the simulation
        self.params = ModelParameters(self.ds.model_name)
        
        # Customizing frame
        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
        plt.rcParams['xtick.top']    = plt.rcParams['xtick.labeltop'] = False
        plt.rcParams['ytick.left']   = plt.rcParams['ytick.labelleft'] = True
        plt.rcParams['ytick.right']  = plt.rcParams['ytick.labelright'] = False
        
        # Creating figure
        self.fig, self.ax = plt.subplots(figsize=(width, height), constrained_layout=True)
        
        self.ax.set_xlim(self.xmin, self.xmax)
        self.ax.set_ylim(self.ymin, self.ymax)
        self.ax.set_xlabel('Distance (m)', fontsize=self.fs)
        self.ax.set_ylabel('Elevation (m)', fontsize=self.fs)
        
#*********************************************************************

    def PlotProfile(self,
                    cmap = cm.roma_r, 
                    prof_dt = 1000,
                    dispage = False,
                    hmax = False,
                   ):
        """ 
        Plots cross-section of the simulated profile at t=tmax, with global sediment age (par ky) 
        and final topographic profiles (per ky) by stacking total sediment thickness. 
        Clastic sediment resulting from wave erosion are shaded areas
        
        
        Parameters
        ----------
        cmap: colormap, optional
            Default is Roma from Crameri's package
        prof_dt: Int, Float, optional
            Timestep for the profile (default value = 1ky)
            Must fit simulation timestep
        """

        t0 = dtime.today()
        # Variables initialisation
#         self.cmap = cmap
        # Profile timestep = 1ky or else...
        self.nky = int(prof_dt/self.dt)
        # Vertical deformation through time
        self.ds['dU'] = xr.Variable(data = self.ds.vertical__u.values * self.ds.out.values, dims = ('out'))
        # Displaced profiles
        self.ds['z'] = xr.Variable(data = np.zeros(np.shape(self.ds.profile__z)), dims=('out', 'x'))        
        # Displaced final (with erosion) profiles
        self.ds['zf'] = xr.Variable(data = np.zeros(np.shape(self.ds.z)), dims=('out', 'x'))
        # Total sediment thickness to plot
        self.ds['sed'] = xr.Variable(data = np.zeros((len(self.ds.out), len(self.ds.init__x))), dims = ('out', 'x'))
        # Clastic sediments after erosion
        self.ds['clasts'] = xr.Variable(data = np.zeros((len(self.ds.out), len(self.ds.init__x))), dims = ('out', 'x'))
        
        # Preparing colormap
        self.SL_Colormap(cmap)    
        
        ## Pre-Processing
        # Displacing profiles
        self.ds['z'] = self.ds.profile__z + (self.ds.dU[-1]-self.ds.dU)
        
        # Eroding displaced profiles and clastic sediments
        for t in range(len(self.ds.out)):
            self.ds.zf[t] = self.ds.z[t:].min(dim = 'out')
        
        self.ds.clasts[0:-1] = self.ds.depot__dS[0:-1] + np.ma.masked_where(
                self.ds.zf[1:] < self.ds.z[1:], 
                self.ds.zf[1:] - self.ds.z[1:]
            )
        
        print('Erosion done')

        # Sediment thickness per time step
        self.ds.sed[1:] = np.diff(self.ds.zf, axis=0)
        self.ds.sed[0] = self.ds.zf[0]    
                
        # Stacking by ky
        self.subds = self.ds.isel(out=slice(None, None, self.nky)) #, x=slice(None, None), time=slice(None, None, self.nky))
        if self.subds.out[-1] != self.ds.out[-1]:
            self.subds = xr.concat([self.subds, self.ds.isel(out = -1, x=slice(None, None))], dim = 'out')

        self.subds['sedky'] = xr.Variable(
            data = np.zeros((len(self.subds.out), len(self.ds.init__x))),
            dims = ('out', 'x')
        )
        print('sedky 0', np.shape(self.subds.sedky))
        
        for i in range(1, len(self.subds.out)):
            self.subds.sedky[i,:] = np.sum(self.ds.sed[(i-1)*self.nky+1:i*self.nky+1, :], axis=0)
        
        print('Sediment done')
        # Substrate
        self.subds.sedky[0] = self.ds.sed[0]
        self.subds.color[1] = np.zeros(4)+0.2
        
        # Plot
        # Plot total sediment thickness
        self.ax.stackplot(
            self.ds.init__x.values, 
            self.subds.sedky,
            colors = self.subds.color.values,
            linewidth = 0,
            rasterized = True
        )
        
        # Plot horizons
        self.ax.plot(
            self.subds.init__x.values,
            self.subds.zf.values, 
            color='black', alpha=0.1, lw=0.05, rasterized = True)
        print('Horizons done')
        
        # Plot clastics
        for t in range(1, len(self.ds.out)-1):
#             if np.any(self.ds.depot__dS[t].values!=0):
#                 deb = np.argmax(self.ds.depot__dS[t].values!=0)
#                 fin = len(self.ds.init__x) - np.argmax(self.ds.depot__dS[t, -1::-1].values!=0)
#             else:
#                 print('No sed - t', t, self.ds.out[t].values, self.ds.depot__dS[t].values.sum())

            self.ax.fill_between(
                self.ds.init__x, 
                self.ds.zf[t+1], 
                self.ds.zf[t+1]-self.ds.clasts[t], 
                color='white', alpha=0.6,
                lw = 0,
                rasterized = True,
            )
            
        # Plot substrate and last profile
#         self.ax.plot(self.ds.init__x, self.ds.zf[0], color='black', lw=0.4)
#         self.ax.plot(self.ds.init__x, self.ds.zf[-1], color='red', lw=10)
        
        # Plot SL
        self.ax.hlines(
            y=self.ds.sealevel__asl[-1].values,
            xmin = self.xmin,
            xmax = shore(
                self.ds.sealevel__asl[-1].values, 
                self.ds.profile__z[-1].values)+self.xmin,
            color = 'tab:blue'
        )
        
        # Plot hmax 
        if hmax:
            
            if self.ds.profile__z[-1, 0] >= self.ds.sealevel__asl[-1].values - self.ds.grid__hmax.values:
                xmax_hmax = 0
            else:
                xmax_hmax = shore(
                        self.ds.sealevel__asl[-1].values - self.ds.grid__hmax.values,
                        self.ds.profile__z[-1].values
                        )
            riv = shore(
                    self.ds.sealevel__asl[-1].values, 
                    self.ds.profile__z[-1].values
                )

            y2 = np.zeros(riv+1) + self.ds.sealevel__asl[-1].values - self.ds.grid__hmax.values
            y2[xmax_hmax:] = self.ds.profile__z[-1, xmax_hmax:riv+1]
            y2 = np.maximum(y2, np.zeros(len(y2)) + (self.ds.sealevel__asl[-1].values - self.ds.grid__hmax.values))
        

            self.ax.fill_between(
                np.arange(riv+1) + self.xmin,
                np.zeros(riv +1) + self.ds.sealevel__asl[-1].values,
                y2,
                color='lightblue', alpha = 0.3, lw=0
            )
        
        if dispage:
            self.DispAge()

#*********************************************************************

    def SL_Colormap(self, cmap = cm.roma_r):
        """ 
        Creates the profile colormap based on the type of ASL scenario, 
        MIS style (ASL reconstructions) or cyclic SL curves
        with cold colors for glacial periods
             hot colors for interglacial periods
        """
        
        # Initialisation
        self.cmap = cmap
        self.ds['color'] = xr.Variable(
            data = np.zeros((len(self.ds.out), 4)),
            dims = ('out', 'colo')
        )
        self.ds.color[:] = self.cmap(np.linspace(0, 1, len(self.ds.out)))
        
        # To run MIS based colormap
        if self.SL in ['Holocene', 'Waelbroeck2002', 'Waelbroeck2002b', 'Waelbroeck2002-137ky', 'Waeltanja-1500k', 'Bintanja2008-1000k']:
            self.MIS_Colormap()

        # To run Cyclic colormap
        elif self.SL in ['Sinus-asym']:
            self.Cyclic_Colormap()
            
        else:
            print('WORK !!!')
            
        # Sets the substrate in white
        self.ds.color[0] = np.ones(4)
            
#********************************************************************************************

    def MIS_Colormap(self):
        """ 
        Creates the profile colormap based on ages for MIS boundaries defined in SL_MIS.py
        - Colors are picked in a diverging colormap
        - MIS boundaries are defined for each ASL reconstruction in SL_MIS.py
            or must be added there once for all by someone
        """

        # Variables initialisation
        tmp = MISage(self.SL)
        MIS = tmp.ds
        age_max = MIS.iloc[-1].ages[0] / self.dt
        
        # Creating a diverging colormap
        colors0 = self.cmap(np.linspace(1, 0, int(age_max*2+1)))
        
        # Loop on MIS
        for i in range(len(MIS)-1, -1, -1):  
            print(int(age_max - MIS.ages[i][0]/ self.dt), int(age_max - MIS.ages[i][1] / self.dt))
            self.ds.color[ int(age_max - MIS.ages[i][0]/ self.dt) : int(age_max - MIS.ages[i][1] / self.dt)] = \
                colors0[ int( age_max + MIS.gla[i] * (age_max - (MIS.ages[i][0] + MIS.ages[i][1]) /2 /self.dt)) ]
        
#*********************************************************************

    def AxeBoundaries(self, xmin, xmax, ymin, ymax):
        """ 
        Defines the final axis limits
        - Called automatically at initiation 
        """
    
        if xmin==0:
            self.xmin = int(max(self.ds.profile__xmin.values-100, 0))
        else:
            self.xmin=int(xmin)
        if xmax==0:
            self.xmax = int(min(self.ds.profile__xmax.values+100, self.ds.init__x[-1]))
        else:
            self.xmax = int(xmax)
        if ymin==0:
            self.ymin = int(self.ds.profile__z[0, self.xmin] + self.ds.vertical__u.values * self.tmax)
        else:
            self.ymin = int(ymin)
        if ymax==0:
            self.ymax = int(self.ds.profile__z[0, self.xmax] + self.ds.vertical__u.values * self.tmax+10)
        else:
            self.ymax = int(ymax)
                        
#*********************************************************************

    def PlotASL(self, 
                time = False,
                loc='bottom', 
                s=8):
        """
        Plots ASL reconstruction used in the simulation
        Colored by MIS 
        Needs to run after PlotProfile() because of the colormap
        
        Parameters:
        -----------
        
        - loc: str, optional
            Sets the location ('top', 'bottom') of the inset. Default is 'bottom'
        - s: int, float, optional
            Size of markers, corresponds to line thickness. 
            Default is 8.
        """
        
        # Set up inset location
        if loc == 'bottom':
            plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
            plt.rcParams['xtick.top']    = plt.rcParams['xtick.labeltop'] = True

            axSL= self.ax.inset_axes([0.65, 0., 0.35, 0.25])
            secax = axSL.secondary_xaxis('top')
            secax.set_xlabel('Age (ka)')
            axSL.set_ylabel('Elevation (m)')
            
        elif loc == 'top':
            plt.rcParams['ytick.left']  = plt.rcParams['ytick.labelleft'] = False
            plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
            
            axSL= self.ax.inset_axes([0, 0.75, 0.35, 0.25])
            secax = axSL.secondary_yaxis('right')
            secax.set_ylabel('Elevation (m)')
            axSL.set_xlabel('Age (ka)')
            
        # Plots ASL
        axSL.scatter(
            (self.time_ASL[::-1]-self.ds.out[0].values)/1000, 
            self.ASL,
            c=self.ds.color[1:], s=s, marker='o') 
                    
        # Plot tmax, for animations
        if time:
            axSL.plot(
                np.ones(2) * (self.time_ASL[-1] - self.ds.time[-1].values)/1000, 
                np.linspace(-140, 15, 2), 
                color='darkgrey'
            )
            
        axSL.set_xlim(0, self.ds.out[0]/1000+1)
        axSL.set_ylim(-140, 15)
        
        axSL.set_title(self.SL)
            
#*********************************************************************

    def Savefig(self, dpi=300):
        """ 
        Saves final figure
        - Automatically creates output filename
        - Creates the directory if needed 
        - Outputs directory can be modified in Dicts.direOuts 
        
        Parameters:
        -----------
        
        - ext: str, optional
            Extension of output image filename
            Default is .png
        - dpi: int, float, optional
            Definition of output image
            Default is 300
        """
        
        # Path and filename for image output
        if self.namepath in globals() == False:
            self.ProfileName()
        self.fig.savefig(self.namepath+self.namefig, dpi=dpi)

#*********************************************************************

    def ProfileName(self, ext='.png', anim=False, savedir='default'):
        """ 
        Defines the file name and path for Profile output
        Creates the directory if it doesn't exist 
        """

        # Basis for the image and directory (if needed) names
        self.namefig = ('Profile-'+self.ds.model_name+'_'+self.SL)
        extime='-tmax'+str(self.tmax/1000)+'ky-dt'+str(int(self.dt))+'y'

        # Loop on parameters to be added to the filename
        for p in self.params:
            self.namefig += '-' +  \
                self.dico.abbrev[p]+ \
                str(np.round(
                    self.ds[p].values*self.dico.factors[p], 
                    self.dico.rounds[p]-int(np.log10(self.dico.factors[p]))
                    )
                )
        self.namefig += extime     

        # Path of outputs directory
        if savedir == 'default':
            if anim:
                self.namepath='/Users/pastier/REEF/Anim/'+self.namefig+'/'
            else:
                self.namepath = self.dico.direOuts[self.ds.model_name]   
        else:
            self.namepath = savedir
            
        # Create directory if needed
        if not path.exists(self.namepath):
            makedirs(self.namepath)
            
        self.namefig += ext


#***************************************************************************************    

    def WriteParams(self, 
                    xtext=0.22,    # Text location along x-axis
                    ytext=0.17,    # Text location along y-axis 
                    ncol=2,        # Number of columns to write parameters
                    npmax = 0,      # Number of parameters per column
                    dx=0.17,       # Space between text colums
                    dy=0.038,      # Space between text lines
                    fs=0           # Fontsize
                   ):
        """ 
        Writes parameters values on the profile 
        
        Parameters
        ----------
        - xtext, ytext: Floats, optional
            Location of the first column and line of text
            Input values are proportion of the figure size
        - dx, dy: Floats, optional
            Space between columns/lines
            Input values are proportion of the figure size
        - ncol: Int, float, optional
            Number of colums to write parameters (2 by default)
        - fs: Float, optional
            Fontsize (based on labels by default)
            
        """
            
        # Initialisation
        if fs == 0:
            self.fsp = self.fs - 2
        else:
            self.fsp = fs
        # Lines of parameters in each column
        if npmax == 0:
            npmax = int(self.dico.nparams[self.ds.model_name])
        else:
            npmax = int(npmax)
        ddy = 0
        
        # Loop on parameters
        for n_p in range(len(self.params)):
            p = self.params[n_p]
            self.ax.text(
                x=xtext, 
                y=ytext - ddy, 
                s=self.dico.abbrev[p]+' = '+str(
                         np.round(
                             self.ds[p].values * self.dico.factors[p], 
                             self.dico.rounds[p] - int(np.log10(self.dico.factors[p]))
                         )
                         ) + self.dico.units[p],
                fontsize=self.fsp,
                transform = self.ax.transAxes
            )
            
            if ((n_p+1) % npmax == 0) & (n_p > 0):
                xtext += dx                                                   
                ddy = 0
            else:
                ddy += dy

#***************************************************************************************    

    def DispAge(self):
        """ 
        Displays the duration of plotted simulation as figure title
        - default is no 
        """
        
        self.ax.set_title('t = '+str(self.tmax/1000)+'ky', fontsize=self.fs+2)
        
#*********************************************************************

    def PlotWidth(self, dsl = 10,            # Height above SL for reef width
                  sw = 20,                   # Size of markers
                  ywid = 0.92):              # Location of text

        """ Plots reef crest, inner limit of a barrier and shore, and writes widths values """

        deb, rwid, twid, typ = ReefWidths(2, self.ds.sealevel__asl[-1].values, self.ds.profile__z[-1].values, self.ds.grid__spacing.values)
        xwid = (deb-self.xmin)/(self.xmax-self.xmin)
        
        # Plotting reef crest
        self.ax.scatter(
            self.ds.init__x[deb], dsl, c='darkorange', s=sw)
        # Plotting shore
        self.ax.scatter(
            self.ds.init__x[shore(self.ds.sealevel__asl[-1].values, self.ds.profile__z[-1].values)], dsl, 
            c='darkorange', s =sw
        )
        # Writing total width
        self.ax.text(
            x=xwid, y=ywid, 
            s='Total width = '+str(twid)+' m',
            transform=self.ax.transAxes, fontsize=self.fs, ha='center')

        if typ == 'bar':
        # Plotting inner reef limit
            self.ax.scatter(
                self.ds.init__x[deb+int(rwid)], dsl, 
                c='darkorange', s=sw)

            self.ax.text(
                x=xwid, y=ywid-0.09,
                s='Reef width = '+str(rwid)+' m',
                fontsize = self.fs, transform=self.ax.transAxes, ha='center')

#*********************************************************************

    def PlotData(self, namefile, color='black', lw=1, dx=0):
        """ 
        Plots a topographic profile to compare 
        
        PARAMETERS:
        -----------
        
        namefile: str
        Path and filename for the file containing the data to plot 
        Data are x and z columns
        
        kwargs:
        - color: str or RGB value...
        - lw : float 
            Linewidth
        - dx : float
            Distance to shift data horizontally
        
        """
        
        # Read data
        x, z = readfile(namefile)
        
        # Shift data to simulation profile
        x0 = dx+np.argmax(self.ds.profile__z[-1].values>=z[0])+self.xmin
        self.ax.plot(x+x0, z, color=color, linewidth=lw)
        
#*********************************************************************

    def PlotLines(self, heights, lw=1, color='grey'):
        """ 
        Plots a series of horizontal lines
        
        PARAMETERS:
        -----------
        
        heights: list or array of floats
            Elevations to plot
            
        kwargs:
        - color: str or RGB value...
        - lw : float 
            Linewidth
        """
        
        for h in heights:
            self.ax.hlines(y=h, xmin=self.xmin, xmax=self.xmin+shore(h, self.ds.profile__z[-1].values), color=color, linewidth = lw)