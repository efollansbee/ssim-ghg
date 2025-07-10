import numpy as np
import sys,pdb
import yaml

def read_settings(file='../site_settings.ini'):
    #--------------------------------------------------------------------------
    with open(file,'r') as file:
        _s = yaml.safe_load(file)
        settings = _s['plumes']
        global_paths = _s['global_paths']
    data_root = global_paths['input_folder']

    # Spatial Grid
    dxy=settings['horizontal_resolution'];          # resolution of the model in both x and y directions
    dz=settings['vertical_resolution'];
    ztop=settings['ztop']
    x=np.mgrid[-100:5000+dxy:dxy]
    y=np.mgrid[-2500:2500+dxy:dxy]
    z=np.mgrid[0:ztop+dz:dz]
    
    grid={}
    grid['advect_axis']=x; # solve on a 5 km domain
    grid['crosswind_axis']=y;              # x-grid is same as y-grid
    grid['vertical_axis']=z;
    grid['dxy']=dxy
    grid['dz']=dz

    #--------------------------------------------------------------------------
    # Thermodynamic variables
    ps = settings['surface_pressure'] # dry surface pressure
    dp = -dz/1000.;  # approximate drop in pressure with altitude (Pa) per 10m
    p = np.arange(ps,ps+grid['vertical_axis'].shape[0]*dp,dp)
    grid['p'] = p
    grid['dp'] = dp
    T = settings['T']
    rho_a = p/(287.058*T)
    rho_a = rho_a.mean() #density of dry air near surface

    #--------------------------------------------------------------------------
    # Gas Constants
    Ms=16e-3; #molar mass of CH4
    Ma=29e-3; #molar mass of Dry Air

    #--------------------------------------------------------------------------
    # Emissions
    source = {}
    source['xo']=settings['source_x'];
    source['yo']=settings['source_y'];
    source['zo']=settings['source_z']; # stack height, m
    source['emis_rate']=settings['source_rate']; # mass emitted per unit time

    #--------------------------------------------------------------------------
    # Meteorological Fields
    atm={}
    atm['advection_wind_speed']=settings['wind_speed']; # m/s
    atm['Dy']=settings['Dy'];
    atm['Dz']=settings['Dz'];
    atm['scaling']= Ms/Ma / rho_a * 1e6
    atm['background'] = settings['mf_bckgd']

    #--------------------------------------------------------------------------
    # Instrument
    obs={}
    obs['col_noise'] = 5.0 #percent of bckgd
    obs['is_noise'] = 0.05 #percent of bckgd

    return grid,source,atm,obs
