import numpy as np
import sys,pdb
from scipy.special import erfcinv as erfcinv

def gaus_pdf(x,mu=0,sig=1):
    return np.exp(-(x-mu)**2/sig**2/2)/sig/np.sqrt(2*np.pi)

def gauss_func_point(emis_rate=0.,u=0.,x=0.,y=0.,z=0.,xs=0.,ys=0.,H=0.,Dy=1.,Dz=1.):
    x_c=x-xs; # shift the coordinates so that stack is centre point
    y_c=y-ys; 
    uh=u*3600.
    dot_product=u*x_c
    magnitudes=u*np.sqrt(x_c**2+y_c**2)
    subtended=np.arccos(dot_product/(magnitudes+1e-15))
    hypotenuse=np.sqrt(x_c**2+y_c**2)
    downwind=np.cos(subtended)*hypotenuse;
    if downwind < 0.: return 0.
   # Now calculate distance cross wind.
    crosswind=np.sin(subtended)*hypotenuse;

    sig_y=np.sqrt(2.*Dy*downwind/u);
    sig_z=np.sqrt(2.*Dz*downwind/u);
    
    C=emis_rate/(2.*np.pi*uh*sig_y*sig_z) \
       * np.exp(-crosswind**2./(2.*sig_y**2.))  \
       *(np.exp(-(z-H)**2./(2.*sig_z**2.)) + \
       np.exp(-(z+H)**2./(2.*sig_z**2.)) );
   
    return C

def gauss_func(Q,u,x,y,z,xs,ys,H,Dy,Dz):
   x_c=x-xs; # shift the coordinates so that stack is centre point
   y_c=y-ys; 
   uh = u*3600.
   x1,y1 = np.meshgrid(x_c,y_c)
   x1 = x1.flatten()
   y1 = y1.flatten()

   # components of u in x and y directions
   wx=u;
   wy=0;

   # Need angle between point x, y and the wind direction, so use scalar product:
   dot_product=wx*x1+wy*y1;
   # product of magnitude of vectors:
   magnitudes=u*np.sqrt(x1**2.+y1**2.); 

   # angle between wind and point (x,y)
   subtended=np.arccos(dot_product/(magnitudes+1e-15));
   # distance to point x,y from stack
   hypotenuse=np.sqrt(x1**2.+y1**2.);

   # distance along the wind direction to perpendilcular line that intesects
   # x,y
   downwind=np.cos(subtended)*hypotenuse;

   # Now calculate distance cross wind.
   crosswind=np.sin(subtended)*hypotenuse;

   ind=np.where(downwind>0.);

   sig_y=np.sqrt(2.*Dy*downwind/u);
   sig_z=np.sqrt(2.*Dz*downwind/u);
   
   # calculate sigmas based on stability and distance downwind
   #(sig_y,sig_z)=calc_sigmas(STABILITY,downwind);
   Sig_y = np.zeros((len(y),len(x)))
   Sig_z = np.zeros((len(y),len(x)))
   for ix in range(Sig_y.shape[1]):
       Sig_y[:,ix] = sig_y[ix]
       Sig_z[:,ix] = sig_z[ix]
   Sig_y = Sig_y.flatten()
   Sig_z = Sig_z.flatten()

   #pdb.set_trace()
   C=np.zeros((len(y),len(x))).flatten();
   C[ind]=Q/(2.*np.pi*uh*Sig_y[ind]*Sig_z[ind]) \
       * np.exp(-crosswind[ind]**2./(2.*Sig_y[ind]**2.))  \
       *(np.exp(-(z-H)**2./(2.*Sig_z[ind]**2.)) + \
       np.exp(-(z+H)**2./(2.*Sig_z[ind]**2.)) );
   
   return C.reshape((y.shape[0],x.shape[0]))

def gp_1D_solution(source={},atm={},loc={}):
    xi,yi,zi=loc['x'],loc['y'],loc['z']
    source_x,source_y,H,emis_rate = source['xo'],source['yo'],source['zo'],source['emis_rate']
    Dy,Dz,u = atm['Dy'],atm['Dz'],atm['advection_wind_speed']
    mass_density=np.zeros(xi.shape[0])
    for i in range(xi.shape[0]):
        mass_density[i] = gauss_func_point(emis_rate=emis_rate,u=u,x=xi[i],y=yi[i],z=zi[i],xs=0.,ys=0.,H=H,Dy=Dy,Dz=Dz)
    return mass_density*atm['scaling']
    
def gp_3D_solution(source={},atm={},grid={}):
    x,y,z = grid['advect_axis'],grid['crosswind_axis'],grid['vertical_axis']
    source_x,source_y,H,emis_rate = source['xo'],source['yo'],source['zo'],source['emis_rate']
    Dy,Dz,u = atm['Dy'],atm['Dz'],atm['advection_wind_speed']
    mass_density = np.zeros((z.shape[0],y.shape[0],x.shape[0]))
    for iz,zi in enumerate(z[:]):
        mass_density[iz,:,:]=gauss_func(emis_rate,u,x,y,zi,source_x,source_y,H,Dy,Dz);
    return mass_density*atm['scaling']

def apply_is_obs_operator(loc={},atm={},source={}):
    sim_obs = gp_1D_solution(source=source,atm=atm,loc=loc)*atm['scaling']
    return sim_obs

def calculate_cost(observation={},emission={},source={},atm={}):
    
    obs_mf = observation['mole_fraction'][:]
    obs_uncert = observation['uncert']
    
    emis = emission['current']
    prior_emis = emission['prior']
    emis_uncert = emission['prior_uncert']
    
    sim_mf = gp_1D_solution(loc=observation,source=source,atm=atm)

    cost_obs = 0.5*np.dot((obs_mf-sim_mf).T,np.dot(1/obs_uncert**2*np.eye(len(obs_mf)),obs_mf-sim_mf))
    cost_prior = 0.5*(emis-prior_emis)**2/emis_uncert**2
    return cost_obs + cost_prior

def calculate_dcost(observation={},emission={},source={},atm={}):
    obs_mf = observation['mole_fraction'][:]
    obs_uncert = observation['uncert']
    
    emis = emission['current']
    prior_emis = emission['prior']
    emis_uncert = emission['prior_uncert']
    
    sim_mf = gp_1D_solution(loc=observation,source=source,atm=atm)

    dcost_obs = np.dot((obs_mf-sim_mf)/obs_uncert**2,-sim_mf/emis)
    dcost_prior = -(emis-prior_emis)/emis_uncert**2

    return dcost_obs + dcost_prior
    