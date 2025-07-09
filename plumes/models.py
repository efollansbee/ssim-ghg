import numpy as np
import sys,pdb
from scipy.special import erfcinv as erfcinv

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
   
def gp_forward_model(source={},atm={},grid={}):
    x,y,z = grid['advect_axis'],grid['crosswind_axis'],grid['vertical_axis']
    source_x,source_y,H,emis_rate = source['xo'],source['yo'],source['zo'],source['emis_rate']
    Dy,Dz,u = atm['Dy'],atm['Dz'],atm['advection_wind_speed']
    mass_density = np.zeros((z.shape[0],y.shape[0],x.shape[0]))
    for iz,zi in enumerate(z[:]):
        mass_density[iz,:,:]=gauss_func(emis_rate,u,x,y,zi,source_x,source_y,H,Dy,Dz);
    return mass_density