################################################################################
#
# Orbital library for gnuplot
#   - Version     :  2.12 beta
      HSG_libver  = "2.12 beta";
#   - Author      : Helton da Silva Gaspar
#   - Mail        : helton.s.gaspar@ufsc.br
#   - Date        : 04 / 06 / 2021
      HSG_libdate = "2021-06-07";
#
################################################################################

################################################################################
# 
# DESCRIPTION
# 
#   This library gathers useful functions for plotting orbits. It contains 
#   four main packages:
#     1. rv2** whose functions convert fixed cartesian coordinates to 
#              keplerian elements 
#     2. eo2** whose functions convert keplerian elements to fixed cartesian 
#              coordinates
#     3. polar and perifocal frame functions
#     4. time/position functions
#     5. constants
#
# DEFINITIONS AND UNIFORMITY
#    p   stands for the semilatus rectum
#    a   stands for the semimajor axis ( a==p case e==1 )
#    e   stands for the eccentricity
#    I   stands for the inclination
#    O   stands for the longitude of the ascending node
#    w   stands for the argument of periapsis
#    f   stands for the true anomaly
#    M   stands for the mean anomaly
#    E   stands for the eccetric anomaly
#    n   stands for the mean motion
#    T   stands for the period
#    t   stands for the flight time since periapsis
#    C 3 stands for the specific mechanical energy
#    mu  stands for the mass parameter (G*mass)
#
# FUNCTIONS WITHIN THE rv2** PACKAGE
#   In order to keep uniformity, all the package functions must have the 
#   following parameters as argument: mu, x, y, z, vx, vy and vz
#     rv2**(mu,x,y,z,vx,vy,vz)
#   despite not all of them are necessary for some of the functions.
#   
#   SYNOPSIS rv2**(mu,x,y,z,vx,vy,vz), where ** may be one of the following:
#               ** =  p,a,e,I,O,w,f,n,T,C3
#
# FUNCTIONS WITHIN THE eo2** PACKAGE
#   In order to keep uniformity, all the package functions must have the 
#   following parameters as argument: mu, a, e, I, O, w and  f.
#     rv2**(mu,a,e,I,O,w,f)
#   despite not all of them are necessary for some of the functions.
#   
#   SYNOPSIS eo2**(mu,a,e,I,O,w,f), where ** may be one of the following:
#               ** = x,y,z,vx,vy,vz,n,T,t
#   
# FUNCTIONS WITHIN THE POLAR PACKAGE
#   - polp(a,e)      returns semilatus rectum     given a and e
#   - polr(a,e,f)    returns range from focus     given a, e and f
#   - poldrdf(a,e,f) returns the derivative dr/dr given a, e and f
#   
# FUNCTIONS OF THE PERIFOCAL FRAME PACKAGE
#   - x_pq(a,e,w,f)     returns the x coordinate within the perifocal frame
#   - y_pq(a,e,w,f)     returns the y coordinate within the perifocal frame
#   - vr(mu,a,e,f)      returns the radial      component of the velocity
#   - vt(mu,a,e,f)      returns the transversal component of the velocity
#   - vx_pq(mu,a,e,w,f) returns the x           component of the velocity
#   - vy_pq(mu,a,e,w,f) returns the y           component of the velocity
#   
# FUNCTIONS WITHIN THE TIME/POSITION PACKAGE
#   - f2E(e,f)    returns eccetric anomaly            given e and f
#   - E2M(e,E)    returns mean     anomaly            given e and E
#   - f2M(e,f)    returns mean     anomaly            given e and f
#   - M2E(e,M)    returns eccetric anomaly            given e and M
#   - E2f(e,E)    returns true     anomaly            given e and E
#   - M2f(e,M)    returns true     anomay             given e and M
#   - t2E(n,e,t)  returns eccetric anomaly            given n, e and t
#   - t2f(n,e,t)  returns true     anomaly            given n, e and t
#   - M2t(n,M)    returns flight time since periapsis given n and M
#   - f2t(n,e,f)  returns flight time since periapsis given n, e and f
#   - t2M(n,t)    returns flight time since periapsis given n and t
# 
################################################################################



################################################################################
#
## Package constants

  liborb_nrmax = 100;
  liborb_nrdE  = 0.0;
  liborb_zero  = 1.000000000000000000e-08;
  liborb_deg   = 1.745329251994329547e-02;
  liborb_hpi   = 1.570796326794896558e+00;
  liborb_dpi   = 6.283185307179586232e+00;
  liborb_au    = 1.495978707000000000e+11;  #m
  liborb_G     = 6.674080000000000335e-11;  #m^3/kg s^2
  liborb_musun = 1.327541252800000000e+20;  #m^3/s^2
  liborb_mumer = 2.203900081014675781e+13;  #m^3/s^2
  liborb_muven = 3.249606376090141250e+14;  #m^3/s^2
  liborb_muear = 3.987256375298450625e+14;  #m^3/s^2
  liborb_mumoon = 4.904338898288770508e+12; #m^3/s^2
  liborb_mumar = 4.284176672342150781e+13;  #m^3/s^2
  liborb_mujup = 1.267525686099165120e+17;  #m^3/s^2
  liborb_musat = 3.795254329314348000e+16;  #m^3/s^2
  liborb_muura = 5.796369087341472000e+15;  #m^3/s^2
  liborb_munep = 6.838681433981858000e+15;  #m^3/s^2
  
  # Alias to short names of the constants
  deg = liborb_deg;
  dpi = liborb_dpi;
  hpi = liborb_hpi;



################################################################################
#
## math sub-pack

  liborb_mag2(Ux,Uy,Uz) = Ux*Ux + Uy*Uy + Uz*Uz;
  liborb_mag (Ux,Uy,Uz) = sqrt( liborb_mag2(Ux,Uy,Uz) );

  liborb_wrap(x,xmin,xmax)   = x<xmin\
                             ? liborb_wrap(x+xmax-xmin,xmin,xmax)\
                             : x>xmax\
                               ? liborb_wrap(x-xmax+xmin,xmin,xmax)\
                               : x;



  # linear interpolation of the points (0,a) (1,b)
  lin(a,b,u)          = (1-u)*a+u*b;

  # Euler angles rotations
  liborb_Rx(x,y,z,a,b,c)     = x*( cos(a)*cos(c)-sin(a)*cos(b)*sin(c) ) - y*( cos(a)*sin(c) + sin(a)*cos(b)*cos(c) ) + z*sin(a)*sin(b);
  liborb_Ry(x,y,z,a,b,c)     = x*( sin(a)*cos(c)+cos(a)*cos(b)*sin(c) ) - y*( sin(a)*sin(c) - cos(a)*cos(b)*cos(c) ) - z*cos(a)*sin(b);
  liborb_Rz(x,y,z,a,b,c)     = x*  sin(b)*sin(c)                        + y*  sin(b)*cos(c)                          + z*cos(b);

  # Alias to short name functions
  mag2(Ux,Uy,Uz)             = liborb_mag2(Ux,Uy,Uz);
  mag (Ux,Uy,Uz)             = liborb_mag (Ux,Uy,Uz);
  wrap(x,xmin,xmax)          = liborb_wrap(x,xmin,xmax);


################################################################################
#
## conversion from "fixed" cartesian coordinates to keplerian elements

  hx(x,y,z,vx,vy,vz)       = y*vz-z*vy;
  hy(x,y,z,vx,vy,vz)       = z*vx-x*vz;
  hz(x,y,z,vx,vy,vz)       = x*vy-y*vx;
  h2(x,y,z,vx,vy,vz)       = liborb_mag2( y*vz-z*vy, z*vx-x*vz, x*vy-y*vx );

  rv2C3(mu,x,y,z,vx,vy,vz) = 0.5*liborb_mag2(vx,vy,vz) - mu*1.0/liborb_mag(x,y,z);
  rv2p (mu,x,y,z,vx,vy,vz) = h2(x,y,z,vx,vy,vz)*1.0/mu;
  rv2a (mu,x,y,z,vx,vy,vz) = ( aux1=rv2C3(mu,x,y,z,vx,vy,vz) )**2 < liborb_zero**2\
                           ? rv2p(mu,x,y,z,vx,vy,vz)\
                           : abs(0.5*mu/aux1);
  rv2e (mu,x,y,z,vx,vy,vz) = ( aux1=rv2C3(mu,x,y,z,vx,vy,vz) ) == 0.0\
                           ? 1.0\
                           : sqrt(2.0*aux1*rv2p(mu,x,y,z,vx,vy,vz)/mu+1.0);
  rv2I (mu,x,y,z,vx,vy,vz) = ( aux1=liborb_mag(hx(x,y,z,vx,vy,vz),hy(x,y,z,vx,vy,vz),0.0) ) == 0.0\
                           ? atan2(0,hz(x,y,z,vx,vy,vz))\
                           : atan2(aux1,hz(x,y,z,vx,vy,vz));
  rv2O (mu,x,y,z,vx,vy,vz) = (z+vz)==0\
                           ? 0.0\
                           : liborb_wrap( atan2(+hx(x,y,z,vx,vy,vz),-hy(x,y,z,vx,vy,vz)),0.0,liborb_dpi );
  rv2f (mu,x,y,z,vx,vy,vz) = (aux1=rv2e(mu,x,y,z,vx,vy,vz))>1e-6\
                           ? acos( (rv2p(mu,x,y,z,vx,vy,vz)/liborb_mag(x,y,z)-1.0)/aux1)*((x*vx+y*vy+z*vz)<0.0?-1:+1)\
                           : 0*(aux3=abs(z+vz)<liborb_zero?1.0:-hy(x,y,z,vx,vy,vz))\
                              *(aux4=abs(z+vz)<liborb_zero?0.0:+hx(x,y,z,vx,vy,vz))\
                              *(aux5=-hz(x,y,z,vx,vy,vz)*aux4)\
                              *(aux6=+hz(x,y,z,vx,vy,vz)*aux3)\
                              *(aux7=liborb_mag2(hx(x,y,z,vx,vy,vz),hy(x,y,z,vx,vy,vz),0.0))\
                           +  atan2((x*aux5+y*aux6+z*aux7)*1.0/liborb_mag(aux5,aux6,aux7),(x*aux3+y*aux4)*1.0/liborb_mag(aux3,aux4,0.0));
  rv2w (mu,x,y,z,vx,vy,vz) = (aux1=rv2e(mu,x,y,z,vx,vy,vz))==0\
                           ? 0.0\
                           : 0*(aux3=(z+vz)==0?1.0:-hy(x,y,z,vx,vy,vz))\
                              *(aux4=(z+vz)==0?0.0:+hx(x,y,z,vx,vy,vz))\
                              *(aux5=-hz(x,y,z,vx,vy,vz)*aux4)\
                              *(aux6=+hz(x,y,z,vx,vy,vz)*aux3)\
                              *(aux7=liborb_mag2(hx(x,y,z,vx,vy,vz),hy(x,y,z,vx,vy,vz),0.0))\
                           +  liborb_wrap( atan2((x*aux5+y*aux6+z*aux7)*1.0/liborb_mag(aux5,aux6,aux7),(x*aux3+y*aux4)*1.0/liborb_mag(aux3,aux4,0.0))-rv2f(mu,x,y,z,vx,vy,vz),0.0,liborb_dpi);
  rv2n (mu,x,y,z,vx,vy,vz) = 0*(aux1=rv2e(mu,x,y,z,vx,vy,vz))\
                              *(aux2=aux1==1.0?1.0:abs(aux1**2-1))\
                           +  mu**2*(aux2/h2(x,y,z,vx,vy,vz))**(1.5);
  rv2T (mu,x,y,z,vx,vy,vz) = liborb_dpi/rv2n (mu,x,y,z,vx,vy,vz);



################################################################################
#
## polar functions

  polp(a,e)      = e==1.0\
                 ? a\
                 : a*abs(1.0-e*e);
  polr(a,e,f)    = (1.0+e*cos(f))==0\
                 ? 1.0e1000\
                 : polp(a,e)*1.0/(1.0+e*cos(f));
  poldrdf(a,e,f) = polp(a,e)*e*sin(f)/(1.0+e*cos(f))**2;



################################################################################
#
## perifocal frame

  x_pq(a,e,w,f)     = polr(a,e,f)*cos(w+f);
  y_pq(a,e,w,f)     = polr(a,e,f)*sin(w+f);
  vr(mu,a,e,f)      = sqrt(mu/polp(a,e))*e*sin(f);
  vt(mu,a,e,f)      = sqrt(mu/polp(a,e))*(1.0+e*cos(f));
  vx_pq(mu,a,e,w,f) = vr(mu,a,e,f)*cos(w+f) - vt(mu,a,e,f)*sin(w+f);
  vy_pq(mu,a,e,w,f) = vr(mu,a,e,f)*sin(w+f) + vt(mu,a,e,f)*cos(w+f);
  flightangle(e,f)  = atan2(e*sin(f),1.0+e*cos(f));


################################################################################
#
## conversion from keplerian elements to "fixed" cartesian coordinates

  eo2x (mu,a,e,I,O,w,f) = liborb_Rx( x_pq(   a,e,0.0,f), y_pq(   a,e,0.0,f),0.0,O,I,w);
  eo2y (mu,a,e,I,O,w,f) = liborb_Ry( x_pq(   a,e,0.0,f), y_pq(   a,e,0.0,f),0.0,O,I,w);
  eo2z (mu,a,e,I,O,w,f) = liborb_Rz( x_pq(   a,e,0.0,f), y_pq(   a,e,0.0,f),0.0,O,I,w);
  eo2vx(mu,a,e,I,O,w,f) = liborb_Rx(vx_pq(mu,a,e,0.0,f),vy_pq(mu,a,e,0.0,f),0.0,O,I,w);
  eo2vy(mu,a,e,I,O,w,f) = liborb_Ry(vx_pq(mu,a,e,0.0,f),vy_pq(mu,a,e,0.0,f),0.0,O,I,w);
  eo2vz(mu,a,e,I,O,w,f) = liborb_Rz(vx_pq(mu,a,e,0.0,f),vy_pq(mu,a,e,0.0,f),0.0,O,I,w);

  eo2n (mu,a,e,I,O,w,f) = sqrt( mu/a**3 );
  eo2T (mu,a,e,I,O,w,f) = liborb_dpi/eo2n(mu,a,e,I,O,w,f);

  # Alias in english (e)lementos (o)rbitais -> (o)rbital (e)lements
  oe2x (mu,a,e,I,O,w,f) = eo2x (mu,a,e,I,O,w,f);
  oe2y (mu,a,e,I,O,w,f) = eo2y (mu,a,e,I,O,w,f);
  oe2z (mu,a,e,I,O,w,f) = eo2z (mu,a,e,I,O,w,f);
  oe2vx(mu,a,e,I,O,w,f) = eo2vx(mu,a,e,I,O,w,f);
  oe2vy(mu,a,e,I,O,w,f) = eo2vy(mu,a,e,I,O,w,f);
  oe2vz(mu,a,e,I,O,w,f) = eo2vz(mu,a,e,I,O,w,f);
  oe2n (mu,a,e,I,O,w,f) = eo2n (mu,a,e,I,O,w,f);
  oe2T (mu,a,e,I,O,w,f) = eo2T (mu,a,e,I,O,w,f);



################################################################################
#
## time and position functions

  f2E(e,f)              = e>1.0\
                        ? 2.0*atanh( sqrt((e-1)/(e+1))*tan(0.5*f) )\
                        : 2.0*atan2( sqrt(1-e)*sin(0.5*f),sqrt(1+e)*cos(0.5*f) );
  E2M(e,E)              = e>1.0\
                        ? e*sinh(E)-E\
                        : E-e*sin(E);
  f2M(e,f)              = e==1.0\
                        ? 0.5*(aux1=tan(0.5*f))*( 1.0 + 0.33333333*aux1*aux1 )\
                        : E2M(e,f2E(e,f));
  M2t(n,M)              = M*1.0/n;
  f2t(n,e,f)            = M2t(n,f2M(e,f));
  eo2t(mu,a,e,I,O,w,f)  = M2t(eo2n(mu,a,e,I,O,w,f),f2M(e,f));

  newtonraph(e,M,E0,itr)= 0*(dE=e<1.0\
                             ? (E0-e*sin(E0)-M)/(1.0-e*cos(E0))\
                             : (e*sinh(E0)-E0-M)/(e*cosh(E0)-1.0))\
                        + itr>0\
                          ? dE**2>liborb_zero**2\
                            ? newtonraph(e,M,E0-dE,itr-1)\
                            : E0-dE\
                          : 0*(liborb_nrdE=dE)\
                            + E0-dE;

  t2M(n,t)              = n*t;
  M2E(e,M)              = e==1.0\
                        ? 1/0\
                        : newtonraph(e,M,M,liborb_nrmax);
  M2Err(e,M)            = e==1.0\
                        ? 1/0\
                        : M-E2M(e,newtonraph(e,M,M,liborb_nrmax));
  E2f(e,E)              = e>1.0\
                        ? 2.0*atan2( sqrt(1+e)*sinh(0.5*E) , sqrt(e-1)*cosh(0.5*E) )\
                        : 2.0*atan2( sqrt(1+e)*sin( 0.5*E) , sqrt(1-e)*cos( 0.5*E) );
  M2f(e,M)              = e==1.0\
                        ? 2.0*atan( (aux1=3.0*M+sqrt(9.0*M*M+1))**(1.0/3)-aux1**(-1.0/3) )\
                        : E2f(e,M2E(e,M));
  t2E(n,e,t)            = M2E(e,t2M(n,t));
  t2f(n,e,t)            = M2f(e,t2M(n,t));



################################################################################
#
## Help "functions"

  more_orbital_anomalies=sprintf("\
  Anomalies function are:\n\
      f2E(e,f)    returns eccentric anomaly            given e and f\n\
      E2M(e,E)    returns mean      anomaly            given e and E\n\
      f2M(e,f)    returns mean      anomaly            given e and f\n\
      M2E(e,M)    returns eccentric anomaly            given e and M\n\
      E2f(e,E)    returns true      anomaly            given e and E\n\
      M2f(e,M)    returns true      anomay             given e and M\n\
      t2E(n,e,t)  returns eccentric anomaly            given n, e and t\n\
      t2f(n,e,t)  returns true      anomaly            given n, e and t\n\
      M2t(n,M)    returns flight time since periapsis given n and M\n\
      f2t(n,e,f)  returns flight time since periapsis given n, e and f\n\
      t2M(n,t)    returns flight time since periapsis given n and t\n\n\
  Please, type:\n\
     print more_orbital_nr\n\
  to learn more about accuracy in Newton-Raphson method\n");

  more_orbital_nr=sprintf("\
  Newton-Raphson iteration methods are employed to solve Kepler equations:\n\
     M = E - e*sin(E)\n\
     M = e*sinh(E) - E\n\
  for eccentric anomaly E, for eliptic and hiperbolic cases, respectively.\n\
  Convergence may fail depending upon the value of the eccentricity. The closer\n\
  to 1 is the eccentricity the harder is the convergence. Iterations are carried\n\
  out until:\n\
     | E_n - E_{n-1} | < liborb_zero (%.0e by default)\n\
  or\n\
     number of iterations = liborb_nrmax (%d by default)\n\n\
  When recursion limit liborb_nrmax is reached, the \"error\" is saved as:\n\
     liborb_nrdE = E_n - E_{n-1}\n\n\
  The library also provides the function:\n\
     M2Err(e,M)\n\
  which allows one to infer the resulting error of Newton-Raphson. Such error is the\n\
  difference of the resulting E converted back to M' to the original value of M.\n", liborb_zero, liborb_nrmax );

  more_orbital_perifocal=sprintf("\
  2D perifocal system of reference:\n\
      x_pq(a,e,w,f)     returns the x coordinate within the perifocal frame\n\
      y_pq(a,e,w,f)     returns the y coordinate within the perifocal frame\n\
      vr(mu,a,e,f)      returns the radial      component of the velocity\n\
      vt(mu,a,e,f)      returns the transversal component of the velocity\n\
      vx_pq(mu,a,e,w,f) returns the x           component of the velocity\n\
      vy_pq(mu,a,e,w,f) returns the y           component of the velocity\n\
      flightangle(e,f)  returns the flight path angle (gamma = atan(vr/vt) )\n");



  more_orbital_eo2rv=sprintf("\
  General conversion FROM orbital elements TO \"fixed\" cartesian frame of reference: (eo2...)\n\
      eo2**(mu,a,e,I,O,w,f) returns ** given the mu, a, e, I, O, w and f\n\
                              where ** may be x, y, z, vx, vy, vz, n, T or t\n\
  In order to keep uniformity, all the eo2** functions must have the seven parameters\n\
    mu, a, e, I, O, w and f\n\
  as argument despite not all of them are necessary for some of the functions. Type:\n\
    print more_orbital_definitions\n\
  for further information\n");



  more_orbital_rv2eo=sprintf("\
  General conversion FROM \"fixed\" cartesian coordinates TO orbital elements\n\
      rv2**(mu,x,y,z,vx,vy,vz) returns ** given the mu, x, y, z, vx, vy and vz\n\
                              where ** may be p, a, e, I, O, w, f, n, T or t\n\
  In order to keep uniformity, all the rv2** functions must have the seven parameters\n\
    mu, x, y, z, vx, vy and vz\n\
  as argument despite not all of them are necessary for some of the functions. Type:\n\
    print more_orbital_definitions\n\
  for further information\n");



  more_orbital_definitions=sprintf("\
  Definitions and uniformity:\n\
   - p  semilatus rectum\n\
   - a  semimajor axis. In particular case parabolic orbits \"a\" stands for the\n\
        semilaltus rectum while for the circular ones, the circumference radius.\n\
   - e  eccentricity\n\
   - I  inclination of the orbital plane\n\
   - O  longitude of the ascending node\n\
   - w  argument of periapsis\n\
   - f  true anomaly\n\
   - M  mean anomaly\n\
   - E  eccetric anomaly\n\
   - n  mean motion\n\
   - T  period\n\
   - t  flight time since periapsis\n\
   - C3 specific mechanical energy\n\
   - mu mass parameter (G*mass)\n");



  more_orbital_constants=sprintf("\
  Builtin constants are:\n\
      liborb_zero   is the threshold value considered as zero within some\n\
                    functions such as Newton-Raphson iteration or the planar\n\
                    orbits constrained to |z+vz|<liborb_zero (%.0e by default)\n\
      liborb_nrmax  is the maximum number of iterations performed in Newton-Raphson\n\
                    method (%d by default).\n\
      liborb_nrdE   the resulting error of Newton-Raphson when it doesn't converge\n\
                    after %d iterations\n\
      liborb_au     is the astronomical unit (%.6e m by default)\n\
      liborb_G      is the Gravitational constant (%.6e m^3/kg s^2 by default)\n\
      liborb_muxxx  is G times mass for xxx planet where:\n\
                      xxx = (sun, mer, ven, ear, moon, mar, jup, sat, ura, nep, plu)\n\
      liborb_dpi    is (d)ouble pi\n\
      liborb_hpi    is (h)alf pi\n\
      liborb_deg    is 1 (deg)ree in radians\n\
  ",liborb_zero, liborb_nrmax, liborb_nrmax, liborb_au, liborb_G )



  more_orbital=sprintf("\n Orbital library %s (%s) for gnuplot \n\
  This pack gathers a set of functions for orbital motion plotting.\n\n\
   Please, type:\n\
    print more_orbital_all       #to see ALL the functions\n\
    print more_orbital_anomalies #to see the Anomalies functions\n\
    print more_orbital_perifocal #to see the Perifocal frame functions\n\
    print more_orbital_eo2rv     #to see the Keplerian to Cartesian conversion functions\n\
    print more_orbital_rv2eo     #to see the Cartesian to Keplerian conversion functions\n\
    print more_orbital_constants #to see the the builtin constants\n\
  \n\
  print more_orbital_all",\
  HSG_libver, HSG_libdate);



  more_orbital_all=sprintf("\n Orbital library %s (%s) for gnuplot \n\
     This pack gathers a set of functions for orbital motion.\n\n\
  %s\n\
  %s\n\
  %s\n\
  %s\n\
  %s\n",\
  HSG_libver, HSG_libdate,\
  more_orbital_anomalies,\
  more_orbital_perifocal,\
  more_orbital_eo2rv,\
  more_orbital_rv2eo,\
  more_orbital_constants)



################################################################################
#
## Package Header

  print sprintf("\n        HSG_lib orbital");
  print sprintf("        Version %-20slast modified %-10s\n", HSG_libver, HSG_libdate);
  print sprintf("        Author:             Helton da Silva Gaspar");
  print sprintf("        Further info:       type print more_orbital");
  print sprintf("        HSG_lib home:       http://helton.paginas.ufsc.br/");
  print sprintf("        Author's contact:   helton.s.gaspar@ufsc.br");
  print sprintf("\n");
  print sprintf("print more_orbital");
