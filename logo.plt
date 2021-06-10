load 'lib_orbital.plt'

dpi=2.0*pi;
hpi=0.5*pi;

mu=2.0;
a0=1.50;            e0=.33;             I0=pi/6; O0=hpi;  w0=hpi;
a1=polr(a0,e0,hpi); e1=0.0;             I1=0.0;  O1=0.0;  w1=0.0;
a2=3.50;            e2=0.0;             I2=0.0;  O2=0.0;  w2=0.0;
a3=(a1+a2)/2.0;     e3=(a2-a1)/(a2+a1); I3=0.0;  O3=0.0;  w3=0.0;

clr4='yellow'
clr2='orange'
clr3='spring-green'
clr1='blue'


set label 1 "A" at eo2x( mu, a0,e0,I0,O0,w0,hpi), eo2y( mu, a0,e0,I0,O0,w0,hpi), eo2z( mu, a0,e0,I0,O0,w0,hpi) point pt 7 ps 2 lc rgb 'red' tc rgb 'red' front offset -0.5 right
set label 2 "B" at eo2x( mu, a1,e1,I1,O1,w1,dpi), eo2y( mu, a1,e1,I1,O1,w1,dpi), eo2z( mu, a1,e1,I1,O1,w1,dpi) point pt 7 ps 2 lc rgb 'red' tc rgb 'red' front offset 0,-1
set label 3 "C" at eo2x( mu, a3,e3,I3,O3,w3, pi), eo2y( mu, a3,e3,I3,O3,w3, pi), eo2z( mu, a3,e3,I3,O3,w3, pi) point pt 7 ps 2 lc rgb 'red' tc rgb 'red' front offset -0.8,0.8

set label 4 "{/*1.8 HSGlib}"  at screen 0.95,0.87 tc rgb "orange-red" right front
set label 5 "{/*4.5 Orbital}" at screen 0.95,0.08 tc rgb "orange-red" right front

set arrow 10 from x0=eo2x( mu, a0,e0,I0,O0,w0,hpi), y0=eo2y( mu, a0,e0,I0,O0,w0,hpi), z0=eo2z( mu, a0,e0,I0,O0,w0,hpi)\
               to x0+eo2vx(mu, a0,e0,I0,O0,w0,hpi), y0+eo2vy(mu, a0,e0,I0,O0,w0,hpi), z0+eo2vz(mu, a0,e0,I0,O0,w0,hpi)\
               lc rgb clr1 lw 2.5

set arrow 20 from x1=eo2x( mu, a1,e1,I1,O1,w1,1.5*pi), y1=eo2y( mu, a1,e1,I1,O1,w1,1.5*pi), z1=eo2z( mu, a1,e1,I1,O1,w1,1.5*pi)\
               to x1+eo2vx(mu, a1,e1,I1,O1,w1,1.5*pi), y1+eo2vy(mu, a1,e1,I1,O1,w1,1.5*pi), z1+eo2vz(mu, a1,e1,I1,O1,w1,1.5*pi)\
               lc rgb clr2 lw 2.5 front

set arrow 12 from x0+eo2vx(mu, a0,e0,I0,O0,w0,hpi),    y0+eo2vy(mu, a0,e0,I0,O0,w0,hpi),    z0+eo2vz(mu, a0,e0,I0,O0,w0,hpi)\
               to x1+eo2vx(mu, a1,e1,I1,O1,w1,1.5*pi), y1+eo2vy(mu, a1,e1,I1,O1,w1,1.5*pi), z1+eo2vz(mu, a1,e1,I1,O1,w1,1.5*pi)\
               lc rgb 'red' lw 2.5 front

set arrow 40 from x0=eo2x( mu, a1,e1,I1,O1,w1,dpi), y0=eo2y( mu, a1,e1,I1,O1,w1,dpi), z0=eo2z( mu, a1,e1,I1,O1,w1,dpi)\
               to x0+eo2vx(mu, a1,e1,I1,O1,w1,dpi), y0+eo2vy(mu, a1,e1,I1,O1,w1,dpi), z0+eo2vz(mu, a1,e1,I1,O1,w1,dpi)\
               lc rgb clr2 lw 2.5

set arrow 30 from x1=eo2x( mu, a3,e3,I3,O3,w3,0.0), y1=eo2y( mu, a3,e3,I3,O3,w3,0.0), z1=eo2z( mu, a3,e3,I3,O3,w3,0.0)\
               to x1+eo2vx(mu, a3,e3,I3,O3,w3,0.0), y1+eo2vy(mu, a3,e3,I3,O3,w3,0.0), z1+eo2vz(mu, a3,e3,I3,O3,w3,0.0)\
               lc rgb clr3 lw 4.5
               
set arrow 34 from x0+eo2vx(mu, a1,e1,I1,O1,w1,dpi), y0+eo2vy(mu, a1,e1,I1,O1,w1,dpi), z0+eo2vz(mu, a1,e1,I1,O1,w1,dpi)\
               to x1+eo2vx(mu, a3,e3,I3,O3,w3,0.0), y1+eo2vy(mu, a3,e3,I3,O3,w3,0.0), z1+eo2vz(mu, a3,e3,I3,O3,w3,0.0)\
               lc rgb 'red' lw 2.5
               
set arrow 60 from x1=eo2x( mu, a3,e3,I3,O3,w3, pi), y1=eo2y( mu, a3,e3,I3,O3,w3, pi), z1=eo2z( mu, a3,e3,I3,O3,w3, pi)\
               to x1+eo2vx(mu, a3,e3,I3,O3,w3, pi), y1+eo2vy(mu, a3,e3,I3,O3,w3, pi), z1+eo2vz(mu, a3,e3,I3,O3,w3, pi)\
               lc rgb clr3 lw 4.5

set arrow 50 from x0=eo2x( mu, a2,e2,I2,O2,w2, pi), y0=eo2y( mu, a2,e2,I2,O2,w2, pi), z0=eo2z( mu, a2,e2,I2,O2,w2, pi)\
               to x0+eo2vx(mu, a2,e2,I2,O2,w2, pi), y0+eo2vy(mu, a2,e2,I2,O2,w2, pi), z0+eo2vz(mu, a2,e2,I2,O2,w2, pi)\
               lc rgb clr4 lw 4.5
               
set arrow 56 from x0+eo2vx(mu, a3,e3,I3,O3,w3, pi), y0+eo2vy(mu, a3,e3,I3,O3,w3, pi), z0+eo2vz(mu, a3,e3,I3,O3,w3, pi)\
               to x1+eo2vx(mu, a2,e2,I2,O2,w2, pi), y1+eo2vy(mu, a2,e2,I2,O2,w2, pi), z1+eo2vz(mu, a2,e2,I2,O2,w2, pi)\
               lc rgb 'red' lw 2.5
               
set arrow 70 from x0=eo2x( mu, a2,e2,I2,O2,w2,1.75*pi), y0=eo2y( mu, a2,e2,I2,O2,w2,1.75*pi), z0=eo2z( mu, a2,e2,I2,O2,w2,1.75*pi)\
               to x0+eo2vx(mu, a2,e2,I2,O2,w2,1.75*pi), y0+eo2vy(mu, a2,e2,I2,O2,w2,1.75*pi), z0+eo2vz(mu, a2,e2,I2,O2,w2,1.75*pi)\
               lc rgb clr4 lw 4.5 front
               

set parametric
set trange [0:1]
set urange [0:1]
set vrange [0:1]

set isosample 13,25

set xyplane at 0
set origin 0.18,0.00
set view 60,51.8,12.0,0.78

# Detele the following backslash to set png output\
set term png size 600,600 background rgb 'grey10' dashlength 1 font "Times,30" interlace;\
set output 'logo.png';\

# Detele the following backslash to set png output
set term jpeg size 400,400 background rgb 'grey20' dashlength 1 font "Times,20" interlace;\
set output 'logo.jpg';\


unset border
unset xtics
unset ytics
unset ztics
               
splot [][][-15:15][-15:15][-10:10]\
    lin(1.0,15,u)*cos(liborb_dpi*v),\
    lin(1.0,15,u)*sin(liborb_dpi*v),\
    0.0\
      notitle w l lt 0 lw 1 lc rgb 'olive'\
,\
    eo2x(mu,a0,e0,I0,O0,w0,dpi*u),\
    eo2y(mu,a0,e0,I0,O0,w0,dpi*u),\
    eo2z(mu,a0,e0,I0,O0,w0,dpi*u)\
      notitle 'insertionLEO' w l lw 1 lc rgb 'grey50' dt 0\
,\
    eo2x(mu,a0,e0,I0,O0,w0,hpi*u),\
    eo2y(mu,a0,e0,I0,O0,w0,hpi*u),\
    eo2z(mu,a0,e0,I0,O0,w0,hpi*u)\
      notitle 'transferorbit' w l lw 4.5 lc rgb clr1 \
,\
    sin(pi*u)*cos(liborb_dpi*v),\
    sin(pi*u)*sin(liborb_dpi*v),\
    cos(pi*u)\
      notitle w l lw 1 lc rgb 'steelblue'\
,\
    eo2x(mu,a1,e1,I1,O1,w1,dpi*u),\
    eo2y(mu,a1,e1,I1,O1,w1,dpi*u),\
    eo2z(mu,a1,e1,I1,O1,w1,dpi*u)\
      notitle 'equatorialLEO' w l lw 1.0 lc rgb 'grey50' lt 0 \
,\
    eo2x(mu,a2,e2,I2,O2,w2,dpi*u),\
    eo2y(mu,a2,e2,I2,O2,w2,dpi*u),\
    eo2z(mu,a2,e2,I2,O2,w2,dpi*u)\
      notitle 'targetorbit'   w l lw 1.0 lc rgb 'grey50' lt 0\
,\
    eo2x(mu,a3,e3,I3,O3,w3,dpi*u),\
    eo2y(mu,a3,e3,I3,O3,w3,dpi*u),\
    eo2z(mu,a3,e3,I3,O3,w3,dpi*u)\
      notitle 'transferorbit' w l lw 1.0 lc rgb 'grey50' lt 0 \
,\
    eo2x(mu,a1,e1,I1,O1,w1,lin(1.5,2.0,u)*pi),\
    eo2y(mu,a1,e1,I1,O1,w1,lin(1.5,2.0,u)*pi),\
    eo2z(mu,a1,e1,I1,O1,w1,lin(1.5,2.0,u)*pi)\
      notitle 'transferorbit' w l lw 4.5 lc rgb clr2 \
,\
    eo2x(mu,a3,e3,I3,O3,w3,pi*u),\
    eo2y(mu,a3,e3,I3,O3,w3,pi*u),\
    eo2z(mu,a3,e3,I3,O3,w3,pi*u)\
      notitle 'transferorbit' w l lw 4.5 lc rgb clr3 \
,\
    eo2x(mu,a2,e2,I2,O2,w2,lin(1.0,1.75,u)*pi),\
    eo2y(mu,a2,e2,I2,O2,w2,lin(1.0,1.75,u)*pi),\
    eo2z(mu,a2,e2,I2,O2,w2,lin(1.0,1.75,u)*pi)\
      notitle 'transferorbit' w l lw 4.5 lc rgb clr4 \
,\
    eo2x(mu,a2,e2,I2,O2,w2,1.75*pi)+0.2*sin(pi*u)*cos(dpi*v),\
    eo2y(mu,a2,e2,I2,O2,w2,1.75*pi)+0.2*sin(pi*u)*sin(dpi*v),\
    eo2z(mu,a2,e2,I2,O2,w2,1.75*pi)+0.2*cos(pi*u)\
      notitle 'transferorbit' w l lw 0.01 lc rgb 'orange-red';\
set output
