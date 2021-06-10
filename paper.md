---
title: liborbital: A gnuplot library to deal with keplerian orbits, and more
tags:
  - gnuplot
  - astronomy
  - dynamics
  - keplerian orbits
  - orbital elements
authors:
  - name: Helton da Silva Gaspar^[first author]
    orcid: 0000-0002-0032-0966
    affiliation: 1
  - name: Ernesto Vieira Neto^[co author]
    orcid: 0000-0002-0032-0966
    affiliation: 2
affiliations:
 - name: UFSC - Universidade Federal de Santa Catarina
   index: 1
 - name: UNESP - Universidade Estadual Paulista "Júlio de Mesquita Filho"
   index: 2
date: 09 june 2021
bibliography: references.bib

---

# Summary

Most of science developed on dynamic astronomy is performed via numerical simulations. 
Numerical simulations main output data are position and velocity as function of the evolution time. 
However, the main analysis are performed over the oscullating orbital elements [@orbitalelements]. 
Thus, converting the Cartesian coordinates into orbital elements, and vice-versa, is a common task for the dynamicists. 
Additional algorithms use to be necessary to convert such output data into orbital osculating elements, and the developers of N-Body simulation software use to provide them along with their main software [@mercury],[@swift],[@rebound]. 
Normally, the main output data are converted into orbital elements and written in a secondary file, from which further analyses are dore such as the plotting of graphics, for instance. 
Sometimes, having to generate such secondary file may be inconvenient. 
The library allows one to plot the orbital elements converted directly from the main output data within the gnuplot. 
Furthermore, other physical quantities can also be evaluated such as orbital period, the flight time, radial and transversal components of the velocity vector, flight path angle and etc. 

Sometimes data must be duplicated into the second file due to the difficulty of plotting data from two different files. 
For instance, let us assume that one wants to plot the quotient of the current speed and the periapsis speed as function of the time. 
In such case, the current velocity components $v_x$, $v_y$ and $v_z$ stored into the original output data file should be copied into a secondary file along with the converted osculating elements. 

The <code>liborbital</code> library was not only designed to assist dynamic astronomers with plotting data, but it has other conveniences as it is described within the next section.



# Statement of need

The <code>liborbital</code> library is just a gnuplot script that can be easily loaded by typing the following command at the gnuplot prompt:
~~~
load `liborbital.plt'
~~~
Thus, it provides a set of functions which makes gnuplot a powerful tool to deal with many regular, but not trivial, tasks of the the celestial mechanics universe. 
The library allows one to easily accomplish tasks like:
 - Evaluating the position in the orbit for a given orbital flight time lapse, and vice-versa.
 - Evaluating any orbital element given the parameter of mass $\mu$ and the cartesian coordinates of the relative position and velocity vectors, vice-versa
 - Plotting an arc of 3D conic orbit given the orbital elements ($a$,$e$,$I$,$\Omega$,$\omega$) and the initial ($f_A$) and final ($f_B$) anomalies, or even the flight times ($t_A$) and ($t_B$) since the epoch of periapsis passage.
   - *This feature makes the library very useful for teaching of the orbital mechanics*
 - Plotting easily any of the keplerian orbital elements directly from a data file with the already mentioned state vectors components, and vice-versa.

Tasks like evaluating orbital elements given the state vectors, and vice-versa, use to be performed by writing short algorithms that must be compiled and run since the conversion routines are frequently provided along with the numerical simulation package. 
But the <code>liborbital</code> pack make such kind of tasks trivial. 
Although gnuplot was designed with the aim of plot graphics, it also allows one to evaluate functions on prompt console by using the gnuplot command <code>print</code>. 

# Applied example with Apophis close encounter
The purpose of this section is to exemplify how the use of the <code>liborbital</code> library avoids the need of generating two additional files when analyzing the results of a numerical simulation: 
 - <code>T\_Sun.dat</code> default ASCII output file which contains the Cartesian coordinates of the position vector $x$, $y$ and $z$ and the velocity vector $v_x$, $v_y$ and  $v_z$ for each body for each output timestep. %
 - <code>K\_Sun.dat</code> additional file which contains the *heliocentric* Keplerian osculating elements converted from <code>T\_Sun.dat</code> via auxiliary routines provided along with the numerical simulation package.
 - <code>K\_Earth.dat</code> additional file similar to the <code>K\_Sun.dat</code>, but with **geocentric** Keplerian osculating elements instead, also obtained via routines provided along with the numerical simulation package.

This applied example presents the results of a numerical simulation of the predicted close passage of the NEA 99942 Apophis (2004 MN4) by the Earth-Moon vicinity [@apophis]. 
Such astronomical phenomenon must occur in April 13$^\mathrm{th}$ 2029 when Apophis' closest distance to the surface of our planet will be about 5 Earth radii. 

Accurate data about the heliocentric configuration for the Earth, Moon and Apophis at March 1$^\mathrm{st}$ 2029 were obtained from JPL Horizons Web-Interface [@horizons]. 
Such initial configuration is dynamically evolved through a numerical simulation and the default output is written into <code>T\_Sun.dat</code> every 15 minutes. 
Figure \autoref{fig:ex1a} presents the trajectories of the three bodies during the close approach April 13$^\mathrm{th}\,0^\mathrm{h}$, with the Sun fixed at the origin of the coordinate system $(0\,\mathrm{AU},0\,\mathrm{AU},0\,\mathrm{AU})$. 
By using the functions <code>rv2...</code> of the <code>liborbital</code> library, the mean osculating elements $\left(\langle a\rangle,\langle e\rangle,\langle I\rangle,\langle\Omega\rangle,\langle\omega\rangle,\langle f\rangle\right)$ were evaluated just before and just after the closest approach from the Cartesian coordinates. 
From such mean values, the <code>liborbital</code> functions <code>eo2...</code> were used to plot the prior and posterior mean orbits with continuous line in peach color, and their respective extrapolations with dashed lines after and before the closest passage, respectively. 
The trajectories plotted with the <code>liborbital</code> illustrate very well how the the geocentric gravity field deviates the asteroid from its original orbit. 
Performing such task of plotting the mean orbits without the library, would demand the use of the additional file <code>K\_Sun.dat</code>. 

Figure \autoref{fig:ex1b} presents the osculating semimajor axis, perihelion and aphelion for the three bodies. 
An heuristic ellipse to explaing the lay readers the meaning of each parameters was drawn beside the key by using the <code>liborbital</code> library functions. 
Triangles and circles stand for the osculating elements from <code>K\_Sun.dat</code>, but the lines in peach color were evaluated from the default output file <code>T\_Sun.dat</code> via <code>liborbital</code> library. 

Figure \autoref{fig:ex1c} is similar to the figure \autoref{fig:ex1b}, but with respect to the geocentric frame of reference instead. 
Analogously to the figure \autoref{fig:ex1b}, the circles stand for data from file <code>K\_Earth.dat</code> while the back curves in peach color were plotted directly from <code>T\_Sun.dat</code> with <code>liborbital</code> functions. 
By observing the scale, one infers that the variations are very small $\sim10^{-3}$ either in semimajor axis as in eccentricity. 
Even for such small variations the library functions are very accurate. 

![Simulated trajectories are plotted with circles. A fit of the prior and posterior orbits to the closest approach were performed with the <code>liborbital</code> library, and are represented in peach color. \label{fig:ex1a}](ex1b.pdf){widht=70%}
![Circles and triangles stand for the osculating elements from <code>K\_Sun.dat</code> while the peach color curves were evaluated from <code>T\_Sun.dat</code> by using the <code>liborbital</code> functions. \label{fig:ex1b}](ex1.pdf){widht=70%}
![Analogous to figure \autoref{fig:ex1b} but for the geocentric frame of reference, instead. Osculating parameters computed with the library functions are so accurate as those computed with traditional routines even for such small variations. \label{fig:ex1c}](ex1c.pdf){widht=81%}

# Applied examples on orbital mechanics teaching
Figure \autoref{fig:ex2} is a classic sketch of orbital mechanics teaching that illustrates correlations between the mechanic energy and the types of conic orbits:
 - Circle $e=0$;
 - Ellipse $0<e<1$;
 - Parabola $e=1$;
 - Hyperbole $e>1$;

![Classic sketch of the orbital mechanics that exemplifies the correlation between the mechanic energy and the type of the conic orbit. The library functions provided not only the facility of plotting the orbits but also aided to align the eccentricity labels. \label{fig:ex2}](ex2.pdf){widht=81%}

The variety of orbits are easily plotted by using the library functions. 
Each orbit is labelled with its respective eccentricity that is properly aligned with the tangent direction to the orbit at the anchor point. 
Such proper alignment were easily achieved with the library functions provided to evaluate the perifocal velocity. 

Figure \autoref{fig:ex3} is a sketch of the three dimensional perspective of the spatial orientation of the orbit. 
Such plot is yielded with the library functions that convert Keplerian elements into Cartesian coordinates. 
Those same functions aided setting the many other elements represented in the figure, such as the perifocal and normal axes, the ascending and descending nodes as well as the nodes line, the velocity vector at the descending node point. 

![Classic sketch that illustrates the three Euler angles $\Omega$, $I$ and $\omega$ of sequential rotations that sets the spatial orientation of the orbit. \label{fig:ex3}](ex3.pdf){widht=81%}


# Acknowledgements
Thanks to Dr. Pablo Andretta Jaskowiak due the advises about software distribution platforms and suggestion of publishing at this journal. 
Dr. Ernesto Vieira Neto who developed the most user friendly N-Body numeric integrator in C language, that one employed to simulate the close encounter with Apophis. 

# References
