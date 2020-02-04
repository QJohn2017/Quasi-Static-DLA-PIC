# Quasi-Static-DLA-PIC


**Q**ausi-**S**tatic **D**irect **L**aser **A**cceleration -- **P**article **i**n **C**ell Simulation.

**QS-DLA-PIC** is a **serial** quasi-static particle-in-cell simulation code, written in the forgotten language **Fortran**. Developed by Tianhong Wang in 2017, at Dr. Gennady Shvets' Group, Cornell University.

You are welcome to clone or download this repository. :smile: We would also like to recommend you to try our new code [WAND-PIC](https://github.com/tianhongg/WAND-PIC). Please also send an email to the author Tianhong Wang(tw474@cornell.edu). We'd like to keep track of user numbers and affiliations. 



## Features
* 2D Planar Geometry or Cylindrical-Symmetric Geometry
* Adjustive longitudinal step according to the speed of plasma trajectories.
* Advanced quasi-static equations.
* Sub-cycling of macro beam particles.



## Prerequisites

* A NetCDF Library
```
netcdf
```


## Compiling
Just Make It! 

A simple & working Makefile is included [Makefile](QS_DLA_2D/Makefile). 


## Running
```
./QSDLA
```



## Simulation Example
**Plasma Wave Driven by Point-like Charge in Large-Bubble Regime** 
(Left: electron density; right: wake potential. Charge moving from right to left)
![logo](https://github.com/tianhongg/Quasi-Static-DLA-PIC/blob/master/Resource/Point_Driver.png)



**Intense Laser Pulse Propagates in the Tenuous Plasma** 
(Charge moving from left to right)
![logo](https://github.com/tianhongg/Quasi-Static-DLA-PIC/blob/master/Resource/Laser_Driver.png)



## Developing
QS-DLA-PIC is a FORTRAN-based serial simulation tool developed during my research. I first developed this code with Dr. Vladimir N. Khudik to verify some new equations we derived and also to simulate plasma wave driven by the point-charge. Then this code is upgraded to simulate laser-driven and beam driven plasma wakefield acceleration.  Although  QS-DLA is a 'single-core' simulation code, it's still very fast on most desktops due to the new equations we used. And you don't need to wait in the queue! It serves as a great tool if you just want to do some small-size simulations, or do a quick check of simulation parameters before you run larger 3D simulation on the supercomputers. More-powerful WAND-PIC has been developed by me later, but this is still a irreplaceable tool in my research.



## Author
* **Tianhong Wang (Cornell University)**(tw474@cornell.edu) 


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details



#### Reference
[About some of the algorithms](https://aip.scitation.org/doi/abs/10.1063/1.4999629)

[About DLA_1](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.184801); [DLA_2](https://aip.scitation.org/doi/abs/10.1063/1.5036967)

[About our group](https://shvets.aep.cornell.edu)
