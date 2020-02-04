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

A simple & working Makefile is included [Makefile](Makefile). 


## Running
```
./QSDLA
```



## Simulation Example
**Plasma Wave Driven by Point-like Charge in Large-Bubble Regime** 
![logo](https://github.com/tianhongg/Quasi-Static-DLA-PIC/blob/master/Resource/Point_Driver.png)



**Intense Laser Pulse Propagates in the Tenuous Plasma** 
![logo](https://github.com/tianhongg/Quasi-Static-DLA-PIC/blob/master/Resource/Laser_Driver.png)



## Developing




## Author
* **Tianhong Wang (Cornell University)**(tw474@cornell.edu) 


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details



#### Reference
[About some of the algorithms](https://aip.scitation.org/doi/abs/10.1063/1.4999629)

[About DLA_1](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.184801); [DLA_2](https://aip.scitation.org/doi/abs/10.1063/1.5036967)

[About our group](https://shvets.aep.cornell.edu)
