# 2D FLIP Simulation
2D simulation of an incompressible and inviscid fluid, based on the _Fluid Implicit Particle_ (FLIP) method<sup><a href="#1">1</a></sup>. The code implements an Eulerian approach over a Marker-and-Cell (MAC) grid to solve the [incompressible Navier-Stokes equations](https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations#Incompressible_flow):

![](https://latex.codecogs.com/svg.image?\frac{\partial&space;\vec{u}}{\partial&space;t}&space;&plus;&space;\vec{u}\cdot\nabla\vec{u}&space;&plus;&space;\frac{1}{\rho}\nabla{p}&space;=&space;\vec{g}&space;&plus;&space;\nu\nabla\cdot\nabla\vec{u})

![](https://latex.codecogs.com/svg.image?\nabla\cdot\vec{u}&space;=&space;0)

In order to keep track of the air-liquid surface, the method uses a [level set](https://en.wikipedia.org/wiki/Level_set) representation of the interface.

### Compilation
```
$ qmake .
$ make
```

### References
<span id="1">[1]</a> R. Bridson, _Fluid Simulation for Computer Graphics_, A.K.Peters (2008).
