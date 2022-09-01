# 2D FLIP Simulation
2D simulation of an incompressible and inviscid fluid, based on the _Fluid Implicit Particle_ (FLIP) method<sup><a href="#1">1</a></sup>. The code implements an Eulerian approach over a Marker-and-Cell (MAC) grid to solve the [incompressible Navier-Stokes equations](https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations#Incompressible_flow):

$$ \frac{\partial\vec{u}}{\partial t} + \vec{u}\cdot\nabla\vec{u} + \frac{1}{\rho}\nabla{p} = \vec{g} + \nu\nabla\cdot\nabla\vec{u} $$

$$ \nabla\cdot\vec{u} = 0 $$

To solve these equations, the code uses an Eulerian approach over a 2D MAC (Marker and Cell) grid. The intermediate values of the physical magnitudes are obtained by either bilinear or Catmull-Rom interpolation.

### Level Set
In order to keep track of the air-liquid surface, the code implements a [level set](https://en.wikipedia.org/wiki/Level_set) representation of the interface. The interface is implicitly described by means of a function $\psi(\vec{r})$. Its value is defined in the cell centers along the grid as the signed distance to the interface’s closest point. Thus, the interface is defined as the set of points where $\psi(\vec{r}) = 0$. The time evolution of the surface results from the simple advection of the function $\psi$ in the velocity field.

$$ \displaystyle \frac{\partial\psi}{\partial t} = -\vec{u}\cdot\nabla\psi $$

### Requirements
```
qmake, libQt
```

### Compilation
```
$ qmake .
$ make
```

### References
<span id="1">[1]</a> R. Bridson, _Fluid Simulation for Computer Graphics_, A.K.Peters (2008).
