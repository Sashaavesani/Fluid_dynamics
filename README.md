# Fluid_dynamics
NUMERICAL SOLUTION OF BOUNDARY LAYER OF A BOUNDING SURFACE

 es5
 
 NUMERICAL SOLUTION OF BOUNDARY LAYER OF A BOUNDING SURFACE (2D)
 Solve the PDE problem:
 
 { ∂u/∂x + ∂v/∂y = 0;                      CONTINUITY EQUATION
 
 { ∂u^2/∂x + ∂uv/∂y = nu*∂u^2/∂y^2         MOMENTUM EQUATION
 
 { u(x,0) = v(x,0) = 0; u(x,inf) = U;      Boundary Conditions
 
 { u(0,y) = U; v(0,y) = 0;                 Initial Conditions
 
 
es6

BLASIUS EQUATION:
NUMERICAL SOLUTION OF BOUNDARY LAYER OF A BOUNDING SURFACE (2D)
During the calculation of Blasius equation, in momentum equation 
we obtain:

U/nu*g(x)*g'(x) = - f'"(eta)/(f(eta)*f"(eta))

and this is satisfied IFF both
are equal to a positive constant called alpha.

Let alpha be equal to 1/2, then:

eta = y * sqrt(U/(nu*xb));

f'" + 1/2 f*f" = 0 (from which we get the diffeq problem below through the 
sostitution f' = u) with its boundary conditions

%%
Solve the diffeq problem

{ fu' + 2u" = 0

{ f' - u = 0

{ f(0) = 0; f'(inf) = 1

{ u(0) = 0; u(inf) = 1

using discretisation with finite difference
