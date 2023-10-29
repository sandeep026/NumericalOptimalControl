# Problem description

The optimal control problem for the cart system [^1] is provided below. z~1 and z~2 are the position and velocity of the cart and they comprise the states. f is the force applied and there is a drag force which is proportional to the velocity of the cart. The system starts from rest and additional boundary condition is placed at the end of the trajectory. Along the trajectory the control effort is minimized from time 0 to 2.

![image](https://user-images.githubusercontent.com/16457676/236567436-9d87b891-e74f-4299-802c-a394693c1f60.png)

# Analytical solution

The system admits the following analytical solution, which can be later used to verify the numerical solution and its accuracy.

![image](https://user-images.githubusercontent.com/16457676/236629178-b6da4837-b1d8-454d-9ec4-2d67fb1abeba.png)

# Numerical optimal control
Direct methods [^2] are used to solve optimal control problem for a cart system. The optimal control problem is discretized and transformed to a nonlinear optimization problem and solved using the state of the art solver IPOPT. The numerical methods used for discretization are

1. single shooting with piecewise constant control (RK45)
2. multiple shooting with piecewise constant control (RK45)
3. trapezoidal with piecewise linear/constant control 
4. hermite simpson with piecewise linear control
5. LGL pseudospectral

# Requirements
- MATLAB/[OCTAVE](https://octave.org/)
- [Casadi](https://web.casadi.org/)

# References

[^1]: Conway, B. A. and K. Larson (1998). Collocation versus differential inclusion in direct optimization. Journal of Guidance, Control, and Dynamics, 21(5), 780–785
[^2]: Diehl, Moritz, and Sébastien Gros. "Numerical optimal control." Optimization in Engineering Center (OPTEC) (2011).
