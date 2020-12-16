# IntegrationMethodSim
Test the order of integration methods by finding the error(Euler's method (k=1), Leapfrog (k=2), and Runge-Kutta 4 (k=4)). Written in Python


To demonstrate this we will use a system where an object with mass = 1 is attached to a spring with stiffness ks = 1. The demonstration will start at t=0 where the position is x0 = 1. The initial velocity (v0) is 0. Each method will estimate the position of the object until the end of the simulation (H=end of simulation time). Meaning they will return the position of the object at xH. 

The force is calculated with Hooke’s law where F = -ks*x and the acceleration is calculated with Newton’s second law F = ma. To view the order we must compare the error between the estimation and an exact calculation. The exact calculation is found using the general solution for a simple oscillator: x(H) = x0 cos(sqrt(ks/m)H) + (v0/sqrt(ks/m))*sin(sqrt(ks/m)H). 
Solved for example H=1:
x(1) = cos(1) = .5403...


Specify e, l, or r for one iteration of total time H=1 of Euler's, Leapfrog, or Runge-Kutta 4.

ex: $python main.py e



Specify e i, l i, or r i for iterations by of desired method

ex: $python main.py e i



Specify e t, l t, r t, or c t for duration time H=1,000. (c t for exact solution). Data will be written to Euler.txt, Leap_Frog.txt, RK4.txt, or Exact.txt

ex: $python main.py l t



Specify e e y, l e y, or r e y for energy conservation check of duration time=10^y where y is an integer.

ex: $python main.py r e 2
