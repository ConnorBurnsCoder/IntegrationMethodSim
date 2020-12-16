# IntegrationMethodSim
Test the order of integration methods by finding the error(Euler's method (k=1), Leapfrog (k=2), and Runge-Kutta 4 (k=4)). Written in Python

Specify e, l, or r for one iteration of total time H=1 of Euler's, Leapfrog, or Runge-Kutta 4.
ex: $python main.py e

Specify e i, l i, or r i for iterations by of desired method
ex: $python main.py e i

Specify e t, l t, r t, or c t for duration time H=1,000. (c t for exact solution). Data will be written to Euler.txt, Leap_Frog.txt, RK4.txt, or Exact.txt
ex: $main.py l t

Specify e e y, l e y, or r e y for energy conservation check of duration time=10^y where y is an integer.
ex: $main.py r e 2
