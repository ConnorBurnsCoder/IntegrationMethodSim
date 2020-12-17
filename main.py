import sys
import math
import time

def get_constants():
    # constants are in order:
    #      mass, x_0, v_0, k_s
    return 1.0, 1.0, 0.0, 1.0

def get_accel(mass, k_s, x):
    F = -k_s*x
    return F/mass

def get_exact_x(mass, k_s, x_0, v_0, H):
    # uses the general solution for a simple oscillator to return exact position
    w = math.sqrt(k_s/mass)
    return x_0*math.cos(w*H)+(v_0/w)*math.sin(w*H)

def get_error(x, exact_x):
    return abs(x-exact_x)

def get_expected_error(H, n, k):
    # returns the expected error for order k method
    return math.pow(float(H)/n, k)

def get_energy(mass, k_s, x, v):
    return .5*mass*v*v + .5*k_s*x*x

def print_percent(t, H, perc):
    perc_new = math.floor(t / H * 10)
    if (perc_new > perc):
        print(str(perc_new * 10) + '%')
        return perc_new
    return perc

def euler(mass, k_s, x_0, v_0, H, n, file=None):
    # does one step of euler's method integration estimation
    # total error should be d_t^k for k=1
    # error per step should be about d_t^2
    d_t = H/n
    x = x_0
    v = v_0
    t = 0.0
    perc = -1
    comp_i = 0
    if file is not None:
        f = open(file, 'a')
    for i in range(int(n)):
        t += d_t
        a = get_accel(mass, k_s, x)
        x = x + d_t * v
        v = v + d_t * a
        perc = print_percent(t, H, perc)
        if file is not None:
            if comp_i == 0:  # print every 100
                f.write(str(t) + ' ' + str(x) + '\n')
                comp_i = 100
            comp_i -= 1
    if file is not None:
        f.write(str(H) + ' ' + str(x) + '\n')
        f.close()
    return x, v

def leap_frog(mass, k_s, x_0, v_0, H, n, file=None):
    # does one step of leap frog integration estimation
    # total error should be d_t^k for k=2
    # error per step should be about d_t^4
    d_t = float(H)/n
    x = x_0
    v_5 = v_0 + .5*d_t*get_accel(mass, k_s, x)
    t=0
    perc = -1
    comp_i = 0
    if file is not None:
        f = open(file, 'a')
    for i in range(int(n)):
        t += d_t
        x = x + d_t * v_5
        a = get_accel(mass, k_s, x)
        v_5 = v_5 + d_t * a
        perc = print_percent(t, H, perc)
        if file is not None:
            if comp_i==0:  # print every 100
                f.write(str(t) + ' ' + str(x) + '\n')
                comp_i = 100
            comp_i -= 1
    if file is not None:
        f.write(str(H) + ' ' + str(x) + '\n')
        f.close()
    return x, v_5

def rk4(mass, k_s, x_0, v_0, H, n, file=None):
    # does one step of Runge Kutta 4 integration estimation
    # total error should be d_t^k for k=4
    # error per step should be about d_t^8
    d_t = H / n
    x = x_0
    v = v_0
    t = 0
    perc = -1
    comp_i = 0
    if file is not None:
        f = open(file, 'a')
    for i in range(int(n)):
        t += d_t
        v1 = v
        a1 = get_accel(mass, k_s, x)

        v2 = v+a1*d_t/2.0
        a2 = get_accel(mass, k_s, x+v1*d_t/2.0)

        v3 = v + a2 * d_t / 2.0
        a3 = get_accel(mass, k_s, x + v2*d_t/2.0)

        v4 = v + a3 * d_t
        a4 = get_accel(mass, k_s, x + v3*d_t)

        x = x + (d_t/6.0)*(v1 + 2.0*v2 + 2.0*v3 + v4)
        v = v + (d_t / 6.0) * (a1 + 2*a2 + 2.0*a3 + a4)
        perc = print_percent(t, H, perc)
        if file is not None:
            if comp_i == 0:  # print every 100
                f.write(str(t) + ' ' + str(x) + '\n')
                comp_i = 100
            comp_i -= 1
    if file is not None:
        f.write(str(H) + ' ' + str(x) + '\n')
        f.close()
    return x, v

def exact_method(mass, k_s, x_0, v_0, H, n, file=None):
    #does integration by steps by getting the exact value
    d_t = H / n
    t = 0.0
    perc = -1
    comp_i = 0
    if file is not None:
        f = open(file, 'a')
    for i in range(int(n)):
        t += d_t
        x = get_exact_x(mass, k_s, x_0, v_0, t)
        perc = print_percent(t, H, perc)
        if file is not None:
            if comp_i == 0:  # print every 100
                f.write(str(t) + ' ' + str(x) + '\n')
                comp_i = 100
            comp_i -= 1
    if file is not None:
        f.write(str(H) + ' ' + str(x) + '\n')
        f.close()
    return x, 0

def integrate(method, method_name, k, H):
    #performs a single integration and prints to console
    mass, x_0, v_0, k_s = get_constants()
    # variable
    n = 10
    x, v = method(mass, k_s, x_0, v_0, H, n)
    exact_x = get_exact_x(mass, k_s, x_0, v_0, H)
    expected_error = get_expected_error(H, n, k)
    print('%s:\nn = %f\nd_t = %.30f\neuler_x = %.30f\nexact_x = %.30f' % (method_name, n, H / n, x, exact_x))
    print('expected_error_scaler (d_t^k) = %.30f' % (expected_error))

def integrate_iter(method, file_name, H, max, start, step):
    # iterates integrating increasing the number of samples by 10s from 10 to 10000
    f = open(file_name, 'w')
    f.close()
    mass, x_0, v_0, k_s = get_constants()
    n = start
    while n <= max:
        x, v = method(mass, k_s, x_0, v_0, H, n)
        exact_x = get_exact_x(mass, k_s, x_0, v_0, H)
        end_error = get_error(x, exact_x)
        with open(file_name, 'a') as f:
            f.write(str(n) + ' ' + str(end_error) + '\n')
        n += step
    print('iterated integration written to file: ' + file_name)

def integrate_steps(method, file_name, H):
    #performs a single integration and writes it to file
    mass, x_0, v_0, k_s = get_constants()
    # variable
    n = 1e7
    f = open(file_name, 'w')
    f.close()
    start_time = time.time()
    x, v = method(mass, k_s, x_0, v_0, H, n, file=file_name)
    run_time = time.time() - start_time
    exact_x = get_exact_x(mass, k_s, x_0, v_0, H)
    print('n = %d\nd_t = %.10f\nx = %.20f\nexact_x = %.20f\nerror = %.20f\nruntime(seconds) %s' % (n, H / n, x, exact_x, get_error(x, exact_x), run_time))
    print("{:e}".format(get_error(x, exact_x)))

def integrate_energy(method, H):
    # performs a single integration and prints the theoretical energy, actual energy, and error energy
    mass, x_0, v_0, k_s = get_constants()
    d_t = .01
    n = int(H/d_t)
    initial_energy = get_energy(mass, k_s, x_0, v_0)
    x, v = method(mass, k_s, x_0, v_0, H, n)
    final_energy = get_energy(mass, k_s, x, v)
    energy_error = get_error(initial_energy, final_energy)
    print('%s\nH = %d\nn = %d\nd_t = %.10f\nx = %.20f\nv = %.20f' % (method, H, n, d_t, x, v))
    print('init_E = %.10f\nfinal_E = %.10f\nE_error = %.10f' % (initial_energy, final_energy, energy_error))

def print_invalid_param():
    print("Invalid Parameters, specify e, l, or r for one iteration of Euler's, Leapfrog, or Runge-Kutta 4.")
    print('ex: $python main.py e')
    print("Specify e i, l i, or r i for iterations by 10.")
    print('ex: $python main.py e i')
    print("Specify e t, l t, or r t for duration time=1,000. (c t for exact solution)")
    print('ex: $python main.py l t')
    print("Specify e e y, l e y, or r e y for energy conservation check of duration time=1ey where y is an integer.")
    print('ex: $python main.py r e 2')

if __name__ == '__main__':
    if len(sys.argv) == 2:
        if sys.argv[1] == 'e':
            integrate(euler, 'Euler', 1, 1.0)
        elif sys.argv[1] == 'l':
            integrate(leap_frog, 'Leap_Frog', 2, 1.0)
        elif sys.argv[1] == 'r':
            integrate(rk4, 'RK4', 4, 1.0)
        else:
            print_invalid_param()
            sys.exit(-1)
    elif len(sys.argv) >= 3:
        if sys.argv[2] == 'i':
            if sys.argv[1] == 'e':
                integrate_iter(euler, 'Euler.txt', 1.0, 10000, 10, 10)
            elif sys.argv[1] == 'l':
                integrate_iter(leap_frog, 'Leap_Frog.txt', 1.0, 10000, 10, 10)
            elif sys.argv[1] == 'r':
                integrate_iter(rk4, 'RK4.txt', 1.0, 1000, 5, 5)
            else:
                print_invalid_param()
                sys.exit(-1)
        elif sys.argv[2] == 't':
            if sys.argv[1] == 'e':
                integrate_steps(euler, 'Euler.txt', 1000.0)
            elif sys.argv[1] == 'l':
                integrate_steps(leap_frog, 'Leap_Frog.txt', 1000.0)
            elif sys.argv[1] == 'r':
                integrate_steps(rk4, 'RK4.txt', 1000.0)
            elif sys.argv[1] == 'c':
                integrate_steps(exact_method, 'Exact.txt', 1000.0)
            else:
                print_invalid_param()
                sys.exit(-1)
        elif sys.argv[2] == 'e' and len(sys.argv) == 4:
            H = math.pow(10, int(sys.argv[3]))
            if sys.argv[1] == 'e':
                integrate_energy(euler, H)
            elif sys.argv[1] == 'l':
                integrate_energy(leap_frog, H)
            elif sys.argv[1] == 'r':
                integrate_energy(rk4, H)
            else:
                print_invalid_param()
                sys.exit(-1)
        else:
            print_invalid_param()
            sys.exit(-1)
    else:
        print_invalid_param()
        sys.exit(-1)