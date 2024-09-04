import sympy as sp
from sympy import *
import multiprocessing
import numpy as np
from functools import reduce
from sympy import cos, sin, sqrt, exp, pi

# Function to calculate the least common multiple of two numbers
def lcm(a, b):
    return abs(a * b) // sp.gcd(a, b)

# Function to calculate the LCM of multiple numbers
def lcm_multiple(numbers):
    return reduce(lcm, numbers)

accuracy = 60

# Initialize variables
sigma1 = sp.Rational('75')
sigma2 = sp.Rational('5')
sigma3 = sp.Rational('15')
sigma4 = sp.Rational('10')
sigma5 = sp.Rational('6')
sigma6 = sp.Rational('450')

# Convert sigmas to integers for LCM calculation
sigma1_int = int(sigma1)
sigma2_int = int(sigma2)
sigma3_int = int(sigma3)
sigma4_int = int(sigma4)
sigma5_int = int(sigma5)
sigma6_int = int(sigma6)

# Calculate rho values for each sigma
rho_1 = sp.N(1/(2 * sigma1) * 9, accuracy)
rho_2 = sp.N(1/(2 * sigma2) * 9, accuracy)
rho_3 = sp.N(1/(2 * sigma3) * 9, accuracy)
rho_4 = sp.N(1/(2 * sigma4) * 9 * 3, accuracy)
rho_5 = sp.N(1/(2 * sigma5) * 9, accuracy)
rho_6 = sp.N(1/(2 * sigma6) * 9, accuracy)

# Calculate total rho
rho_total = sp.N(rho_1 + rho_2 + rho_3 + rho_4 + rho_5 + rho_6, accuracy)
print("zcdp rho:", rho_total)

# Specify delta
delta = sp.Rational('1E-10')
print("zcdp Delta:", delta)

# Calculating epsilon
eps = sp.N(rho_total + 2 * sqrt(-rho_total * sp.log(delta)), accuracy)
print("zcdp Espilon", eps)

# Perform calculations
L = lcm_multiple([sigma1_int, sigma2_int, sigma3_int, sigma4_int, sigma5_int, sigma6_int])

teps_first = sp.N(eps - 4.5 * (1/sigma1 + 1/sigma2 + 1/sigma3 + 3/sigma4 + 1/sigma5 + 1/sigma6), accuracy)
teps_second = sp.N(eps + 4.5 * (1/sigma1 + 1/sigma2 + 1/sigma3 + 3/sigma4 + 1/sigma5 + 1/sigma6), accuracy)

N = sp.ceiling(teps_first * L)
N_2 = sp.ceiling(teps_second * L)
print("N:", N)
print("N_2:", N_2)

def char_func(sigma, hatsigma, t):
    u_range = range(1, int(sp.ceiling(20 * sqrt(9 * sigma))) + 1)
    sum_exp_cos = sum(sp.N(exp(-u**2 / (2 * 9 * sigma)) * cos(hatsigma * t * u), accuracy) for u in u_range)
    factor = sp.N(1 / sqrt(2 * pi * 9 * sigma), accuracy)
    return (1 + 2 * sum_exp_cos) * factor

def char_func_4(t):
    u_range = range(1, int(sp.ceiling(20 * sqrt(9 * 3 * sigma4))) + 1)
    sum_exp_cos = sum(sp.N(exp(-u**2 / (2 * 9 * 3 * sigma4)) * cos(L/sigma4 * t * u), accuracy) for u in u_range)
    factor = sp.N(1 / sqrt(2 * pi * 9 * 3 * sigma4), accuracy)
    return (1 + 2 * sum_exp_cos) * factor

# Precomputed weighted variances
sigmas = [sigma1, sigma2, sigma3, sigma5, sigma6]

# Precomputed weighted variances
hatsigmas = [L / sigma for sigma in [sigma1, sigma2, sigma3, sigma5, sigma6]]

# Define CHAR
def CHAR(t):
    product = 1
    for sigma in sigmas:
        product *= char_func(sigma, L / sigma, t)
    product *= char_func_4(t)
    return product

# Define weight functions
def weight_first(t):
    if t != 0:
        return (cos((3 * N_2 + 1) * t/2) * sin(3 * N_2 * t/2) / sin(t/2) - cos(N * t/2) * sin((N-1) * t/2) / sin(t/2))
    else:
        return 3 * N_2 - N + 1

def weight_second(t):
    if t != 0:
        return (cos((3 * N_2 + 1) * t/2) * sin(3 * N_2 * t/2) / sin(t/2) - cos(N_2 * t/2) * sin((N_2-1) * t/2) / sin(t/2))
    else:
        return 2 * N_2 + 1

# Define the expression to be integrated
def tobeint_first(t):
    return sp.N(1/(pi) * CHAR(t) * weight_first(t), accuracy)

# Define the expression to be integrated
def tobeint_second(t):
    return sp.N(1/(pi) * CHAR(t) * weight_second(t), accuracy)


# Number of intervals - change this to 10**16 if you must, but it's not recommended
n = 10**(2) + 1

h = (1/20)/(n - 1)
print("h:", h)

print(multiprocessing.cpu_count())

def boole(x, h, tobeint_function):
    return sp.N(2/45 * h * (7 * tobeint_function(x) + 32 * tobeint_function(x + h) + 12 * tobeint_function(x + 2*h) + 32 * tobeint_function(x + 3*h) + 7 * tobeint_function(x + 4*h)), accuracy)

# Define the worker function for parallel computation
def compute_boole_sum_segment(start, end, h, tobeint_function):
    total = 0
    for x in np.linspace(start, end, int((end - start) / (4 * h)) + 1)[:-1]:
        total += boole(x, h, tobeint_function)
    return total

def parallel_boole_sum(tobeint_function):
    num_segments = multiprocessing.cpu_count()
    segment_length = (1/20) / num_segments
    segments = [(i * segment_length, (i + 1) * segment_length, h, tobeint_function) for i in range(num_segments)]
    num_processes = multiprocessing.cpu_count()

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.starmap(compute_boole_sum_segment, segments)

    return sum(results)

# Perform the Boole sum
prob_first = sp.N(parallel_boole_sum(tobeint_first), accuracy)

# Print the result with high precision
print("New Delta First Term:", sp.N(prob_first, accuracy))

# Perform the Boole sum
prob_second = sp.N(parallel_boole_sum(tobeint_second), accuracy)

# Print the result with high precision
print("New Delta Second Term:", sp.N(prob_second, accuracy))
