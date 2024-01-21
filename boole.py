import mpmath
from mpmath import mp
import multiprocessing
import math
import numpy as np
from functools import reduce

# Function to calculate the least common multiple of two numbers
def lcm(a, b):
    return abs(a * b) // math.gcd(a, b)

# Function to calculate the LCM of multiple numbers
def lcm_multiple(numbers):
    return reduce(lcm, numbers)

# Set precision for mpmath
mp.dps = 70  # Set the decimal places

# Initialize variables using mpf
sigma1 = mp.mpf('75')
sigma2 = mp.mpf('5')
sigma3 = mp.mpf('15')
sigma4 = mp.mpf('10')
sigma5 = mp.mpf('6')
sigma6 = mp.mpf('450')

# Convert sigmas to integers for LCM calculation
sigma1_int = int(sigma1)
sigma2_int = int(sigma2)
sigma3_int = int(sigma3)
sigma4_int = int(sigma4)
sigma5_int = int(sigma5)
sigma6_int = int(sigma6)

# Calculate rho values for each sigma
rho_1 = mp.mpf('1')/(mp.mpf('2') * sigma1) * mp.mpf('9')
rho_2 = mp.mpf('1')/(mp.mpf('2') * sigma2) * mp.mpf('9')
rho_3 = mp.mpf('1')/(mp.mpf('2') * sigma3) * mp.mpf('9')
rho_4 = mp.mpf('1')/(mp.mpf('2') * sigma4) * mp.mpf('9') * mp.mpf('3')
rho_5 = mp.mpf('1')/(mp.mpf('2') * sigma5) * mp.mpf('9')
rho_6 = mp.mpf('1')/(mp.mpf('2') * sigma6) * mp.mpf('9')

# Calculate total rho
rho_total = rho_1 + rho_2 + rho_3 + rho_4 + rho_5 + rho_6
print("Rho_zcdp:", mp.nstr(rho_total, n=70))

# Specify delta
delta = mp.mpf('1E-10')
print("Delta_zcdp:", mp.nstr(delta, n=70))

# Calculating epsilon
eps = rho_total + mp.mpf('2') * mp.sqrt(-rho_total * mp.log(delta))
print("eps_zcdp:", mp.nstr(eps, n=70))

# Perform calculations
L = lcm_multiple([sigma1_int, sigma2_int, sigma3_int, sigma4_int, sigma5_int, sigma6_int])
print("L:", L)

teps_first = eps - mp.mpf('4.5') * (1/sigma1 + 1/sigma2 + 1/sigma3 + 3/sigma4 + 1/sigma5 + 1/sigma6)
teps_second = eps + mp.mpf('4.5') * (1/sigma1 + 1/sigma2 + 1/sigma3 + 3/sigma4 + 1/sigma5 + 1/sigma6)

N = mp.ceil(teps_first * L)
N_2 = mp.ceil(teps_second * L)
print("N:", N)
print("N_2:", N_2)

# Define the characteristic functions using mpf
def char_func(sigma, hatsigma, t):
    u_range = range(1, int(mp.ceil(20 * mp.sqrt(9 * sigma))) + 1)
    sum_exp_cos = sum(mp.exp(-u**2 / (2 * 9 * sigma)) * mp.cos(hatsigma * t * u) for u in u_range)
    factor = 1 / mp.sqrt(2 * mp.pi * 9 * sigma)
    return (1 + 2 * sum_exp_cos) * factor

def char_func_4(t):
    u_range = range(1, int(mp.ceil(20 * mp.sqrt(9 * 3 * sigma4))) + 1)
    sum_exp_cos = sum(mp.exp(-u**2 / (2 * 9 * 3 * sigma4)) * mp.cos(L/sigma4 * t * u) for u in u_range)
    factor = 1 / mp.sqrt(2 * mp.pi * 9 * 3 * sigma4)
    return (1 + 2 * sum_exp_cos) * factor

# Precomputed weighted variances
sigmas = [sigma1, sigma2, sigma3, sigma5, sigma6]

# Precomputed weighted variances
hatsigmas = [L / sigma for sigma in [sigma1, sigma2, sigma3, sigma5, sigma6]]

# Define CHAR using mpf
def CHAR(t):
    t = mp.mpf(t)
    product = mp.mpf('1')
    for sigma in sigmas:
        product *= char_func(sigma, L / sigma, t)
    product *= char_func_4(t)
    return product

# Define weight function
def weight_first(t):
  t = mp.mpf(t)
  if t != 0:
    return mp.cos((3 * N_2 + 1) * t/2) * mp.sin(3 * N_2 * t/2) / mp.sin(t/2) - mp.cos(N * t/2) * mp.sin((N-1) * t/2) / mp.sin(t/2)
  if t == 0:
    return 3 * N_2 - N + 1

# Define weight function
def weight_second(t):
  t = mp.mpf(t)
  if t != 0:
    return  mp.cos((3 * N_2 + 1) * t/2) * mp.sin(3 * N_2 * t/2) / mp.sin(t/2) - mp.cos(N_2 * t/2) * mp.sin((N_2-1) * t/2) / mp.sin(t/2)
  if t == 0:
    return 2 * N_2 + 1

# Define the expression to be integrated
def tobeint_first(t):
    return 1/(mp.pi) * CHAR(t) * weight_first(t)

# Define the expression to be integrated
def tobeint_second(t):
    return 1/(mp.pi) * CHAR(t) * weight_second(t)

# Number of intervals - change this to 10**16 if you must, but it's not recommended
n = 10**(7) + 1

h = (1/20)/(n - 1)
print("h:", h)

print(multiprocessing.cpu_count())

print("Maximal Approximate Error:", (10**5)**6 * h**6)

# Define the worker function for parallel computation
def compute_boole_sum_segment(start, end, h, tobeint_function):
    return sum([2/45 * h * ( 7 * tobeint_function(x) + 32 * tobeint_function(x + h) + 12 * tobeint_function(x + 2*h) + 32 * tobeint_function(x + 3*h) + 7 * tobeint_function(x + 4*h) )
      for x in mp.linspace(start, end, int((end - start) / (4 * h)) + 1)[:-1]])

# Main computation with a pool of workers
def parallel_boole_sum(tobeint_function):

    # Define how to split the range [0, mp.pi]
    num_segments = multiprocessing.cpu_count()  # Number of segments to split the task into
    segment_length = (1/20) / num_segments

    # Create tuples representing each segment
    segments = [(i * segment_length, (i + 1) * segment_length, h, tobeint_function) for i in range(num_segments)]

    # Set the desired number of parallel processes
    num_processes = multiprocessing.cpu_count()  # Change this number to your desired value

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.starmap(compute_boole_sum_segment, segments)

    return sum(results)


# Perform boole Sum
prob_first = parallel_boole_sum(tobeint_first)

# Print the result with high precision
print("New Delta First Term:", mp.nstr(prob_first, n=70))


# Perform boole Sum
prob_second = parallel_boole_sum(tobeint_second)

# Print the result with high precision
print("New Delta Second Term:", mp.nstr(prob_second, n=70))


