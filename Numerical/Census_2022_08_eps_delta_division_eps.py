import mpmath as mp
import time
import sys
start_time = time.time()



############# Precision #############
accuracy = 70
mp.dps = accuracy



############# Loop over all k #############
if len(sys.argv) > 1:
    k = int(sys.argv[1])
else:
    k = 0 



############# Parameters #############
N_CPU = mp.mpf('96')
N_CPU_int = int(N_CPU)

rho_base = mp.mpf('3.65')
percentage = mp.mpf('0.9485')
eps_base = mp.mpf('21.97')

eps = percentage * eps_base
rho = rho_base

hundred = mp.mpf('100')
n = mp.mpf('10')
L = mp.mpf('1000')

zero = mp.mpf('0')
one = mp.mpf('1')
two = mp.mpf('2')
three = mp.mpf('3')
four = mp.mpf('4')
six = mp.mpf('6')
seven = mp.mpf('7')
ten = mp.mpf('10')
twelve = mp.mpf('12')
twenty_five = mp.mpf('25')
thirty_two = mp.mpf('32')
forty_five = mp.mpf('45')
two_over_45 = two / forty_five



############# Privacy Budget #############
aL_list = [
    mp.mpf('20'), mp.mpf('274'), mp.mpf('85'), mp.mpf('131'), 
    mp.mpf('238'), mp.mpf('118'), mp.mpf('3')
]

percentages = [
    mp.mpf('2'), mp.mpf('27.40'), mp.mpf('8.5'), mp.mpf('13.10'),
    mp.mpf('23.80'), mp.mpf('11.80'), mp.mpf('0.3')
]

rho_list = [rho * percentage / hundred for percentage in percentages]

n_list = [
    mp.mpf('10'), mp.mpf('10'), mp.mpf('10'), mp.mpf('20'),
    mp.mpf('10'), mp.mpf('10'), mp.mpf('10')
]

sigma2_list = []
for rhoo in rho_list:
    sigma2 = one / (two * rhoo / n)
    sigma2_list.append(sigma2)



############# Summation Parameters #############
t_eps = n / two * (eps / rho - one)
T_eps = n / two * (eps / rho + one)
k_comp = seven  # Number of independent compositions

N = mp.ceil(t_eps * L)
print('Summation Lower Bound in first term:', N)
N_2 = mp.ceil(T_eps * L)
print('Summation Lower Bound in second term:', N_2)
U = mp.ceil(six * t_eps * L)
print('Summation Upper Bound:', U)



############# Characteristic Functions #############
def char_func(sigma, aL, n, t):
    # Compute the maximum value of u
    u_max = int(mp.ceil(twenty_five * mp.sqrt(n * sigma)) + one)
    sum_exp_cos = zero
    for u in range(1, u_max):
        u_mpf = mp.mpf(u)
        term = mp.exp(-u_mpf ** 2 / (two * n * sigma)) * mp.cos(aL * t * u_mpf)
        sum_exp_cos += term
    factor = one / mp.sqrt(two * mp.pi * n * sigma)
    result = (one + two * sum_exp_cos) * factor
    return result

def char_func_prod(t):
    product = one
    for i in range(len(sigma2_list)):
        sigma = sigma2_list[i]
        aL = aL_list[i]
        n_i = n_list[i]
        product *= char_func(sigma, aL, n_i, t)
    return product



############# Weight #############
def weight_first(t):
    stable_cut = mp.mpf('1e-20')
    if abs(t) >= stable_cut:
        first = mp.cos(N * t) + mp.cos(U * t)
        factor = mp.cos(t / two) / mp.sin(t / two)
        second = mp.sin(U * t) - mp.sin(N * t)
        return (one / two) * (first + factor * second)
    
    elif 0 < abs(t) < stable_cut:
        first = mp.cos(N * t) + mp.cos(U * t)
        factor = mp.cos(t / two) / (t / two)
        second = mp.sin(U * t) - mp.sin(N * t)
        return (one / two) * (first + factor * second)
    
    else:
        return U - N + one

def weight_second(t):
    stable_cut = mp.mpf('1e-20')
    if abs(t) >= stable_cut:
        first = mp.cos(N_2 * t) + mp.cos(U * t)
        factor = mp.cos(t / two) / mp.sin(t / two)
        second = mp.sin(U * t) - mp.sin(N_2 * t)
        return (one / two) * (first + factor * second)
    
    elif 0 < abs(t) < stable_cut:
        first = mp.cos(N_2 * t) + mp.cos(U * t)
        factor = mp.cos(t / two) / (t / two)
        second = mp.sin(U * t) - mp.sin(N_2 * t)
        return (one / two) * (first + factor * second)
    
    else:
        return U - N_2 + one



############# Integrant #############
def tobeint_first(t):
    return (one / mp.pi) * char_func_prod(t) * weight_first(t)

def tobeint_second(t):
    return (one / mp.pi) * char_func_prod(t) * weight_second(t)



############# Boole Sum #############
# Number of intervals
n_intervals = int(ten ** 5)
print("Number of intervals:", n_intervals)

# Define the integration limits for each k
lower_limit = mp.mpf(k) * (one / (N_CPU * hundred))
upper_limit = mp.mpf(k + 1) * (one / (N_CPU * hundred))
int_limit = upper_limit - lower_limit
h = int_limit / n_intervals
print(f"Integration limits for k={k}: [{lower_limit}, {upper_limit}]")
print("h:", h)

def boole(x, h, tobeint_function):
    return two_over_45 * h * (
        seven * tobeint_function(x)
        + thirty_two * tobeint_function(x + h)
        + twelve * tobeint_function(x + two * h)
        + thirty_two * tobeint_function(x + three * h)
        + seven * tobeint_function(x + four * h)
    )

def boole_sum_first():
    total = zero
    num_steps = int(n_intervals / four)
    for i in range(num_steps):
        x = lower_limit + i * four * h
        total += boole(x, h, tobeint_first)
    return total

def boole_sum_second():
    total = zero
    num_steps = int(n_intervals / four)
    for i in range(num_steps):
        x = lower_limit + i * four * h
        total += boole(x, h, tobeint_second)
    return total



############# Print #############
prob_first = boole_sum_first()
print(f"New Delta First Probability for k={k}:", mp.nstr(prob_first, n=70))

# Save outputs to a file
with open(f'/home/ec2-user/output_{k}.txt', 'w') as f:
    f.write(f"epsilin: {mp.nstr(eps, n=70)}\n")
    f.write(f"rho: {mp.nstr(rho, n=70)}\n")

# Save outputs to a file
with open(f'/home/ec2-user/output_{k}.txt', 'a') as f:
    f.write(f"New Delta First Probability for k={k}: {mp.nstr(prob_first, n=70)}\n")
    f.write(f"Process finished --- {time.time() - start_time} seconds ---\n")

# # Perform the Boole sum for the second term
# prob_second = boole_sum_second()
# weighted_prob_second = prob_second * mp.exp(eps)
# print(f"New Delta Second Probability for k={k}:", mp.nstr(prob_second, n=70))
# print(f"New Delta Second Term for k={k}:", mp.nstr(weighted_prob_second, n=70))

# print("Process finished --- %s seconds ---" % (time.time() - start_time))

# # Save outputs to a file
# with open(f'/home/ec2-user/output_{k}.txt', 'a') as f:
#     f.write(f"New Delta Second Probability for k={k}: {mp.nstr(prob_second, n=70)}\n")
#     f.write(f"New Delta Second Term for k={k}: {mp.nstr(weighted_prob_second, n=70)}\n")

