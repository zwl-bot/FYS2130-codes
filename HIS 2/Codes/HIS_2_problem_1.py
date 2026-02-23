import numpy as np
import matplotlib.pyplot as plt


## initial terms
L = 2.5
g = 9.81
m = 30                          # mass (kg)
F0 = 15                         # driving force amplitude
k = 0.8                         # coupling strength
gamma = 0.15                    # damping coefficient
w0 = np.sqrt(g/L)               # natural frequency
wd = np.sqrt(w0**2 - gamma**2)  # damped natural frequancy
wr = np.sqrt(w0**2 - gamma**2/2)# resonance frequency

T = 60
dt = 0.001
N = int(T/dt)


## initial conditions
theta1_0 = 0.3; d_theta1_0 = 0 
theta2_0 = 0; d_theta2_0 = 0


## initiate arrays
t_s = np.linspace(0, T, N+1)
E1_s = np.zeros(N+1); E2_s = np.zeros(N+1)



## numerical scheme function
def numerical_Euler(t_s, theta1_0, d_theta1_0, theta2_0, d_theta2_0, 
                    w0, L, m, F0, k, gamma, w):
    N = len(t_s)-1; dt = t_s[1]
    theta1 = np.zeros(N+1); d_theta1 = np.zeros(N+1)
    theta2 = np.zeros(N+1); d_theta2 = np.zeros(N+1)
    
    theta1[0] = theta1_0; d_theta1[0] = d_theta1_0
    theta2[0] = theta2_0; d_theta2[0] = d_theta2_0
    
    for n in range(N):
        theta1[n+1] = theta1[n] + dt * d_theta1[n]
        theta2[n+1] = theta2[n] + dt * d_theta2[n]
        driv_F = F0/(m*L) * np.cos(w*t_s[n])
        d_theta1[n+1] = d_theta1[n] + dt*(driv_F - 2*gamma*d_theta1[n] - w0**2 * theta1[n]
                                          - k*(theta1[n] - theta2[n]))
        d_theta2[n+1] = d_theta2[n] - dt*(2*gamma*d_theta2[n] + w0**2 * theta2[n]
                                          + k*(theta2[n] - theta1[n]))
    
    return theta1, theta2, d_theta1, d_theta2


## find max resonance amplitude
def find_amplitude_resonance(theta_lst):
    start = int(0.8*len(theta_lst))     # Searching in the list starts from start
    theta_max = max(theta_lst[start:]); theta_min = min(theta_lst[start:])
    A = (theta_max - theta_min)/2
    return A










## Solve numerical scheme c)
## Calculate total mechanical energy for each system

########## plot c) ##########

ws = [0.5*w0, 1.5*w0, w0, -w0, 2*w0]              # different driving frequency values


theta1_ws = [0]*len(ws); E1_ws = [0]*len(ws)
theta2_ws = [0]*len(ws); E2_ws = [0]*len(ws)


for i in range(len(ws)):
    theta1, theta2, d_theta1, d_theta2 = numerical_Euler(t_s, theta1_0, d_theta1_0, theta2_0, d_theta2_0, w0, L, m, F0, k, gamma, ws[i])
    theta1_ws[i] = theta1; theta2_ws[i] = theta2
    E1_ws[i] = 0.5 * d_theta1**2 + 0.5*w0**2 * theta1**2
    E2_ws[i] = 0.5 * d_theta2**2 + 0.5*w0**2 * theta2**2

for i in range(len(ws)):
    plt.plot(t_s, theta1_ws[i], label="theta1")
    plt.plot(t_s, theta2_ws[i], label="theta2", linestyle="dashed")
    plt.title(f"w={ws[i]:.3f}, dt={dt}, w0={w0:.3f}, gamma={gamma}, k={k}")
    plt.xlabel("t [s]"); plt.ylabel("theta"); plt.legend(); plt.show()

for i in range(len(ws)):
    plt.plot(t_s, E1_ws[i], label="theta1")
    plt.plot(t_s, E2_ws[i], label="theta2", linestyle="dashed")
    #plt.yscale("log")
    plt.title(f"E(t) w={ws[i]:.3f}, dt={dt}, w0={w0:.3f}, gamma={gamma}, k={k}")
    plt.xlabel("t [s]"); plt.ylabel("E(t)"); plt.legend(); plt.show()





#####################################################################




########## d) ##########
w = wr
ks = [0, 0.2, 0.4, k]; gs = [0.05, 0.10, gamma, 0.20] 


"""

theta1_ks = [0]*len(ks); E1_ks = [0]*len(ks)
theta1_gs = [0]*len(gs); E1_gs = [0]*len(gs)

count = 0
for k_val in ks:
    curr_theta1, curr_theta2, curr_dtheta1, curr_dtheta2 = numerical_Euler(t_s, theta1_0, d_theta1_0, theta2_0, d_theta2_0, w0, L, m, F0, k_val, gamma, w)
    theta1_ks[count] = curr_theta1
    E1_ks[count] = 0.5*curr_dtheta1**2 + 0.5*w0**2 * curr_theta1**2
    count += 1

count = 0
for g_val in gs:
    curr_w = np.sqrt(w0**2 - gs[count]**2/2)
    curr_theta1, curr_theta2, curr_dtheta1, curr_dtheta2 = numerical_Euler(t_s, theta1_0, d_theta1_0, theta2_0, d_theta2_0, w0, L, m, F0, k, g_val, curr_w)
    theta1_gs[count] = curr_theta1
    E1_gs[count] = 0.5*curr_dtheta1**2 + 0.5*w0**2 * curr_theta1**2
    count += 1



for i in range(len(ks)):
    plt.plot(t_s, theta1_ks[i], label=f"k={ks[i]}")
plt.title(f"Theta1 w={w:.3f}, dt={dt}, w0={w0:.3f}, gamma={gamma}")
plt.xlabel("t [s]"); plt.ylabel("theta"); plt.legend(); plt.show()
    

for i in range(len(ks)):
    plt.plot(t_s, E1_ks[i], label=f"k={ks[i]}")
plt.title(f"Theta1 E(t) w={w:.3f}, dt={dt}, w0={w0:.3f}, gamma={gamma}")
plt.xlabel("t [s]"); plt.ylabel("E(t)"); plt.legend(); plt.show()

'''Amplitude is smaller when k is larger. Can show this relation with equation
k=0 is the no coupling case and just a single driven motion'''



for i in range(len(gs)):
    curr_w = np.sqrt(w0**2 - gs[i]**2/2)
    plt.plot(t_s, theta1_gs[i], label=f"gm={gs[i]}, w={curr_w:.4f}")
plt.title(f"Theta1 dt={dt}, w0={w0:.3f}, k={k}")
plt.xlabel("t [s]"); plt.ylabel("theta"); plt.legend(); plt.show()


for i in range(len(gs)):
    curr_w = np.sqrt(w0**2 - gs[i]**2/2)
    plt.plot(t_s, E1_gs[i], label=f"gm={gs[i]}, w={curr_w:.4f}")
plt.yscale("log")
plt.title(f"Theta1 E(t) dt={dt}, w0={w0:.3f}, k={k}")
plt.xlabel("t [s]"); plt.ylabel("E(t)"); plt.legend(); plt.show()

'''Maximal amplitude is reached when w=w0, i.e. when gm=0
gm=0 is the no damping case. Because of coupling, the energy will just increase over time.
Having damping will automatically prevent oscillation amplitude to become infinite 
over time'''

"""



#####################################################################






########## d) Find amplitude of resonance ##########

"""
dt = 0.0001
N_new = int(T/dt)

t_s_new = np.linspace(0, T, N_new+1)
theta1, theta2, d_theta1, d_theta2 = numerical_Euler(t_s_new, theta1_0, d_theta1_0, theta2_0, d_theta2_0, w0, L, m, F0, k, gamma, w)
A1 = find_amplitude_resonance(theta1); A2 = find_amplitude_resonance(theta2)

plt.plot(t_s_new, theta1, label="theta1")
plt.plot(t_s_new, theta2, label="theta2", linestyle="dashed")
plt.title(f"w={w:.3f}, dt={dt}, w0={w0:.3f}, gamma={gamma}, k={k}")
plt.xlabel("t [s]"); plt.ylabel("theta"); plt.legend(); plt.show()

print(f"w={w:.3f}, dt={dt}, w0={w0:.3f}, gamma={gamma}, k={k}")
print(f"Amplitude theta1: A1={A1:.3f}")
print(f"Amplitude theta2: A2={A2:.3f}")


new_theta1, new_theta2, new_d_theta1, new_d_theta2 = numerical_Euler(t_s_new, theta1_0, d_theta1_0, theta2_0, d_theta2_0, w0, L, m, F0, 0, gamma, w)
A1 = find_amplitude_resonance(new_theta1) # single driven oscillator
print(f"Amplitude theta1 (single driven): A={A1:.3f}")


'''Its worth noting that as dt gets smaller, its amplitude also converges
towards the theroetical A_max for single driven oscillators'''
"""





#####################################################################



########## e) beat frequency with no external driving ##########

ks = [0.4, k, 1.2]; gs = [0.01, 0.05, gamma]

theta1_ks = [0]*len(ks); E1_ks = [0]*len(ks)
theta1_gs = [0]*len(gs); E1_gs = [0]*len(gs)

theta2_ks = [0]*len(ks); E2_ks = [0]*len(ks)
theta2_gs = [0]*len(gs); E2_gs = [0]*len(gs)


count = 0
for k_val in ks:
    curr_theta1, curr_theta2, curr_dtheta1, curr_dtheta2 = numerical_Euler(t_s, theta1_0, d_theta1_0, theta2_0, d_theta2_0, w0, L, m, 0, k_val, gamma, w)
    theta1_ks[count] = curr_theta1
    theta2_ks[count] = curr_theta2
    E1_ks[count] = 0.5*curr_dtheta1**2 + 0.5*w0**2 * curr_theta1**2
    E2_ks[count] = 0.5*curr_dtheta2**2 + 0.5*w0**2 * curr_theta2**2
    count += 1


count = 0
for g_val in gs:
    curr_w = np.sqrt(w0**2 - gs[count]**2/2)
    curr_theta1, curr_theta2, curr_dtheta1, curr_dtheta2 = numerical_Euler(t_s, theta1_0, d_theta1_0, theta2_0, d_theta2_0, w0, L, m, 0, k, g_val, curr_w)
    theta1_gs[count] = curr_theta1
    theta2_gs[count] = curr_theta2
    E1_gs[count] = 0.5*curr_dtheta1**2 + 0.5*w0**2 * curr_theta1**2
    E2_gs[count] = 0.5*curr_dtheta2**2 + 0.5*w0**2 * curr_theta2**2
    count += 1





## Effetcts of varying k, gamma=0.15

"""
for i in range(len(ks)):
    plt.plot(t_s, theta1_ks[i], label=f"theta 1 k={ks[i]}")
    plt.plot(t_s, theta2_ks[i], label=f"theta 2 k={ks[i]}")
    plt.title(f"Theta(t) w={w:.3f}, dt={dt}, w0={w0:.3f}, gamma={gamma}")
    plt.xlabel("t [s]"); plt.ylabel("theta"); plt.legend(); plt.show()


for i in range(len(ks)):
    plt.plot(t_s, E1_ks[i], label=f"theta1 k={ks[i]}")
    plt.plot(t_s, E2_ks[i], label=f"theta2 k={ks[i]}")
    plt.title(f"Theta E(t) w={w:.3f}, dt={dt}, w0={w0:.3f}, gamma={gamma}")
    plt.xlabel("t [s]"); plt.ylabel("E(t)"); plt.legend(); plt.show()
"""




## Effects of varying gamma, k=0.8
"""
for i in range(len(gs)):
    curr_w = np.sqrt(w0**2 - gs[i]**2/2)
    plt.plot(t_s, theta1_gs[i], label=f"th1, gm={gs[i]}, w={curr_w:.4f}")
    plt.plot(t_s, theta2_gs[i], label=f"th2, gm={gs[i]}, w={curr_w:.4f}")
    plt.title(f"Theta(t) dt={dt}, w0={w0:.3f}, k={k}")
    plt.xlabel("t [s]"); plt.ylabel("theta"); plt.legend(); plt.show()


for i in range(len(gs)):
    curr_w = np.sqrt(w0**2 - gs[i]**2/2)
    plt.plot(t_s, E1_gs[i], label=f"th1, gm={gs[i]}, w={curr_w:.4f}")
    plt.plot(t_s, E2_gs[i], label=f"th2, gm={gs[i]}, w={curr_w:.4f}")
    plt.yscale("log"); plt.title(f"Theta E(t) dt={dt}, w0={w0:.3f}, k={k}")
    plt.xlabel("t [s]"); plt.ylabel("E(t)"); plt.legend(); plt.show()

"""






#####################################################################

"""
for i in range(len(ks)):
    plt.plot(t_s, theta1_ks[i], label=f"k={ks[i]}")
plt.title(f"Theta1 w={w:.3f}, dt={dt}, w0={w0:.3f}, gamma={gamma}")
plt.xlabel("t [s]"); plt.ylabel("theta"); plt.legend(); plt.show()


for i in range(len(ks)):
    plt.plot(t_s, E1_ks[i], label=f"k={ks[i]}")
plt.title(f"Theta1 E(t) w={w:.3f}, dt={dt}, w0={w0:.3f}, gamma={gamma}")
plt.xlabel("t [s]"); plt.ylabel("E(t)"); plt.legend(); plt.show()

"""


"""
for i in range(len(gs)):
    curr_w = np.sqrt(w0**2 - gs[i]**2/2)
    plt.plot(t_s, theta1_gs[i], label=f"gm={gs[i]}, w={curr_w:.4f}")
plt.title(f"Theta1 dt={dt}, w0={w0:.3f}, k={k}")
plt.xlabel("t [s]"); plt.ylabel("theta"); plt.legend(); plt.show()
"""


"""
for i in range(len(gs)):
    curr_w = np.sqrt(w0**2 - gs[i]**2/2)
    plt.plot(t_s, E1_gs[i], label=f"gm={gs[i]}, w={curr_w:.4f}")
plt.yscale("log")
plt.title(f"Theta1 E(t) dt={dt}, w0={w0:.3f}, k={k}")
plt.xlabel("t [s]"); plt.ylabel("E(t)"); plt.legend(); plt.show()
"""


#####################################################################





########## e) beat frequency with no external driving 2 ##########

## Varying k, gamma=0

ks = [0.4, k, 1.2]

theta1_ks = [0]*len(ks); E1_ks = [0]*len(ks)


"""
count = 0
for k_val in ks:
    curr_theta1, curr_theta2, curr_dtheta1, curr_dtheta2 = numerical_Euler(t_s, theta1_0, d_theta1_0, theta2_0, d_theta2_0, w0, L, m, 0, k_val, 0, w0)
    theta1_ks[count] = curr_theta1
    E1_ks[count] = 0.5*curr_dtheta1**2 + 0.5*w0**2 * curr_theta1**2
    count += 1


for i in range(len(ks)):
    color_lst = ["blue", "red", "green"]
    plt.plot(t_s, theta1_ks[i], label=f"k={ks[i]}", color=color_lst[i])
    plt.title(f"Theta1 w=w0, dt={dt}, w0={w0:.3f}, gamma=0")
    plt.xlabel("t [s]"); plt.ylabel("theta"); plt.legend(); plt.show()
"""



##########################################################


########## d) Observation Conclusion ##########

## 1) Varying k, gamma = 0.15 
"""
K = 0 --> Single driven Theta 1 --> Resonance amplitude at maximum
As K increases --> Resonance amplitude decreases
"""


## 2) Varying gamma, k = 0.8
"""
As gamma = 0 --> amplitude increases infinitely
As gamma increases --> amplitude decreases
"""



########## e) Observation Conclusion ##########

## 1) Varying k, gamma = 0
"""
Small K --> smaller beat frequency
Each beat occurs longer, but less frequently

As K gets larger --> larger beat frequency
Each beat occurs shorter, but more frequently
"""



## 2) Varying gamma, k=0.8
"""
Small gamma --> beat pattern decays slowly, and beat amplitude is large

Large gamma --> beat pattern decays quickly, and beat amplitude is small
"""


