import numpy as np
import matplotlib.pyplot as plt


def f_periodic_2pi(x):
    f = 0
    if x >= -np.pi/2 and x <= np.pi: f = 1
    return f


def f_periodic_2pi_periods(N, periods, interval):
    fxs = np.zeros(periods*(N+1))
    xs_final = np.zeros(periods*(N+1))
    #jumps = [0]*periods
    
    x_step = 0
    for per in range(periods):
        curr_xs = np.linspace(interval[0] + per*2*np.pi, interval[1] + per*2*np.pi, N+1)
        for x in curr_xs:
            xs_final[x_step] = x
            fx = f_periodic_2pi(x - per*2*np.pi)
            fxs[x_step] = fx
            x_step += 1
    
    return xs_final, fxs


def plot_f_periodic_2pi_periods(xs_, fxs_, periods, lnstyle="solid"):
    plt.plot(xs_, fxs_, linestyle=lnstyle)
    plt.title(f"2pi-periodic {periods} periods")
    plt.xlabel("x"); plt.ylabel("f(x)"); plt.show()






## 2a)


def plot_f_periodic_2pi_periods_complete(N_t, periods, interval):
    jumps = [0]*(2*periods+1)    
    x_jumps = [0]*len(jumps); x_jumps[0] = interval[0]
    N = int(N_t/2)
    
    
    jump_cnt = 1
    for per in range(periods):
        x_step = 0
        curr_xs_1 = np.linspace(interval[0] + per*2*np.pi, interval[0] + np.pi/2 + per*2*np.pi, N+1)
        curr_xs_2 = np.linspace(interval[0] + np.pi/2 + per*2*np.pi, interval[1] + per*2*np.pi, N+1)
        curr_fxs_1 = np.zeros(N+1); curr_fxs_2 = np.zeros(N+1)
        
        x_jumps[jump_cnt] = interval[0] + np.pi/2 + per*2*np.pi
        x_jumps[jump_cnt+1] = interval[1] + per*2*np.pi
        jump_cnt += 2
        
        for x in curr_xs_2:
            fx = f_periodic_2pi(x - per*2*np.pi)
            curr_fxs_2[x_step] = fx
            x_step += 1
        
        plt.plot(curr_xs_1, curr_fxs_1)
        plt.plot(curr_xs_2, curr_fxs_2)
    
    for i in range(len(jumps)):
        plt.plot([x_jumps[i], x_jumps[i]], [0, 1], color="black", linestyle="dashed")
    
    plt.title(f"2pi-periodic {periods} periods")
    plt.xlabel("x"); plt.ylabel("f(x)"); plt.show()
    
    
    return 
    





## 2b)

def a_n(n):
    a = 1/(n*np.pi) * ((-1)**n - 1)/2 * (-1)**(int((n+1)/2))
    
    return a


def b_n(n):
    b = 1/(n*np.pi) * ((-1)**(n+1) + ((-1)**(n+1) - 1)/2 * (-1)**(int((n+2)/2)) )
    
    return b


def periodic_2pi_Fourier(N, x_s):
    a0 = 3/2
    an_s = np.zeros(N); bn_s = np.zeros(N); cn_s = np.zeros(N)
    ffx_s = np.zeros(len(x_s))
    for n in range(1, N+1):
        an_s[n-1] = a_n(n)
        bn_s[n-1] = b_n(n)
        cn_s[n-1] = np.sqrt(an_s[n-1]**2 + bn_s[n-1]**2)
    
    
    for i in range(len(x_s)):
        curr_sum = 0; curr_x = x_s[i]
        for n in range(1, N+1):
            curr_sum += an_s[n-1] * np.cos(n*curr_x) + bn_s[n-1] * np.sin(n*curr_x)
        ffx_s[i] = a0/2 + curr_sum
    
    print(bn_s)
    
    return ffx_s, cn_s




if __name__ == "__main__":
    X = 100
    periods = 3
    interval = [-np.pi, np.pi]
    #N = 40
    
    Ns = [10, 100]
    
    
    ## a) ##
    #xs_, fxs_ = f_periodic_2pi_periods(X, periods, interval)
    #plot_f_periodic_2pi_periods_complete(X, periods, interval)
    
    
    ## b) d) ##
    periods = 1
    xs_, fxs_ = f_periodic_2pi_periods(X, periods, interval)
    
    
    plt.plot(xs_, fxs_, label="f(x)", linestyle="dashed")
    for N in Ns:
        ffx_s, cn_s = periodic_2pi_Fourier(N, xs_)
        plt.plot(xs_, ffx_s, label=f"N={N}")
    plt.xlabel("x"); plt.ylabel("f(x)")
    plt.legend(); plt.show()
        
    ## c) ##
    for N in Ns:
        ffx_s, cn_s = periodic_2pi_Fourier(N, xs_)
        #print(cn_s)
        n_s = np.arange(1, N+1)
        plt.plot(n_s, cn_s, label=f"N={N}")
        plt.xlabel("n"); plt.ylabel("c_n")
        plt.legend(); plt.show()
    
    


