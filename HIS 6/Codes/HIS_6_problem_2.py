import numpy as np
import matplotlib.pyplot as plt


#### Probability density function ####
def P(theta):
    """ theta : array or scalar"""
    
    #---- Option 1 ----#
    #p = (np.cos(0.5*k*d*np.sin(theta)))**2
    
    
    
    #---- Option 2 ----#
    """Better suited"""
    p = (np.cos(0.5*np.pi*d*np.sin(theta)/lmbda))**2
    
    return p



#### List of probability density ####
def create_P(theta_lst):
    P_lst = P(theta_lst)
    return P_lst



#### Plot of probability density ####
def plot_P(theta_lst):
    """ theta_lst : list of radians, with size N+1 """
    P_lst = create_P(theta_lst)
    
    plt.plot(theta_lst, P_lst, color="r")
    plt.title("Electron probability density angular")
    plt.xlabel("theta [radians]"); plt.ylabel("P(theta)")
    plt.show()
    return



#### Convert degree to radians ####
def convert_deg_rad(deg_lst, steps):
    rad_start, rad_end = deg_lst[0]*np.pi/180, deg_lst[1]*np.pi/180
    rad_lst = np.linspace(rad_start, rad_end, steps + 1)
    return rad_lst



#---------------------------------------------------------#
                ########## ACTION ##########
#---------------------------------------------------------#
if __name__ == "__main__":
    print("Problem 2 - Probability density plot")
    
    ##### Initial conditions #####
    theta = [-5, 5]             # Degree (must convert into radians format)
    N = 100                     # steps
    
    d = 2 * 10**(-6)            # Slit distance (m)
    
    ##### Option 1 #####
    k = 6.27 * 10**(10)         # Wave number (m^(-1))
    
    
    ##### Option 2 #####
    lmbda = 10**(-10)           # Wave length (m)
    
    
    
    
    ##### Make computations #####
    new_theta = convert_deg_rad(theta, N)
    
    
    ##### Plot probability density #####
    plot_P(new_theta)
    
    
    
    ##### Words #####
    """Option 1 (using k) vs Option 2 (using lmbda): 
        Shows difference in the plots
        Option 2 is preferred"""
    