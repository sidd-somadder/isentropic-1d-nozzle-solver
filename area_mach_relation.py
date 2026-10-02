import numpy as np

def area_mach_relation(M, gamma=1.4):
    '''
    Returns the square of the area ratio for a 1D isentropic nozzle using input mach number M
    Assumes M > 0
    Specific heat ratio (gamma) assumed 1.4 for air, can be changed
    '''
    inv_M2 = 1 / (M**2)
    n = (gamma+1)/(gamma-1)
    C = 2 / (gamma+1)
    R = 1 + ((gamma-1)/2) * (M**2)
    return inv_M2 * (C * R)**n

def area_mach_df(M, gamma= 1.4):
    '''
    Used in newton_raphson.py to return the first derivative of f(M) = (A/A*)**2 - target for input mach number M
    Assumes M > 0
    Specific heat ratio (gamma) assumed 1.4 for air, can be changed
    '''
    n = (gamma+1)/(gamma-1)
    C = 2 / (gamma+1)
    R = 1 + ((gamma-1)/2) * (M**2)

    return (C**n) * (M**(-3)) * (R**(n-1)) * 2*(M**2 - 1)