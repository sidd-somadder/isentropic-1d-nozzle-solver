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

def get_stag_properties(mach, gamma=1.4):
    '''
    Uses mach number and isentropic flow properties to return stagnation temperature, pressure, density, and acoustic speed ratios
    Returns dictionary of isentropic flow properties:
    stag_temp_ratio is T0/T
    stag_density_ratio is rho0/rho
    stag_pressure_ratio is p0/p
    stag_acoustic_ratio is a0/a
    '''
    T_stag = 1 + ((gamma-1)/2)*(mach**2)
    p_stag = T_stag**(gamma/(gamma-1))
    rho_stag = T_stag**(1/(gamma-1))
    a_stag = T_stag**(0.5)

    return {"stag_temp_ratio": T_stag, "stag_density_ratio": rho_stag, "stag_pressure_ratio": p_stag, "stag_acoustic_ratio": a_stag}