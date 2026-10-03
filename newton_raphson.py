import numpy as np
import area_mach_relation as amr

def newton_raphson_solver(area_ratio_target, max_steps=100, gamma=1.4, tolerance= 1e-8):
    '''
    Uses the target area ratio and returns the subsonic and supersonic Mach number solution to the area-to-mach relation
    max_steps and tolerance used in newton_iteration(...) to specify max run time for nonconverging solutions and tolerance for convergence.
    specific heat ratio (gamma) assumed 1.4 for air, but may be changed.
    '''
    # Handle sonic flow and no solution cases separately where numerical approximation not applicable
    if area_ratio_target == 1:
        return {"subsonic_mach": 1.0, "supersonic_mach": 1.0,
                "NR_steps_sub": np.array([1.0]), "NR_steps_super": np.array([1.0])}
    if area_ratio_target < 1:
        raise ValueError("A/A* < 1 has no solution.")

    # initial guesses for both subsonic and supersonic root of area-to-mach relation
    M_sub, M_super = get_initial_guess(area_ratio_target, gamma)

    # run iteration for each initial guess to approximate each root within tolerance:
    M_sub, Msub_steps = newton_iteration(M_sub, max_steps, area_ratio_target, gamma, tolerance)
    M_super, Msuper_steps = newton_iteration(M_super, max_steps, area_ratio_target, gamma, tolerance)

    # convert python list quantities into numpy array
    Msub_log = np.array(Msub_steps)
    Msuper_log = np.array(Msuper_steps)

    # returns dictionary of final approximations and array log of each step to get there
    return {"subsonic_mach": M_sub, "supersonic_mach": M_super, 
            "subsonic_log": Msub_log, "supersonic_log": Msuper_log}

def newton_iteration(M_initial, max_steps, A_target, gamma, tolerance):
    '''
    Uses initial Mach number guess, target area ratio, specific heat ratio (gamma), and tolerance to iterate until area-to-Mach relation root is found
    max_steps specified against diverging Newton-Raphson runs; raises RuntimeError if no convergence in those number of steps
    '''
    M_guess = M_initial
    M_log = []
    M_log.append(M_initial)

    for i in range(max_steps):
        Mprev = M_guess
        M_guess = next_guess(M_guess, A_target, gamma)
        M_log.append(M_guess)
    
        if np.abs(M_guess - Mprev) < tolerance:
            return M_guess, M_log

    raise RuntimeError(f"No convergence within {max_steps} iterations from initial guess M0 = {M_initial} and tolerance = {tolerance}")

def next_guess(M_guess, A_target, gamma):
    '''
    Calculates next Mach number approximation using the formula M_{i+1} = M_i - f(M_i)/f'(M_i)
    Uses current Mach guess (M_guess), target area ratio, and specific heat ratio (gamma) to calculate f(M_guess) and f'(M_guess) 
    Raises error if M = 1, and so first derivative is 0
    '''
    # get residual function and derivative values from area_mach_relation.py
    f = amr.area_mach_relation(M_guess, gamma) - A_target**2
    df = amr.area_mach_df(M_guess, gamma)

    if df == 0:
        raise ZeroDivisionError("Derivative is zero; Newton-Raphson reached M = 1.")
    
    # tangent line approximation step    
    h = f / df

    return M_guess - h

def get_initial_guess(area_ratio, gamma):
    '''
    Use limiting behaviour of Mach number and specifc heat ratio (gamma) in area-to-Mach relation to calculate good initial guesses
    Returns two initial guesses: one subsonic and one supersonic 
    '''
    C = 2 / (gamma+1)
    n = (gamma+1)/(gamma-1)

    M0_sub = C**(n/2) / area_ratio
    M0_super = (n**(n/2) * area_ratio)**((gamma-1)/2)

    return M0_sub, M0_super;


    