import numpy as np
import area_mach_relation as amr

# Note: for now use predefined steps, update to tolerance based iteration after confirming initial setup
def newton_raphson_solver(area_ratio_target, max_steps=100, gamma=1.4, tolerance= 1e-8):

    # Handle sonic flow and no solution cases separately where numerical approximation not applicable
    if area_ratio_target == 1:
        return {"M_sub": 1.0, "M_super": 1.0,
                "NR_steps_sub": np.array([1.0]), "NR_steps_super": np.array([1.0])}
    if area_ratio_target < 1:
        raise ValueError("A/A* < 1 has no solution.")

    # initial guesses for both subsonic and supersonic root of area-to-mach relation
    M_sub = 0.6
    M_super = 1.4

    # run iteration for each initial guess to approximate each root within tolerance:
    M_sub, Msub_steps = newton_iteration(M_sub, max_steps, area_ratio_target, gamma, tolerance)
    M_super, Msuper_steps = newton_iteration(M_super, max_steps, area_ratio_target, gamma, tolerance)

    # convert python list quantities into numpy array
    Msub_log = np.array(Msub_steps)
    Msuper_log = np.array(Msuper_steps)

    # returns dictionary of final approximations and array log of each step to get there
    return {"M_sub": M_sub, "M_super": M_super, 
            "NR_steps_sub": Msub_log, "NR_steps_super": Msuper_log}

def newton_iteration(M_initial, max_steps, A_target, gamma, tolerance):
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
    # get residual function and derivative values from area_mach_relation.py
    f = amr.area_mach_relation(M_guess, gamma) - A_target**2
    df = amr.area_mach_df(M_guess, gamma)

    if df == 0:
        raise ZeroDivisionError("Derivative is zero; Newton-Raphson reached M = 1.")
    
    # tangent line approximation step    
    h = f / df

    return M_guess - h


    