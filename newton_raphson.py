import numpy as np
import area_mach_relation as amr;

# Note: for now use predefined steps, update to tolerance based iteration after confirming initial setup
def newton_raphson_solver(area_ratio_target, steps=5, gamma=1.4):

    # initial guesses for both subsonic and supersonic root of area-to-mach relation
    Msub = 0.6
    Msuper = 1.4

    # initialize array of saved step values for each root; includes initial guess
    Msub_steps = np.zeros(steps+1)
    Msub_steps[0] = Msub

    Msuper_steps = np.zeros(steps+1)
    Msuper_steps[0] = Msuper

    # run iteration for each initial guess to approximate each root:
    for i in range(steps):  # subsonic
        Msub = next_guess(Msub, area_ratio_target, gamma)
        Msub_steps[i+1] = Msub
    for i in range(steps):  # supersonic
        Msuper = next_guess(Msuper, area_ratio_target, gamma)
        Msuper_steps[i+1] = Msuper

    # returns dictionary of final approximations and array log of each step to get there
    return {"M_sub": Msub, "M_super": Msuper, "NR_steps_sub": Msub_steps, "NR_steps_super": Msuper_steps}

def next_guess(M_guess, A_target, gamma):
    # get residual function and derivative values from area_mach_relation.py
    f = amr.area_mach_relation(M_guess, gamma) - A_target**2
    df = amr.area_mach_df(M_guess, gamma)

    # tangent line approximation step    
    h = f / df

    return M_guess - h


    