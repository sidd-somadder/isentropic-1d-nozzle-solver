import numpy as np;
from area_mach_relation import *
from newton_raphson import newton_raphson_solver as nrs;
np.set_printoptions(formatter={"float": lambda x: f"{x:.5g}"})

def forward_check(tolerance=5e-5): 
    '''
    Inputs Mach number sweep into the area_mach_relation method and compares output to table data
    Assumes specific heat ratio gamma = 1.4
    Reference values from Table B1 of Gas Dynamics, 3rd Edition by James E.A. John and Theo G. Keith
    '''
    table_mach = np.array([0.2,0.5,0.8,0.96,1.50,2.00,3.00,5.00])
    table_area_ratio = np.array([2.9635,1.3398,1.0382,1.0014,1.1762,1.6875,4.2346,25.000])
    solved_area_ratio = np.sqrt(area_mach_relation(table_mach))

    abs_error = np.abs(solved_area_ratio-table_area_ratio);
    max_error = np.max(abs_error)
    print(f"Forward solution check (tolerance = {tolerance}):") 
    if max_error <= tolerance:
        print(f"PASS; max error =       {max_error:.4e}")
    else: 
        print(f"FAIL; max error =   {max_error}")
    print(f"Mach Numbers used:      {table_mach}")
    print(f"Tabulated A/A*:         {table_area_ratio}")
    print(f"Solved A/A*:            {solved_area_ratio}")
    # print(f"Errors: {abs_error}")

def finite_difference_test(M_max=5, M_min=0.01, dM=1e-6, gamma=1.4, tolerance=1e-6):
    '''
    Inputs a sweep of Mach numbers and calculates the area ratio and derivative at that mach number
    The resulting derivative array is compared to np.gradient applied on the area ratio array
    Since the derivative is blows up as M is close to 0, max relative error in max norm is reported rather than absolute error.
    Assumes specific heat ratio gamma = 1.4
    '''
    print(f"Finite difference test: all inputs nonzero; tolerance = {tolerance}")
    step_count = round((M_max-M_min)/dM) +1
    M_sweep = np.linspace(M_min,M_max,step_count, endpoint=True)
    area_results = area_mach_relation(M_sweep, gamma)
    deriv_results = area_mach_df(M_sweep, gamma)

    # Max relative error in max norm used since bad behaviour as M -> 0
    deriv_numerical = np.gradient(area_results,M_sweep, edge_order=2)
    deriv_error = np.abs(deriv_numerical-deriv_results)
    max_deriv =np.max(np.abs(deriv_results))
    max_rel_error = np.max(deriv_error / max_deriv)

    if max_rel_error <= tolerance:
        print(f"PASS; max relative error =      {max_rel_error:.4e}")
    else: 
        print(f"FAIL; max relative error =      {max_rel_error:.4e}")
    print(f"Mach sweep:                     {M_min} to {M_max}, {step_count} points (dM = {dM:.0e})")
    print(f"Specific Heat Ratio (gamma):    {gamma}")

def inverse_check(table_rounding=5e-5, gamma=1.4):
    '''
    Inputs tabulated area ratios into the Newton-Raphson solver and compares the recovered
    Mach number (on the matching branch) to the tabulated Mach number.
    Each case is judged against the error that 4-decimal table rounding alone would cause,
    since dM/d(A/A*) grows large near M = 1.
    Assumes specific heat ratio gamma = 1.4
    Reference values from Table B1 of Gas Dynamics, 3rd Edition by James E.A. John and Theo G. Keith
    '''
    table_mach = np.array([0.2,0.5,0.8,0.96,1.50,2.00,3.00,5.00])
    table_area_ratio = np.array([2.9635,1.3398,1.0382,1.0014,1.1762,1.6875,4.2346,25.000])
 
    solved_mach = np.zeros(len(table_mach))
    iterations = np.zeros(len(table_mach), dtype=int)
    for i in range(len(table_mach)):
        result = nrs(table_area_ratio[i], gamma=gamma)
        # pick the regime matching the tabulated Mach number; 
        # substract 1 from iterations since initial guess is included
        if table_mach[i] < 1:
            solved_mach[i] = result["subsonic_mach"]
            iterations[i] = len(result["subsonic_log"]) - 1
        else:
            solved_mach[i] = result["supersonic_mach"]
            iterations[i] = len(result["supersonic_log"]) - 1
 
    abs_error = np.abs(solved_mach - table_mach)
 
    # expected error from table rounding: |dM/d(A/A*)| * rounding, with d(A/A*)/dM = f'/(2*sqrt(f))
    dA_dM = area_mach_df(table_mach, gamma) / (2*np.sqrt(area_mach_relation(table_mach, gamma)))
    expected_error = np.abs(table_rounding / dA_dM)
 
    passed = abs_error <= expected_error
    print(f"Inverse Newton-Raphson check (table rounding = {table_rounding}):")
    if np.all(passed):
        print(f"PASS; all {len(table_mach)} cases within expected rounding error")
    else:
        print(f"FAIL; {np.sum(~passed)} of {len(table_mach)} cases exceed expected rounding error")
    print(f"Tabulated A/A* (input): {table_area_ratio}")
    print(f"Tabulated Mach:         {table_mach}")
    print(f"Solved Mach:            {solved_mach}")
    print(f"Errors:                 {abs_error}")
    print(f"Expected errors:        {expected_error}")
    print(f"Iterations:             {iterations}")
    
print(40*"~")
forward_check()
print(40*"~")
finite_difference_test()
print(40*"~")
inverse_check()