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
    M_sweep = np.array([0.2,0.5,0.8,0.96,1.50,2.00,3.00,5.00])
    table_area_ratio = np.array([2.9635,1.3398,1.0382,1.0014,1.1762,1.6875,4.2346,25.000])
    solved_area_ratio = np.sqrt(area_mach_relation(M_sweep))

    abs_error = np.abs(solved_area_ratio-table_area_ratio);
    max_error = np.max(abs_error)
    print(f"Forward solution check (tolerance = {tolerance}):") 
    if max_error <= tolerance:
        print(f"PASS; max error = {max_error:.4e}")
    else: 
        print(f"FAIL; max error = {max_error}")
    print(f"Mach Numbers used: {M_sweep}")
    print(f"Tabulated A/A*: {table_area_ratio}")
    print(f"Errors: {abs_error}")


forward_check()
print(20*"~")