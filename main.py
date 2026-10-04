import numpy as np;
from area_mach_relation import *
from newton_raphson import newton_raphson_solver as nrs;

def print_results(A,gamma):
    results = nrs(area_ratio_target=A,gamma=gamma);
    subsonic_stag = get_stag_properties(mach=results["subsonic_mach"])
    supersonic_stag = get_stag_properties(mach=results["supersonic_mach"])
    print(f"For nozzle area ratio = {A}:")
    print(15*"~")
    print(f"SUBSONIC Mach = {results["subsonic_mach"]:.4f}")
    print(f"Converged in {len(results["subsonic_log"])-1} iterations");
    print(f"   * T0/T = {subsonic_stag["stag_temp_ratio"]:.4f}")
    print(f"   * P0/P = {subsonic_stag["stag_pressure_ratio"]:.4f}")
    print(f"   * rho0/rho = {subsonic_stag["stag_density_ratio"]:.4f}")
    print(f"   * a0/a = {subsonic_stag["stag_acoustic_ratio"]:.4f}")
    print(f"Subsonic solver LOG = {results["subsonic_log"]}");
    print(15*"~")
    print(f"SUPERSONIC Mach = {results["supersonic_mach"]:.4f}")
    print(f"Converged in {len(results["supersonic_log"])-1} iterations");
    print(f"   * T0/T = {supersonic_stag["stag_temp_ratio"]:.4f}")
    print(f"   * P0/P = {supersonic_stag["stag_pressure_ratio"]:.4f}")
    print(f"   * rho0/rho = {supersonic_stag["stag_density_ratio"]:.4f}")
    print(f"   * a0/a = {supersonic_stag["stag_acoustic_ratio"]:.4f}")
    print(f"Supersonic solver LOG = {results["supersonic_log"]}");
    print("Note: LOGs also contain initial guess.")

if __name__ == "__main__":
    try:
        while True:
            area_ratio = float(input("Please input your desired nozzle area ratio (A/A*): "))
            gamma = float(input("Please input your fluid specific heat ratio (calorically perfect air: gamma = 1.4): "))
            if area_ratio >= 1 and gamma > 1:
                print_results(A=area_ratio, gamma=gamma)
                break;
            else:
                print("One or more of inputs are invalid. A/A* must be greater than or equal to 1; gamma must be greater than 1. Please try again.")
    except ValueError:
        print("Nozzle area ratio must be a real number greater or equal to 1.")

