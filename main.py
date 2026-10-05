import numpy as np;
from area_mach_relation import *
from newton_raphson import newton_raphson_solver as nrs;

def print_branch(label, mach, log, gamma):
    stag = get_stag_properties(mach=mach, gamma=gamma)
    print(15 * "~")
    print(f"{label} Mach = {mach:.4f}")
    print(f"Converged in {len(log) - 1} iterations")
    print(f"   * T0/T = {stag['stag_temp_ratio']:.4f}")
    print(f"   * P0/P = {stag['stag_pressure_ratio']:.4f}")
    print(f"   * rho0/rho = {stag['stag_density_ratio']:.4f}")
    print(f"   * a0/a = {stag['stag_acoustic_ratio']:.4f}")
    print(f"{label.capitalize()} solver LOG = {log}")

def print_results(A, gamma):
    results = nrs(area_ratio_target=A, gamma=gamma)
    print(f"For nozzle area ratio = {A}, gamma = {gamma}:")
    print_branch("SUBSONIC", results["subsonic_mach"], results["subsonic_log"], gamma)
    print_branch("SUPERSONIC", results["supersonic_mach"], results["supersonic_log"], gamma)
    print("Note: LOGs also contain the initial guess.")

if __name__ == "__main__":
    while True:
        try:
            area_ratio = float(input("Please input your desired nozzle area ratio (A/A*): "))
            gamma = float(input("Please input your fluid specific heat ratio (air: gamma = 1.4): "))
        except ValueError:
            print("Inputs must be numbers. Please try again.")
            continue
        if area_ratio >= 1 and gamma > 1:
            print_results(A=area_ratio, gamma=gamma)
            break
        print("Invalid input: A/A* must be at least 1, and gamma must be greater than 1. Please try again.")
