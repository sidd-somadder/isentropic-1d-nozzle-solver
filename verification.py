import numpy as np;
from newton_raphson import newton_raphson_solver as nrs;

def target_area(A):
    results = nrs(area_ratio_target=A);
    print(f"M_subsonic = {results["subsonic_mach"]}");
    print(f"M_supersonic = {results["supersonic_mach"]}");
    print(f"M_subsonic LOG = {results["subsonic_log"]}");
    print(f"M_supersonic LOG = {results["supersonic_log"]}");

target_area(10);