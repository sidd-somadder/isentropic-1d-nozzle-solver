import numpy as np
import matplotlib.pyplot as plt
from newton_raphson import newton_raphson_solver as nrs
from area_mach_relation import area_mach_relation as amr

def plot_area_mach_regime(area_ratio, subsonic, supersonic, area_result, mach_sweep):
    '''
    Plots forward calculated values of area ratio from Mach number sweep.
    Plots Newton-Raphson solved Mach numbers from area ratio sweep on top.
    Plot used as leading image of readme.md 
    '''
    fig, ax = plt.subplots()
    ax.plot(area_result, mach_sweep, color='black', lw=1.5, label="Exact Relation")
    ax.plot(area_ratio, subsonic, "o", color='steelblue', markersize=5, label="Newton-Raphson (Subsonic)")
    ax.plot(area_ratio, supersonic, "o",color='red', markersize=5, label="Newton-Raphson (Supersonic)")
    ax.set_xlabel(r"Area Ratio $\left(\frac{A}{A^*}\right)$")
    ax.axvline(x=1, color='grey', lw= 1, linestyle="--")
    ax.axhline(y=1, color='grey', lw= 1, linestyle="--")
    ax.set_ylabel(r"Mach Number $M$")
    ax.set_title(r"Nozzle Area-Mach Relation ($\gamma = 1.4$)")
    ax.set_xlim([0,10])
    ax.legend()
    ax.grid(True, alpha=0.4)
    plt.show()

def plot_residual_roots(mach, f, area_target, gamma=1.4):
    '''
    Plots the residual function f_A as discussed in technical report to demonstrate that there are 4 solutions
    '''
    res = nrs(area_ratio_target=area_target, gamma=gamma)
    roots = [-res["supersonic_mach"], -res["subsonic_mach"],
             res["subsonic_mach"], res["supersonic_mach"]]

    fig, ax = plt.subplots()
    ax.axvspan(mach.min(), 0, color="grey", alpha=0.15)
    ax.plot(mach, f, color="black", lw=1.5, label=r"$f_A(M)$")
    ax.axhline(0, color="grey", lw=1, ls="--")
    ax.plot(roots, np.zeros(4), "o", color="tab:red", ms=6, zorder=3, label="Roots")
    ax.set_ylim(-5, 8)
    ax.set_xlabel(r"Mach Number $M$")
    ax.set_ylabel(r"$f_A(M) = \left(\frac{A}{A^*}\right)^2 - \left(\frac{A}{A^*}\right)^2_{\mathrm{target}}$")
    ax.set_title(rf"Residual Roots for area ratio = {area_target} ($\gamma = {gamma}$)")
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.4)
    plt.show()

# 
N = 41
M = 501
K = 200

mach_theoretical = np.union1d(np.linspace(-3,0, K, endpoint=False),np.linspace(0.015,3,K))
area_sweep = np.union1d(np.linspace(1,10,N, endpoint=True), np.array([1.02, 1.08]))
mach_sweep = np.linspace(0.01,4,M, endpoint=True)

area_result = np.zeros(M)
supersonic_result = np.zeros(N+2)
subsonic_result = np.zeros(N+2)
f = np.zeros(2*K)

for j in range(2*K):
    f[j] = amr(mach_theoretical[j], 1.4) - 2**2

for i in range(N+2):
    result = nrs(area_ratio_target=area_sweep[i])
    subsonic_result[i] = result["subsonic_mach"]
    supersonic_result[i] = result["supersonic_mach"]

for k in range(M):
    area_result[k] =  np.sqrt(amr(mach_sweep[k]))

plot_area_mach_regime(area_sweep, subsonic_result, supersonic_result, area_result, mach_sweep)
plot_residual_roots(mach_theoretical,f,2)

