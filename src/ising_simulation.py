import numpy as np
import matplotlib.pyplot as plt
from numba import jit
from tqdm import tqdm
from astropy.stats import sigma_clipped_stats
from scipy import special
import os

random = np.random.default_rng()
plt.rcParams["axes.labelsize"] = 16
saveopts = dict(bbox_inches="tight",pad_inches=0.1)
prefix = "./out/"
a = 0.3 # global alpha control

os.makedirs(prefix,exist_ok=True)


@jit("i8(i8[:,:],i8,i8)")
def local_energy(grid,i,j):
    """Calculates the energy around the spin on the i,j position on the grid.
    This function is used to calculate the transition probability to avoid innecesary 
    recalculation of the local grid"""
    # periodic boundary conditions
    N = len(grid); ip1=(i+1)%N; im1=(i-1)%N; jp1=(j+1)%N; jm1=(j-1)%N;
    E_local = grid[i,j] * (grid[ip1,j]+grid[im1,j]+grid[i,jp1]+grid[i,jm1])
    return E_local

@jit("i8(i8[:,:])")
def total_energy(grid):
    """Total energy of the grid using the Ising Hamiltonian. 
    Only used once at the beginning of the simulation"""
    n, m = grid.shape
    E = 0
    for i in range(n):
        for j in range(m):
            E -= grid[i, j] * (grid[(i-1)%n, j]+grid[(i+1)%n,j]+grid[i,(j-1)%m]+grid[i,(j+1)%m])
    return E


@jit("f8[:](i8[:,:],f8)")
def update_grid(grid,β):
    """A single step in the Metropolis algorithm.
    1) Selects an i,j position on the grid
    2) Calculates the change of energy and the flip probability
    3) Shoots a random number to decide if the flips happens
    4) Updates the grid in place
    5) returns the change in energy and the change in magnetization as a 2D array"""
    i,j = np.random.randint(0,len(grid),2)
    ΔE = 2.0 * local_energy(grid,i,j)
    if ΔE <= 0:
        flip = True
    else:
        A = np.exp(-β*ΔE)
        u = np.random.rand()
        flip = u < A # wheter to flip or not
    if flip:
        grid[i,j] = -grid[i,j]
        return np.array([ΔE, 2*grid[i,j]])
    else: 
        return np.array([0., 0.])

@jit("f8[:,:](i8[:,:],f8,i8,i8)")
def simulate_and_measure(grid,β,steps,thin):
    """Simulates the Ising square system given an initial grid 
    and temperature for the given number of steps.

    Measure the energy and magnetization each `thin` steps,
    returns a 2 x (steps//thin) matrix with the results.

    This function updates the grid in-place, so the provided 
    initial matrix is modified to show the system at the last iteration.
    
    The Energy and the magnetization are in microscopic units, 
    i.e. not normalized to the total number of spins.
    """
    results = np.zeros((steps//thin,2))
    E = 0.0
    M = grid.sum()
    for i in range(steps):
        ΔE,ΔM = update_grid(grid,β)
        M += ΔM
        E += ΔE
        if i%thin or thin==1:
            results[i//thin,0] = E
            results[i//thin,1] = M
    return results

@jit
def simulate_ising(initial_grid,β,steps,thin):
    """Convenience function around `simulate_and_measure`. 
    Copies the grid instead of modifiying it in place, 
    and returns the results as individual arrays: 
        energy in the correct units, 
        magnetization in the correct units, 
        and the final grid"""
    grid = initial_grid.copy()
    Ei = total_energy(grid)
    Es,Ms = simulate_and_measure(grid,β,steps,thin).T
    Es += Ei
    Es /= 2*len(grid)**2
    Ms /= len(grid.ravel())
    return Es,Ms,grid

def teo_M(β):
    "Theoretical solution for the mean magnetization"
    β = np.asarray(β)
    M = abs(1-np.sinh(2*β)**(-4))**(1/8)
    return np.where(β<=teo_βc,0.0,M)

def teo_E(β):
    "Theoretical solution for the mean energy"
    m = 4 * np.sinh(2*β)**2 / np.cosh(2*β)**4
    return -1/np.tanh(2*β) * (1 + (2/np.pi) * (2*np.tanh(2*β)**2 - 1) * special.ellipk(m))

teo_βc = np.asinh(1)/2

# Show iteration effect on measures

plt.figure()
N = 30
β = 0.8
n_iter = int(2e5)
x_iter = np.arange(n_iter)
spins = random.choice(np.array([-1,1],dtype=int),(N,N))
Es,Ms,spins = simulate_ising(spins,β,int(2e5),1)

plt.plot(Es,label="Mean energy")
plt.plot(Ms,label="Mean magnetization")
plt.xlabel("Iteration number")
plt.ylabel("Measure")
plt.legend(loc=2)

thin = 20
for i in range(20):
    spins = random.choice(np.array([-1,1],dtype=int),(N,N))
    Es,Ms,spins = simulate_ising(spins,β,int(2e5),thin)
    plt.plot(x_iter[::thin],Es,c="C0",alpha=a/2)
    plt.plot(x_iter[::thin],Ms,c="C1",alpha=a/2)
plt.xscale("asinh",linear_width=1000)
plt.xticks([1e3,1e4,1e5]);
plt.savefig(prefix+"iterations.pdf",**saveopts)
plt.close()


# Show planned iterations agains beta

betas = np.linspace(0.1,0.9,20)

def max_iter(β):
    """As the system takes more iterations to reach equilibrium near the critical point,
    We use this function to decide how many iterations we will perform depending on the temperature"""
    return np.round(500_000 * (6-np.log(np.cosh(13*(β-teo_βc))))).astype(int)

plt.figure()
plt.plot(betas,max_iter(betas))
plt.yscale("log")
plt.xlabel("$\\beta$")
plt.ylabel("Planned iterations (empirical)")
plt.savefig(prefix+"iteration_plan.pdf",**saveopts)
plt.close()


# Simulate

N = 30
betas = np.linspace(0.1,0.9,100)
spins = random.choice(np.array([-1,1],dtype=int),(N,N))

mean_energy = np.zeros_like(betas)
mean_magnet = np.zeros_like(betas)
var_energy = np.zeros_like(betas)
var_magnet = np.zeros_like(betas)
for i,β in tqdm(list(enumerate(betas))):
    Es,Ms,spins = simulate_ising(spins,β,max_iter(β),1)
    # Energy statistics
    mean,median,std = sigma_clipped_stats(Es[100:],sigma=5)
    mean_energy[i] = mean
    var_energy[i] = std**2
    # magnetization statistics
    mean,median,std = sigma_clipped_stats(Ms[100:],sigma=5)
    mean_magnet[i] = mean
    var_magnet[i] = std**2

# Save CV plot 

db = 1e-5
dEdb = (teo_E(betas-db)-teo_E(betas+db))/(2*db)
plt.figure()
plt.ylabel("Specific heat")
plt.xlabel("$\\beta$")
plt.plot(betas,(betas/2)**2*dEdb,c='r',alpha=a,lw=3,label="Theory")
plt.axvline(teo_βc,c='b',alpha=0.4,zorder=-10,label="Critical point")
plt.plot(betas,betas**2*var_energy*N**2,c='k')
plt.legend()
plt.savefig(prefix+"Cv.pdf",**saveopts)
plt.close()


# Save E plot 

plt.figure()
plt.errorbar(betas,mean_energy,np.sqrt(var_energy),fmt=" k",capsize=3,alpha=0.3,zorder=-10)
plt.scatter(betas,mean_energy,s=1,c='k',zorder=10)
plt.plot(betas,teo_E(betas),c='r',alpha=a,lw=3,label="Theory")
plt.xlabel("$\\beta$"); plt.ylabel("Mean energy")
plt.axvline(teo_βc,c='b',alpha=0.4,zorder=-10,label="Critical point")
plt.legend(loc=0)
plt.savefig(prefix+"E.pdf",**saveopts)
plt.close()


# Save M plot

plt.figure()
plt.errorbar(betas,mean_magnet/np.sign(mean_magnet[-1]),np.sqrt(var_magnet),fmt=" k",capsize=3,alpha=0.3,zorder=-10)
plt.scatter(betas,mean_magnet/np.sign(mean_magnet[-1]),s=4,c='k')
plt.plot(betas,teo_M(betas),c='r',alpha=a,lw=3,label="Theory")
plt.axvline(teo_βc,c='b',alpha=0.4,zorder=-10,label="Critical point")
plt.xlabel("$\\beta$"); plt.ylabel("Mean magnetization")
plt.legend(loc=0)
plt.savefig(prefix+"m.pdf",**saveopts)
plt.close()
