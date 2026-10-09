#!/usr/bin/env python3
"""Independent 0D check of the coverage and the impingement rate of this test case.

Collisions are switched off (Particles-DSMC-CollisMode = 0), so the gas between the inflow
plane (BC_Xminus) and the catalytic wall (BC_Xplus) is exactly free molecular: every molecule
emitted by the surface flux reaches the catalytic wall exactly once, and molecules leaving the
wall are absorbed at BC_Xminus without returning. The impingement flux on the wall is therefore
the Hertz-Knudsen flux of the surface flux itself,

    Gamma = n * sqrt(k*T / (2*pi*m))                                        [1/(m^2 s)]

independent of what the wall does with the molecules. The coverage then follows the 0D balance
implemented in surfacemodel_chemistry.f90 (adsorption) and particle_surface_chemflux.f90
(desorption):

    dtheta/dt = DissOrder * S(theta) * Gamma / N_s
              - nu(theta) * N_s^(DissOrder-1) * theta^DissOrder * exp(-E(theta)/T_w)

    S(theta)  = S0 / (1 + K*(1/Theta - 1)),  Theta = (1 - theta/theta_max)^DissOrder
    nu(theta) = 10^(Ca + Cb*theta)   (x 1e15 for DissOrder = 2)
    E(theta)  = E_ini + theta * W
    N_s       = MolPerUnitCell / LatticeVector^2

Usage
-----
    ./verify_analytic.py                        # check the files of a finished run
    ./verify_analytic.py SurfaceAnalyze_ref.csv # check the reference files of this case
    ./verify_analytic.py path/to/SurfaceAnalyze.csv -c path/to/coverage_....h5

With no argument the newest non-reference output of the current directory is used. When a
'*_ref.csv' is given, the matching '*_ref.h5' is picked up automatically. Either file may be
absent, the corresponding check is then skipped. Requires numpy/scipy/h5py.
"""
import argparse, glob, math, os, sys
import numpy as np
import h5py
from scipy.integrate import solve_ivp

kB = 1.380649e-23

# --- set-up (must match parameter.ini / hopr.ini) --------------------------------------
NAME     = 'CO'
T_w      = 1116.8                 # Part-Boundary1-WallTemp                 [K]
lattice  = 0.389e-9               # Part-Boundary1-LatticeVector            [m]
molcell  = 1.0                    # Part-Boundary1-NbrOfMol-UnitCell        [-]
A_wall   = 400e-9 * 400e-9        # area of BC_Xplus                        [m^2]
MPF      = 1.0                    # Part-Species$-MacroParticleFactor       [-]
dt       = 1.25e-12               # ManualTimeStep                          [s]
t_end    = 10.0e-9                # tend                                    [s]
n_out    = 1600                   # Surface-AnalyzeStep                     [-]

mass     = 4.6510e-26             # Part-Species1-MassIC                    [kg]
n_in     = 7.03323e23             # Part-Species1-Surfaceflux1-PartDensity  [1/m^3]
S0       = 1.00                   # Surface-Reaction1-StickingCoefficient
K        = 0.60                   # Surface-Reaction1-EqConstant
thmax    = 0.333                  # Part-Boundary1-Species1-MaxCoverage
DO       = 1                      # Surface-Reaction$-DissOrder
Ca       = 16.0                   # Surface-Reaction2-Ca
Cb       = -15.0                  # Surface-Reaction2-Cb
E_ini    = 17688.8                # Surface-Reaction2-Energy                [K]
W        = -18410.8               # Surface-Reaction2-LateralInteraction    [K]
th0      = 0.09                   # Part-Boundary1-Species1-Coverage

# tolerance of the check against the 0D solution
TOL_COV, TOL_COLL = 0.05, 0.01

N_s  = molcell / lattice**2       # surface site density                [1/m^2]
cbar = math.sqrt(kB * T_w / (2.0 * math.pi * mass))


def sticking(th):
    th = min(max(th, 0.0), 1.0)
    Theta = min(max(1.0 - th / thmax, 0.0), 1.0) ** DO
    if Theta <= 0.0: return 0.0
    return min(max(S0 / (1.0 + K * (1.0 / Theta - 1.0)), 0.0), 1.0)


def desorption(th):
    th = min(max(th, 0.0), 1.0)
    nu = 10.0 ** (Ca + Cb * th)
    if DO == 2: nu *= 1.0e15
    return nu * N_s ** (DO - 1) * th ** DO * math.exp(-(E_ini + th * W) / T_w)


def rhs(th):
    return DO * sticking(th) * n_in * cbar / N_s - desorption(th)


def newest(pattern, want_ref):
    """Newest file matching pattern, restricted to reference or non-reference files."""
    hits = [f for f in glob.glob(pattern) if ('_ref.' in os.path.basename(f)) == want_ref]
    return max(hits, key=os.path.getmtime) if hits else None


def parse_args():
    p = argparse.ArgumentParser(description='Check the coverage and the impingement rate against a 0D solution.',
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__.split('Usage\n-----\n')[1])
    p.add_argument('surface_csv', nargs='?', default=None,
                   help='SurfaceAnalyze CSV file (default: the newest non-reference '
                        'SurfaceAnalyze*.csv in the current directory)')
    p.add_argument('-c', '--coverage', default=None, metavar='FILE',
                   help='DSMCSurfState HDF5 file holding the coverage (default: matched to '
                        'surface_csv, i.e. the *_ref.h5 file when a *_ref.csv is given)')
    a = p.parse_args()

    # A '*_ref.csv' argument selects the reference set, anything else the run output
    want_ref = a.surface_csv is not None and '_ref.' in os.path.basename(a.surface_csv)
    base     = os.path.dirname(a.surface_csv) if a.surface_csv else ''
    if a.surface_csv is None:
        a.surface_csv = newest(os.path.join(base or '.', 'SurfaceAnalyze*.csv'), want_ref)
    if a.coverage is None:
        a.coverage = newest(os.path.join(base or '.', 'coverage_DSMCSurfState_*.h5'), want_ref)
    return a


def main():
    args = parse_args()
    if not (args.surface_csv or args.coverage):
        sys.exit('verify_analytic.py: no SurfaceAnalyze*.csv and no coverage_DSMCSurfState_*.h5 '
                 'found - run the case first, or pass the files explicitly (see --help)')

    # --- reference solution ------------------------------------------------------------
    sol       = solve_ivp(lambda t, y: [rhs(y[0])], [0.0, t_end], [th0], rtol=1e-11, atol=1e-15)
    theta_ref = sol.y[0, -1]
    coll_ref  = n_in * cbar * A_wall * dt * n_out / MPF

    rows, scatter, fail = [], None, False

    # --- coverage ----------------------------------------------------------------------
    if args.coverage and os.path.exists(args.coverage):
        with h5py.File(args.coverage, 'r') as f:
            data  = f['SurfaceData'][:]                 # (sides, nSurfSample, nSurfSample, nVar)
            names = [v.decode() for v in f.attrs['VarNamesSurface']]
        theta   = data[..., names.index('Spec001_Coverage')].ravel()   # every sub-surface
        scatter = (len(theta), 100 * theta.std(ddof=1) / theta.mean())
        rows.append(('coverage', theta.mean(), theta.std(ddof=1) / math.sqrt(len(theta)),
                     theta_ref, TOL_COV))
    else:
        print('skipping the coverage check: no DSMCSurfState file')

    # --- number of impacts -------------------------------------------------------------
    if args.surface_csv and os.path.exists(args.surface_csv):
        csv  = np.genfromtxt(args.surface_csv, delimiter=',', names=True)
        cols = csv.dtype.names
        coll = csv[cols[1]][csv[cols[0]] > 0.2 * t_end]   # skip the start-up of the run
        rows.append(('nSurfColl/output', coll.mean(), coll.std(ddof=1) / math.sqrt(len(coll)),
                     coll_ref, TOL_COLL))
    else:
        print('skipping the impact check: no SurfaceAnalyze file')

    if not rows:
        sys.exit('verify_analytic.py: none of the requested files exist')

    print('%-18s %14s %14s %10s' % (NAME, 'simulation', '0D reference', 'deviation'))
    print('-' * 60)
    for label, sim, sd, ref, tol in rows:
        dev  = sim / ref - 1.0
        flag = '' if abs(dev) <= tol else '   <== outside %.0f%%' % (100 * tol)
        if flag: fail = True
        print('%-18s %9.5g+-%-4.2g %14.5g %+9.2f%%%s' % (label, sim, sd, ref, 100 * dev, flag))
    print('-' * 60)
    if scatter:
        print('theta %.3f -> %.5f, averaged over %d catalytic sub-surfaces, per-sub-surface scatter %.1f%%'
              % (th0, theta_ref, scatter[0], scatter[1]))
    print('files: %s' % ', '.join(f for f in (args.coverage, args.surface_csv) if f and os.path.exists(f)))
    return 1 if fail else 0


if __name__ == '__main__':
    sys.exit(main())
