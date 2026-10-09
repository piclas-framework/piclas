import os
import re
import itertools
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =============================================================================== #
# Verify the cross-section based collision/relaxation/chemistry models in PICLas
# -> Plot the input cross-sections sigma(E) from the LXCat database and overlay
#    the cross-sections recovered from the PICLas reservoir results.
#
# Each PartAnalyze column corresponds to a collision pair (i,j); one partner is the
# background gas (~at rest), the other is the projectile with a velocity sweep. The
# relative speed g is the projectile velocity, the collision energy is
# E = 0.5 * m_red(i,j) * g^2. From the per-pair rate the cross-section is recovered:
#   * Rates [1/s] (collision rate, backscatter, electronic relaxation):
#       sigma = Rate / (n_i * n_j * V * g)
#   * Reaction rate coefficient k = sigma*g [m3/s] (e.g. ionization):
#       sigma = k / g
#
# Everything (database name, masses, densities, velocities, mesh volume) is read
# from parameter.ini / hopr.ini; the processes to compare are discovered from the
# PartAnalyze CSV header and mapped to the database by species name / threshold.
# =============================================================================== #

# --- Setup -------------------------------------------------------------------- #
PARAMETER_FILE = 'parameter.ini'
HOPR_FILE = 'hopr.ini'
E_CHG  = 1.602176634E-19   # eV -> J conversion

_FLOAT = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'   # int/float incl. scientific notation


# --- parameter.ini parsing ---------------------------------------------------- #
def read_parameter(key):
    """Return the right-hand side of 'key = value' from parameter.ini (comments stripped)."""
    with open(PARAMETER_FILE) as f:
        for line in f:
            lhs, sep, rhs = line.partition('=')
            if sep and lhs.strip().lower() == key.lower():
                return rhs.split('!')[0].strip()
    raise KeyError(f'{key} not found in {PARAMETER_FILE}')


def read_parameter_list(key):
    """Parse a comma-separated numeric list parameter."""
    return [float(v) for v in read_parameter(key).split(',')]


def read_index_vector(key):
    """Parse a PICLas index vector like '(/3,1,1,0/)' into a list of ints (zeros dropped)."""
    return [int(n) for n in re.findall(r'-?\d+', read_parameter(key)) if int(n) != 0]


def mesh_volume_from_hopr():
    """Bounding-box volume [m3] of the (Cartesian) mesh, from the Corner vector in hopr.ini."""
    with open(HOPR_FILE) as f:
        match = re.search(r'Corner\s*=\s*\(/(.*?)/\)', f.read(), re.IGNORECASE | re.DOTALL)
    if not match:
        raise KeyError(f'Corner not found in {HOPR_FILE}')
    pts = np.array([float(x) for x in re.findall(_FLOAT, match.group(1))]).reshape(-1, 3)
    return float(np.prod(pts.max(axis=0) - pts.min(axis=0)))


def read_species():
    """Read per-species name, mass, density and (sorted) velocity sweep from parameter.ini."""
    species = {}
    for i in range(1, int(read_parameter('Part-nSpecies')) + 1):
        try:
            velo = sorted(read_parameter_list(f'Part-Species{i}-Init1-VeloIC'))
        except KeyError:
            velo = [0.0]
        species[i] = {
            'name': read_parameter(f'Part-Species{i}-SpeciesName'),
            'mass': float(read_parameter(f'Part-Species{i}-MassIC')),
            'dens': float(read_parameter(f'Part-Species{i}-Init1-PartDensity')),
            'velo': velo,
        }
    return species


# --- Database lookup ---------------------------------------------------------- #
def find_group(hfile, species_indices):
    """Find the top-level collision-pair group matching the given species indices."""
    target = sorted(SPECIES[i]['name'].lower() for i in species_indices)
    for grp in hfile:
        if sorted(tok.lower() for tok in grp.split('-')) == target:
            return grp
    raise KeyError(f'No database group for species {species_indices} (names {target})')


def find_elec_dataset(hfile, grp, level):
    """Find the ELECTRONIC dataset in group whose threshold matches the level [eV]."""
    for name, ds in hfile[f'{grp}/ELECTRONIC'].items():
        threshold = float(ds.attrs.get('Threshold [eV]', name))
        if abs(threshold - level) <= 1.0E-6 * max(1.0, level):
            return f'{grp}/ELECTRONIC/{name}'
    raise KeyError(f'No ELECTRONIC dataset at {level} eV in group {grp}')


def find_reaction_dataset(hfile, grp, product_names):
    """Find the REACTION dataset in group whose product species match (order-independent)."""
    target = sorted(n.lower() for n in product_names)
    for name in hfile[f'{grp}/REACTION']:
        if sorted(tok.lower() for tok in name.split('-')) == target:
            return f'{grp}/REACTION/{name}'
    raise KeyError(f'No REACTION dataset for products {product_names} in group {grp}')


def require_dataset(hfile, path):
    if path not in hfile:
        raise KeyError(f'No dataset {path}')
    return path


def effective_xsec(hfile, grp):
    """Effective collision cross-section of a pair: ELASTIC plus every sub-process
    (electronic, reaction, backscatter, ...) summed onto the ELASTIC energy grid.
    This mirrors how PICLas builds the collision cross-section when the database
    provides ELASTIC rather than an EFFECTIVE dataset, and is what the collision
    rate (CollRate) actually samples."""
    elastic = np.array(hfile[require_dataset(hfile, f'{grp}/ELASTIC')])
    energy, total = elastic[:, 0], elastic[:, 1].copy()

    def add(name, obj):
        if isinstance(obj, h5py.Dataset) and name != 'ELASTIC':
            data = np.array(obj)
            total[:] += np.interp(energy, data[:, 0], data[:, 1])

    hfile[grp].visititems(add)
    return np.column_stack([energy, total])


def projectile_velocity(pair, n_runs):
    """Relative speed array: the velocity sweep of the projectile (non-background) species."""
    cand = [s for s in pair if len(SPECIES[s]['velo']) == n_runs]
    if not cand:
        raise KeyError(f'No projectile velocity sweep of length {n_runs} for pair {pair}')
    return np.array(SPECIES[cand[0]]['velo'])


def make_pair_process(pair, n_runs, kind, xsec, label, color):
    """Attach per-run relative speed, energy and density prefactor to a process."""
    i, j = pair
    g = projectile_velocity(pair, n_runs)
    m_red = SPECIES[i]['mass'] * SPECIES[j]['mass'] / (SPECIES[i]['mass'] + SPECIES[j]['mass'])
    prefactor = (0.5 if i == j else 1.0) * SPECIES[i]['dens'] * SPECIES[j]['dens']
    return {'kind': kind, 'label': label, 'color': color, 'xsec': xsec,
            'g': g, 'energy': 0.5 * m_red * g**2 / E_CHG, 'prefactor': prefactor}


def discover_processes(hfile, header_columns, n_runs):
    """Build the list of processes to compare from the PartAnalyze CSV header.

    Each column is mapped to a collision pair and a database dataset. Columns whose
    pair has no matching database group/dataset (e.g. e-e, total collision rate) are
    skipped. 'kind' is 'reaction' (rate coefficient) or 'rate' (raw rate [1/s]).
    """
    colors = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    processes = {}
    for col in header_columns:
        m_back = re.search(r'BackscatterCollRate(\d{3})\+(\d{3})', col)
        m_coll = re.search(r'CollRate(\d{3})\+(\d{3})', col)
        m_elec = re.search(r'ElecRelaxRate(\d{3})\+(\d{3})-([\d.]+)', col)
        m_reac = re.search(r'Reaction(\d+)', col)
        try:
            if m_back:
                pair = (int(m_back.group(1)), int(m_back.group(2)))
                grp = find_group(hfile, pair)
                xsec = np.array(hfile[require_dataset(hfile, f'{grp}/BACKSCATTER')])
                proc = make_pair_process(pair, n_runs, 'rate', xsec, f'{grp} backscatter', next(colors))
            elif m_coll:
                pair = (int(m_coll.group(1)), int(m_coll.group(2)))
                grp = find_group(hfile, pair)
                # CollRate samples the effective cross-section (ELASTIC + all sub-processes)
                xsec = effective_xsec(hfile, grp)
                proc = make_pair_process(pair, n_runs, 'rate', xsec, f'{grp} effective', next(colors))
            elif m_elec:
                pair = (int(m_elec.group(1)), int(m_elec.group(2)))
                grp = find_group(hfile, pair)
                xsec = np.array(hfile[find_elec_dataset(hfile, grp, float(m_elec.group(3)))])
                proc = make_pair_process(pair, n_runs, 'rate', xsec,
                                         f'Excitation {float(m_elec.group(3)):g} eV', next(colors))
            elif m_reac:
                idx = int(m_reac.group(1))
                reactants = read_index_vector(f'DSMC-Reaction{idx}-Reactants')
                products = read_index_vector(f'DSMC-Reaction{idx}-Products')
                grp = find_group(hfile, reactants)
                product_names = [SPECIES[k]['name'] for k in products]
                label = ' + '.join(SPECIES[k]['name'] for k in reactants) + ' → ' + ' + '.join(product_names)
                xsec = np.array(hfile[find_reaction_dataset(hfile, grp, product_names)])
                proc = make_pair_process(tuple(reactants[:2]), n_runs, 'reaction', xsec, label, next(colors))
            else:
                continue  # TIME, TotalCollRate or any non-process column
        except KeyError:
            continue       # pair without database cross-section -> nothing to compare
        processes[col] = proc
    return processes


def rate_to_xsec(proc, rate):
    """Recover the cross-section [m2] from the PICLas rate output for one run."""
    g = proc['g']
    if proc['kind'] == 'reaction':
        return rate / g                              # rate coefficient k = sigma*g
    return rate / (proc['prefactor'] * VOLUME * g)   # raw rate [1/s]


# --- Read setup from the parameter file --------------------------------------- #
DATABASE = read_parameter('Particles-CollXSec-Database')
# VOLUME = mesh_volume_from_hopr()
VOLUME = 1.57079632679490E-013
SPECIES = read_species()

# --- Reference files (sorted by ascending energy code = ascending velocity) ---- #
directory = os.getcwd()
ref_files = sorted(f for f in os.listdir(directory) if f.startswith('PartAnalyze'))
n_runs = len(ref_files)

# --- Discover processes from the CSV header and recover the cross-sections ----- #
with h5py.File(DATABASE, 'r') as hfile:
    header = list(pd.read_csv(ref_files[0], nrows=0).columns)
    processes = discover_processes(hfile, header, n_runs)

# sigma per process: array over runs (recovered all at once), energies stored per process
recovered = {col: rate_to_xsec(p, np.array([pd.read_csv(f)[col].tail(1).values[0]
                                            for f in ref_files]))
             for col, p in processes.items()}

# --- Plot --------------------------------------------------------------------- #
fig, ax = plt.subplots(figsize=(16, 6))
xmax = max(p['energy'].max() for p in processes.values()) * 1.1

for col, p in processes.items():
    ref = p['xsec']
    mask = ref[:, 0] <= xmax                      # restrict to simulated range
    ax.plot(ref[mask, 0], ref[mask, 1], '-', color=p['color'], label=f"{p['label']} (database)")
    ax.plot(p['energy'], recovered[col], 'o', color=p['color'], markeredgecolor='k',
            label=f"{p['label']} (PICLas)")

ax.set_xlabel('Collision energy [eV]')
ax.set_ylabel(r'Cross-section [m$^2$]')
ax.set_title('PICLas cross-section verification')
ax.set_xlim(left=10, right=xmax)
ax.set_yscale('log')
ax.set_ylim(bottom=1e-23)
ax.grid(True, which='both', ls=':', alpha=0.6)
ax.legend(fontsize=7, loc='center left', bbox_to_anchor=(1.02, 0.5),ncol=2)
fig.tight_layout()
fig.savefig('Figure_CrossSection_Verification.jpg', dpi=150)
print('Wrote Figure_CrossSection_Verification.jpg')
