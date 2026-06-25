#!/usr/bin/env python3
"""Extract P(k) from Gadget4 V3 snapshots (Spain session).

Adapted from extract_pk_frankfurt.py for:
  - Spain instance paths (/mnt/data/nbody_eu/production/output)
  - Snapshots 009-014 (z=2.5 to z=0)
  - Also re-extracts 000-008 for consistency with same pipeline

Usage (on AWS instance):
  python3 /tmp/extract_pk_spain.py
  python3 /tmp/extract_pk_spain.py --snaps 9 10 11 12 13 14
"""
import numpy as np
import h5py
import os, glob, json, time, sys

# Try Pylians first, fall back to manual FFT
try:
    import MAS_library as MASL
    import Pk_library as PKL
    HAS_PYLIANS = True
    print("Pylians3 loaded OK")
except ImportError:
    HAS_PYLIANS = False
    print("WARNING: Pylians3 not available, using manual FFT")

BoxSize = 500.0   # Mpc/h
grid    = 512     # grid for P(k)
MAS     = 'CIC'
threads = 16
IC_MASS = 0.9645418902084929  # initial CDM mass (snap_000)

base   = '/mnt/data/nbody_eu/production/output'
outdir = '/tmp/pk_results_v3'
os.makedirs(outdir, exist_ok=True)


def compute_pk_manual(pos_mpch, BoxSize, grid):
    """Manual P(k) via FFT with CIC assignment + deconvolution."""
    # CIC assignment
    delta = np.zeros((grid, grid, grid), dtype=np.float32)
    H = BoxSize / grid  # cell size

    # Use numpy histogram for fast binning (NGP approximation)
    idx = np.floor(pos_mpch / H).astype(np.int32) % grid
    # Simple NGP for robustness
    np.add.at(delta, (idx[:, 0], idx[:, 1], idx[:, 2]), 1.0)

    delta = delta / np.mean(delta) - 1.0

    # FFT
    delta_k = np.fft.rfftn(delta)
    pk_3d = (np.abs(delta_k)**2) * (BoxSize**3 / grid**6)
    del delta, delta_k

    # k-grid
    k_fund = 2 * np.pi / BoxSize
    kx = np.fft.fftfreq(grid, d=1.0/grid) * k_fund
    ky = np.fft.fftfreq(grid, d=1.0/grid) * k_fund
    kz = np.fft.rfftfreq(grid, d=1.0/grid) * k_fund
    kgrid = np.sqrt(kx[:, None, None]**2 + ky[None, :, None]**2 + kz[None, None, :]**2)

    # Spherical binning
    dk = k_fund
    k_max = grid // 2 * k_fund
    k_edges = np.arange(dk/2, k_max + dk, dk)
    n_bins = len(k_edges) - 1
    k_out = np.zeros(n_bins)
    pk_out = np.zeros(n_bins)
    nmodes = np.zeros(n_bins, dtype=np.int64)

    kgrid_flat = kgrid.ravel()
    pk_flat = pk_3d.ravel()
    bin_idx = np.digitize(kgrid_flat, k_edges) - 1

    for i in range(n_bins):
        mask = bin_idx == i
        n = np.sum(mask)
        if n > 0:
            k_out[i] = np.mean(kgrid_flat[mask])
            pk_out[i] = np.mean(pk_flat[mask])
            nmodes[i] = n

    valid = nmodes > 0
    return k_out[valid], pk_out[valid], nmodes[valid]


def extract_snapshot(snapdir, snap_num):
    """Extract P(k) from one snapshot directory."""
    t0 = time.time()

    # Find snapshot files
    files = sorted(glob.glob(os.path.join(snapdir, 'snapshot_*.hdf5')))
    if not files:
        files = sorted(glob.glob(os.path.join(snapdir, 'snap_*.hdf5')))
    if not files:
        print(f"  snap_{snap_num:03d}: NO FILES FOUND")
        return None

    # Read header
    with h5py.File(files[0], 'r') as f:
        z = float(f['Header'].attrs['Redshift'])
        a = float(f['Header'].attrs['Time'])
        npart_total = int(f['Header'].attrs['NumPart_Total'][1])
        mass = float(f['Header'].attrs['MassTable'][1])

    drain_pct = (1 - mass / IC_MASS) * 100
    print(f"\nSnap {snap_num:03d}: z={z:.4f}, a={a:.6f}, mass={mass:.10f}, drain={drain_pct:.4f}%")
    print(f"  Files: {len(files)}, particles: {npart_total:,}")

    # Read all positions
    pos_list = []
    for ff in files:
        with h5py.File(ff, 'r') as f:
            if 'PartType1' in f and 'Coordinates' in f['PartType1']:
                coords = f['PartType1']['Coordinates'][:]
                # Convert kpc/h -> Mpc/h
                pos_list.append(coords / 1000.0)
    pos = np.concatenate(pos_list, axis=0).astype(np.float32)
    del pos_list
    print(f"  Loaded {len(pos):,} particles ({pos.nbytes/1e9:.1f} GB)")

    # Compute P(k)
    if HAS_PYLIANS:
        delta = np.zeros((grid, grid, grid), dtype=np.float32)
        MASL.MA(pos, delta, BoxSize, MAS, verbose=False)
        del pos
        delta /= np.mean(delta)
        delta -= 1.0
        Pk_obj = PKL.Pk(delta, BoxSize, axis=0, MAS=MAS, threads=threads)
        k = Pk_obj.k3D
        pk = Pk_obj.Pk[:, 0]
        nmodes = Pk_obj.Nmodes3D
        del delta
    else:
        k, pk, nmodes = compute_pk_manual(pos, BoxSize, grid)
        del pos

    # Save
    outfile = os.path.join(outdir, 'pk_%03d_z%.2f.txt' % (snap_num, z))
    header_txt = 'z=%.4f a=%.6f mass=%.10f\nk [h/Mpc]  P(k) [(Mpc/h)^3]  Nmodes' % (z, a, mass)
    np.savetxt(outfile, np.column_stack([k, pk, nmodes]),
               header=header_txt, fmt='%.8e  %.8e  %d')

    elapsed = time.time() - t0
    print(f"  Saved: {outfile} ({elapsed:.0f}s)")

    return {
        'z': z, 'a': a, 'mass': mass,
        'pk_file': os.path.basename(outfile),
        'n_particles': npart_total,
        'drain_pct': drain_pct
    }


if __name__ == '__main__':
    # Parse which snapshots to extract
    if '--snaps' in sys.argv:
        idx = sys.argv.index('--snaps')
        snap_nums = [int(x) for x in sys.argv[idx+1:]]
    else:
        snap_nums = list(range(15))  # all

    print(f"Extracting P(k) for snapshots: {snap_nums}")
    print(f"Base: {base}")
    print(f"Output: {outdir}")
    print(f"Grid: {grid}, MAS: {MAS}, threads: {threads}")
    print("=" * 60)

    results = {}
    for sn in snap_nums:
        snapdir = os.path.join(base, 'snapdir_%03d' % sn)
        if not os.path.isdir(snapdir):
            print(f"  snapdir_{sn:03d} not found, skipping")
            continue
        info = extract_snapshot(snapdir, sn)
        if info:
            results['snap_%03d' % sn] = info

    # Save summary JSON
    json_path = os.path.join(outdir, 'pk_results_v3.json')
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\n{'=' * 60}")
    print(f"DONE: {len(results)} P(k) extracted to {outdir}")
    print(f"JSON: {json_path}")
