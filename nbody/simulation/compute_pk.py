#!/usr/bin/env python3
"""
compute_pk.py — Compute P(k) from Gadget-4 EU snapshots
NB10 — Evaporating Universe N-Body Analysis

Uses Pylians3 for power spectrum computation.
Generates P(k) data files + comparison plots for paper.

Usage:
  python3 compute_pk.py --snapshot output/snapshot_006.hdf5
  python3 compute_pk.py --all-snapshots output/
  python3 compute_pk.py --all-snapshots output/ --plot

Output:
  analysis/pk_z{redshift}.dat    — k [h/Mpc], P(k) [(Mpc/h)³]
  analysis/pk_all_redshifts.png  — Combined P(k) plot
  analysis/pk_ratio_vs_lcdm.png  — P_EU/P_LCDM ratio

Author: EU Pipeline (NB10)
Date: 2026-05-31
"""

import numpy as np
import h5py
import argparse
import os
import sys
import glob

# ============================================================
# EU Parameters (must match generate_eu_tables.py)
# ============================================================
H0 = 68.886
h = H0 / 100.0
omega_b = 0.02237
omega_cdm = 0.1193
Omega_b = omega_b / h**2
Omega_cdm_prim = omega_cdm / h**2
Omega_m = Omega_b + Omega_cdm_prim

# Box size in Mpc/h (must match param_eu.txt)
BOX_SIZE_KPCH = 500000.0  # kpc/h
BOX_SIZE_MPCH = BOX_SIZE_KPCH / 1000.0  # 500 Mpc/h

# Grid resolution for P(k) estimation
MESH_PK = 512  # Nyquist is k_Ny = pi * N_mesh / L


def read_snapshot(filename):
    """Read particle positions and header from Gadget-4 HDF5 snapshot."""
    print(f"  Reading {filename}...")
    
    with h5py.File(filename, 'r') as f:
        header = dict(f['Header'].attrs)
        z = header['Redshift']
        a = header['Time']
        boxsize = header['BoxSize']  # kpc/h
        n_part = header['NumPart_Total'][1]  # DM particles
        
        # Read DM positions (Type 1)
        pos = f['PartType1/Coordinates'][:]  # shape (N, 3), kpc/h
        
        # Read mass if available
        if 'Masses' in f['PartType1']:
            masses = f['PartType1/Masses'][:]
        else:
            mass_table = header.get('MassTable', [0]*6)
            masses = np.full(len(pos), mass_table[1])
    
    print(f"  z = {z:.4f}, a = {a:.4f}")
    print(f"  N_particles = {n_part:,}")
    print(f"  Box = {boxsize:.0f} kpc/h = {boxsize/1000:.0f} Mpc/h")
    
    return pos, masses, z, a, boxsize


def compute_pk_pylians(pos, boxsize_mpch, mesh=512):
    """Compute P(k) using Pylians3 MAS + FFT."""
    try:
        import MAS_library as MASL
        import Pk_library as PKL
    except ImportError:
        print("[ERROR] Pylians3 not installed!")
        print("        Install with: pip3 install Pylians3")
        print("        Falling back to simple FFT method...")
        return compute_pk_simple(pos, boxsize_mpch, mesh)
    
    # Convert positions to Mpc/h (Pylians expects Mpc/h)
    pos_mpch = pos / 1000.0  # kpc/h -> Mpc/h
    pos_mpch = pos_mpch.astype(np.float32)
    
    # Create density field using CIC (Cloud-in-Cell) assignment
    delta = np.zeros((mesh, mesh, mesh), dtype=np.float32)
    MASL.MA(pos_mpch, delta, boxsize_mpch, 'CIC', verbose=False)
    
    # Convert to overdensity: delta = rho/rho_mean - 1
    delta /= np.mean(delta)
    delta -= 1.0
    
    # Compute P(k) with Pylians
    Pk = PKL.Pk(delta, boxsize_mpch, axis=0, MAS='CIC', threads=4)
    
    k = Pk.k3D       # k in h/Mpc
    pk = Pk.Pk[:, 0]  # P(k) in (Mpc/h)^3 (monopole)
    nmodes = Pk.Nmodes3D
    
    return k, pk, nmodes


def compute_pk_simple(pos, boxsize_mpch, mesh=512):
    """Simple FFT-based P(k) computation (fallback if Pylians not available)."""
    print("  Using simple FFT P(k) (no deconvolution)...")
    
    pos_mpch = pos / 1000.0
    
    # CIC assignment to grid
    delta = np.zeros((mesh, mesh, mesh), dtype=np.float32)
    cell_size = boxsize_mpch / mesh
    
    # Simple NGP (Nearest Grid Point) for speed
    idx = ((pos_mpch / boxsize_mpch * mesh) % mesh).astype(int)
    for i in range(len(pos_mpch)):
        delta[idx[i, 0], idx[i, 1], idx[i, 2]] += 1.0
    
    # Overdensity
    delta = delta / np.mean(delta) - 1.0
    
    # FFT
    delta_k = np.fft.rfftn(delta)
    pk_3d = np.abs(delta_k)**2 * (boxsize_mpch / mesh)**3 / boxsize_mpch**3
    
    # Spherical averaging
    k_fund = 2.0 * np.pi / boxsize_mpch
    kx = np.fft.fftfreq(mesh, d=1.0/mesh) * k_fund
    ky = np.fft.fftfreq(mesh, d=1.0/mesh) * k_fund
    kz = np.fft.rfftfreq(mesh, d=1.0/mesh) * k_fund
    kgrid = np.sqrt(kx[:, None, None]**2 + ky[None, :, None]**2 + kz[None, None, :]**2)
    
    # Bin in k shells
    k_bins = np.arange(k_fund, mesh//2 * k_fund, k_fund)
    k_centers = 0.5 * (k_bins[:-1] + k_bins[1:])
    pk_binned = np.zeros(len(k_centers))
    nmodes = np.zeros(len(k_centers), dtype=int)
    
    for i in range(len(k_centers)):
        mask = (kgrid >= k_bins[i]) & (kgrid < k_bins[i+1])
        if np.any(mask):
            pk_binned[i] = np.mean(pk_3d[mask])
            nmodes[i] = np.sum(mask)
    
    # Remove empty bins
    valid = nmodes > 0
    return k_centers[valid], pk_binned[valid], nmodes[valid]


def compute_pk_lcdm_halofit(k_array, z):
    """Compute LCDM P(k) using simple fitting formula (Eisenstein & Hu 1998).
    
    This is an approximation for the ratio plot.
    For precise comparison, use CLASS or CAMB output.
    """
    # Transfer function (Eisenstein & Hu 1998, zero-baryon approximation)
    Omega_m_z0 = Omega_m
    Gamma = Omega_m_z0 * h  # shape parameter
    
    q = k_array / (Gamma * h)  # Note: k is in h/Mpc
    T_k = np.log(1 + 2.34 * q) / (2.34 * q) * \
          (1 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4)**(-0.25)
    
    # Linear P(k) ~ k^n_s * T(k)^2
    n_s = 0.9649
    pk_lcdm = k_array**n_s * T_k**2
    
    # Normalize to sigma_8 = 0.8275 at z=0 (approximate)
    # This is rough — for the paper, use CLASS output
    sigma8_target = 0.8275
    
    # Growth factor D(z) for LCDM
    Omega_L = 1.0 - Omega_m_z0
    az = 1.0 / (1.0 + z)
    Omega_m_z = Omega_m_z0 / (Omega_m_z0 + Omega_L * az**3)
    D_z = az * (5.0/2.0 * Omega_m_z) / \
          (Omega_m_z**(4.0/7.0) - Omega_L + (1 + Omega_m_z/2.0) * (1 + Omega_L/70.0))
    D_0 = (5.0/2.0 * Omega_m_z0) / \
          (Omega_m_z0**(4.0/7.0) - Omega_L + (1 + Omega_m_z0/2.0) * (1 + Omega_L/70.0))
    growth = D_z / D_0
    
    pk_lcdm *= growth**2
    
    return pk_lcdm


def save_pk(k, pk, nmodes, z, outdir):
    """Save P(k) to text file."""
    outfile = os.path.join(outdir, f"pk_z{z:.2f}.dat")
    header = (f"# P(k) from EU N-Body Simulation (NB10)\n"
              f"# Redshift: z = {z:.4f}\n"
              f"# Columns: k [h/Mpc], P(k) [(Mpc/h)^3], N_modes\n"
              f"# Box: {BOX_SIZE_MPCH:.0f} Mpc/h, Mesh: {MESH_PK}\n")
    
    data = np.column_stack([k, pk, nmodes])
    np.savetxt(outfile, data, header=header, fmt='%.6e  %.6e  %d')
    print(f"  Saved: {outfile}")
    return outfile


def plot_pk_all(pk_data, outdir):
    """Plot P(k) at all redshifts."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARNING] matplotlib not available, skipping plots.")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(pk_data)))
    
    # Left panel: P(k) at each redshift
    for i, (z, k, pk) in enumerate(sorted(pk_data, key=lambda x: -x[0])):
        ax1.loglog(k, pk, color=colors[i], label=f'z = {z:.1f}', alpha=0.8, lw=1.5)
    
    ax1.set_xlabel(r'$k$ [$h$/Mpc]', fontsize=14)
    ax1.set_ylabel(r'$P(k)$ [(Mpc/$h$)$^3$]', fontsize=14)
    ax1.set_title('EU Power Spectrum (N-Body)', fontsize=15)
    ax1.legend(fontsize=11)
    ax1.set_xlim(0.01, 10)
    ax1.grid(True, alpha=0.3)
    
    # Right panel: Ratio P_EU / P_LCDM at z=0
    z0_data = [d for d in pk_data if d[0] < 0.01]
    if z0_data:
        z, k, pk = z0_data[0]
        pk_lcdm = compute_pk_lcdm_halofit(k, 0.0)
        # Normalize both to match at k=0.01
        idx_norm = np.argmin(np.abs(k - 0.01))
        if idx_norm > 0 and pk_lcdm[idx_norm] > 0:
            pk_lcdm *= pk[idx_norm] / pk_lcdm[idx_norm]
            ratio = pk / pk_lcdm
            
            ax2.semilogx(k, ratio, 'b-', lw=2, label='EU / LCDM (N-body)')
            ax2.axhline(y=1.0, color='gray', ls='--', alpha=0.5)
            ax2.fill_between(k, 0.95, 1.0, alpha=0.1, color='red',
                           label=r'EU suppression zone')
            ax2.set_xlabel(r'$k$ [$h$/Mpc]', fontsize=14)
            ax2.set_ylabel(r'$P_{\rm EU}(k) / P_{\Lambda\rm CDM}(k)$', fontsize=14)
            ax2.set_title('Power Spectrum Ratio (z=0)', fontsize=15)
            ax2.set_xlim(0.01, 5)
            ax2.set_ylim(0.9, 1.05)
            ax2.legend(fontsize=11)
            ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    outfile = os.path.join(outdir, 'pk_all_redshifts.png')
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"  Plot saved: {outfile}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Compute P(k) from EU N-Body snapshots")
    parser.add_argument("--snapshot", type=str, help="Single snapshot file")
    parser.add_argument("--all-snapshots", type=str, help="Directory with all snapshots")
    parser.add_argument("--mesh", type=int, default=MESH_PK, help=f"Mesh size for P(k) (default: {MESH_PK})")
    parser.add_argument("--plot", action="store_true", help="Generate plots")
    parser.add_argument("--outdir", type=str, default="analysis", help="Output directory")
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    # Collect snapshots
    snapshots = []
    if args.snapshot:
        snapshots = [args.snapshot]
    elif args.all_snapshots:
        snapshots = sorted(glob.glob(os.path.join(args.all_snapshots, "snapshot_*.hdf5")))
        if not snapshots:
            print(f"[ERROR] No snapshots found in {args.all_snapshots}")
            sys.exit(1)
    else:
        parser.print_help()
        sys.exit(1)
    
    print(f"Found {len(snapshots)} snapshot(s)")
    print(f"Mesh size: {args.mesh}")
    print(f"Output: {args.outdir}/")
    print()
    
    # Process each snapshot
    pk_data = []  # (z, k, pk) for plotting
    
    for snap_file in snapshots:
        print(f"--- Processing {os.path.basename(snap_file)} ---")
        
        # Read snapshot
        pos, masses, z, a, boxsize = read_snapshot(snap_file)
        
        # Compute P(k)
        print(f"  Computing P(k) (mesh={args.mesh})...")
        k, pk, nmodes = compute_pk_pylians(pos, BOX_SIZE_MPCH, mesh=args.mesh)
        
        # Save
        save_pk(k, pk, nmodes, z, args.outdir)
        pk_data.append((z, k, pk))
        
        # Print summary
        k_nl = k[np.argmax(k * pk)]  # non-linear scale
        print(f"  k range: [{k[0]:.4f}, {k[-1]:.2f}] h/Mpc")
        print(f"  P(k=0.1): {pk[np.argmin(np.abs(k-0.1))]:.2e} (Mpc/h)^3")
        print(f"  k_NL ~ {k_nl:.2f} h/Mpc")
        print()
    
    # Generate plots
    if args.plot and pk_data:
        print("--- Generating plots ---")
        plot_pk_all(pk_data, args.outdir)
    
    print("============================================")
    print(f" Analysis complete! {len(snapshots)} snapshot(s) processed.")
    print(f" Results in: {args.outdir}/")
    print("============================================")


if __name__ == "__main__":
    main()
