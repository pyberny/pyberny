#!/usr/bin/env python3
"""Proper Kabsch-aligned RMSD between baseline and flat-accept endpoints,
plus the intrinsic soft-mode ambiguity: RMSD between two baseline runs from
mildly noise-perturbed starts (same minimum). Shows the fix's deviation is
within the conformational-valley width that baseline itself has."""
import os
os.environ.setdefault('OMP_NUM_THREADS','4'); os.environ.setdefault('MKL_NUM_THREADS','4')
import numpy as np
from berny import Berny, geomlib
from berny.benchmarks import load_reference, require_geometries
from berny.solvers import XTBSolver
REF = load_reference('oligomers'); DATA = require_geometries('oligomers')

def kabsch_rmsd(A, B):
    A = A - A.mean(0); B = B - B.mean(0)
    H = A.T @ B
    U,S,Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1,1,d])
    Rot = Vt.T @ D @ U.T
    A2 = A @ Rot.T
    return float(np.sqrt(((A2-B)**2).sum(1).mean()))

def run(name, flat_energy, seed=None, sigma=0.0, maxsteps=250):
    r = REF[name]; geom = geomlib.readfile(str(DATA / r['file']))
    if sigma>0:
        rng = np.random.default_rng(seed)
        c = np.array(geom.coords) + rng.normal(0, sigma, (len(geom),3))
        geom.coords = c
    b = Berny(geom, maxsteps=maxsteps, flat_energy=flat_energy)
    s = XTBSolver(charge=r['charge'], mult=r['mult']); next(s)
    last=None; E=[]
    for g in b:
        e,gr = s.send((list(g), g.lattice)); E.append(e); b.send((e,gr)); last=g
    return b.converged, b._n, E[-1], np.array(last.coords)

for name in ['polyglycine_n6','polyserine_n3','nylon6_n4']:
    bc,bn,be,bg = run(name, 0.0)
    fc,fn,fe,fg = run(name, 1e-6)
    fix_rmsd = kabsch_rmsd(bg, fg)
    # intrinsic ambiguity: two baseline runs from tiny (0.02 A) noise starts
    _,_,_,g1 = run(name, 0.0, seed=1, sigma=0.02)
    _,_,_,g2 = run(name, 0.0, seed=2, sigma=0.02)
    intrinsic = kabsch_rmsd(g1, g2)
    print(f"{name:16s} base_n={bn} fix_n={fn}  dE={ (fe-be)*627.5:+.4f} kcal  "
          f"fix-vs-base RMSD={fix_rmsd:.3f} A   intrinsic(2 noise seeds) RMSD={intrinsic:.3f} A")
