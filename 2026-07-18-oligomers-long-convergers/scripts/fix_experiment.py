#!/usr/bin/env python3
"""Baseline vs flat-surface-accept: step count, final energy, geometry RMSD."""
import json, os, sys
os.environ.setdefault('OMP_NUM_THREADS','4'); os.environ.setdefault('MKL_NUM_THREADS','4')
import numpy as np
from berny import Berny, geomlib
from berny.benchmarks import load_reference, require_geometries
from berny.solvers import XTBSolver

REF = load_reference('oligomers'); DATA = require_geometries('oligomers')
HK = 627.5
FLAT = 1e-6  # Ha ~ 6e-4 kcal/mol

def run(name, flat_energy, maxsteps=200):
    r = REF[name]
    geom = geomlib.readfile(str(DATA / r['file']))
    berny = Berny(geom, maxsteps=maxsteps, flat_energy=flat_energy)
    solver = XTBSolver(charge=r['charge'], mult=r['mult'])
    next(solver)
    last_geom = None; energies=[]
    for g in berny:
        e, grad = solver.send((list(g), g.lattice)); energies.append(e)
        berny.send((e, grad)); last_geom = g
    coords = np.array(last_geom.coords)
    return berny.converged, berny._n, energies[-1], coords

TARGETS = ['polyglycine_n4','polyglycine_n5','polyglycine_n6','polyglycine_n7',
           'polyserine_n3','polyserine_n5','polyserine_n6','polyserine_n7','polyserine_n8',
           'polyalanine_n6','polyalanine_n8','nylon6_n4','nylon6_n7',
           'nylon6_n5','nylon6_n6','polyglycine_n8',
           # controls that already converge fast — must be unaffected
           'polyethylene_n8','thiophene_n4','anthracene','pentacene']

print(f"{'name':18s}{'base_n':>7}{'fix_n':>7}{'saved':>7}{'dE(kcal)':>11}{'RMSD(A)':>10}{'base_conv':>10}{'fix_conv':>9}")
rows=[]
for name in TARGETS:
    bc, bn, be, bg = run(name, 0.0)
    fc, fn, fe, fg = run(name, FLAT)
    # align by centroid only (same optimizer path, so atom order identical)
    rmsd = float(np.sqrt(((bg-bg.mean(0))-(fg-fg.mean(0)))**2).sum(1).mean()**0.5) if bg.shape==fg.shape else float('nan')
    dE = (fe-be)*HK
    saved = bn-fn if (bc and fc) else None
    rows.append(dict(name=name,bn=bn,fn=fn,bc=bc,fc=fc,dE=dE,rmsd=rmsd,saved=saved))
    print(f"{name:18s}{bn:>7}{fn:>7}{str(saved):>7}{dE:>11.4f}{rmsd:>10.4f}{str(bc):>10}{str(fc):>9}")

conv_saved=[r['saved'] for r in rows if r['saved'] is not None]
print(f"\nsaved steps (both converged): total={sum(conv_saved)} on {len(conv_saved)} mols")
print(f"max |dE| final: {max(abs(r['dE']) for r in rows):.4f} kcal/mol")
print(f"max RMSD: {max(r['rmsd'] for r in rows if not np.isnan(r['rmsd'])):.4f} A")
json.dump(rows, open('./data/fix_results.json','w'), indent=2)
