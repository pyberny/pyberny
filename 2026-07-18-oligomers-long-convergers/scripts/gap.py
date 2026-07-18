#!/usr/bin/env python3
"""Quantify how many tail steps are spent with forces already converged but
displacement criteria still tripping — i.e. the savings a force-based
early-accept (Gaussian-style) would yield."""
import json, os
import numpy as np
from berny.berny import BernyParams

TR = os.environ.get('OLIGO_TRACES','./traces')
REF = json.load(open('/home/user/pyberny/src/berny/benchmarks/oligomers/reference.json'))
P = BernyParams()
HK = 627.5

def load(name):
    return json.load(open(f'{TR}/xtb-oligomers-{name}.trace.json'))

def crit(r, nm):
    for c in r['convergence']['criteria']:
        if c['name'] == nm:
            return c
    return None

def step_forces_ok(r, factor=1.0):
    grms = crit(r, 'Gradient RMS'); gmax = crit(r, 'Gradient maximum')
    return (grms['value'] < factor*P.gradientrms) and (gmax['value'] < factor*P.gradientmax)

def analyze(name):
    t = load(name); n = len(t)
    conv = t[-1]['converged']
    # first step from which BOTH force criteria stay satisfied to the end
    force_ok_from = None
    for i in range(n):
        if all(step_forces_ok(t[j]) for j in range(i, n)):
            force_ok_from = i + 1
            break
    # first step where forces are 1/3 below thresholds (Gaussian "well below") and stay
    force_loose_from = None
    for i in range(n):
        if all(step_forces_ok(t[j], factor=1.0) for j in range(i, n)):
            force_loose_from = i + 1
            break
    # energies
    E = [r['energy'] for r in t]; Emin = min(E)
    # max energy change over the "wasted tail" (from force_ok_from to end)
    dE_tail = None
    if force_ok_from:
        seg = E[force_ok_from-1:]
        dE_tail = (max(seg) - min(seg)) * HK
    return dict(name=name, family=REF[name]['family'], atoms=REF[name]['atoms'],
               n=n, conv=conv, force_ok_from=force_ok_from,
               saved=(n - force_ok_from) if force_ok_from else None,
               dE_tail_kcal=dE_tail)

names = sorted(REF)
rows = [analyze(x) for x in names]

print("Molecules where forces converge >=8 steps before the run ends:")
print(f"{'name':20s}{'fam':10s}{'n':>4}{'conv':>6}{'force_ok@':>10}{'saved':>7}{'dE_tail(kcal)':>14}")
big = [r for r in rows if r['saved'] and r['saved'] >= 8]
for r in sorted(big, key=lambda x: -x['saved']):
    print(f"{r['name']:20s}{r['family']:10s}{r['n']:4d}{str(r['conv']):>6}"
          f"{str(r['force_ok_from']):>10}{r['saved']:7d}{r['dE_tail_kcal']:14.4f}")

# aggregate
saved_all = [r['saved'] for r in rows if r['saved'] is not None]
print(f"\nTotal measured steps (converged+ceiling): {sum(r['n'] for r in rows)}")
print(f"Total steps AFTER forces converged (potential early-accept savings): {sum(saved_all)}")
tail_dE = [r['dE_tail_kcal'] for r in rows if r['dE_tail_kcal'] is not None]
print(f"Max energy change anywhere in a force-converged tail: {max(tail_dE):.4f} kcal/mol")
print(f"Median: {np.median(tail_dE):.5f} kcal/mol")

# what about the ceiling non-convergers specifically
print("\n=== Ceiling non-convergers: were forces already converged at step 100? ===")
for name in ['nylon6_n5','nylon6_n6','nylon6_n8','polyglycine_n8']:
    t = load(name); last = t[-1]
    grms = crit(last,'Gradient RMS'); gmax = crit(last,'Gradient maximum')
    srms = crit(last,'Step RMS'); smax = crit(last,'Step maximum')
    def fmt(c):
        return f"{c['value']:.2e}({'OK' if c['matched'] else 'NO'})" if c else "on-sphere"
    print(f"{name:16s} final: gRMS={grms['value']:.2e}/{P.gradientrms:.1e}({'OK' if grms['matched'] else 'NO'}) "
          f"gMax={gmax['value']:.2e}/{P.gradientmax:.1e}({'OK' if gmax['matched'] else 'NO'}) "
          f"sRMS={fmt(srms)} sMax={fmt(smax)}")
    r = [x for x in rows if x['name']==name][0]
    print(f"                 force_ok_from_step={r['force_ok_from']}  (both force criteria met from here to step 100)")
