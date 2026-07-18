#!/usr/bin/env python3
import json, os, glob, os, re
from collections import defaultdict
import numpy as np

TR = os.environ.get('OLIGO_TRACES','./traces')
REF = json.load(open('/home/user/pyberny/src/berny/benchmarks/oligomers/reference.json'))

def family(name):
    return REF[name]['family']

def load(name):
    return json.load(open(f'{TR}/xtb-oligomers-{name}.trace.json'))

HARTREE_KCAL = 627.5

def analyze(name):
    t = load(name)
    n = len(t)
    conv = t[-1]['converged']
    energies = [r['energy'] for r in t]
    Emin = min(energies)
    # crit trajectories
    def crit(r, nm):
        for c in r['convergence']['criteria']:
            if c['name'] == nm:
                return c
        return None
    # which criterion is the bottleneck: for the last step where NOT converged-ish,
    # find how many steps each criterion was already satisfied before the end
    names = ['Gradient RMS', 'Gradient maximum', 'Step RMS', 'Step maximum']
    # step at which each criterion last became (and stayed) satisfied
    last_unsat = {nm: 0 for nm in names}
    for i, r in enumerate(t):
        for nm in names:
            c = crit(r, nm)
            if c is not None and not c['matched']:
                last_unsat[nm] = i + 1  # 1-based step
    # trust radius trajectory
    trusts = [r.get('quadratic_step', {}).get('trust_radius') for r in t]
    trusts = [x for x in trusts if x is not None]
    # negative eigenvalues
    negs = [r.get('quadratic_step', {}).get('n_negative_eigenvalues') for r in t]
    negs = [x for x in negs if x is not None]
    lowev = [r.get('quadratic_step', {}).get('lowest_eigenvalue') for r in t]
    lowev = [x for x in lowev if x is not None]
    # step types
    on_sphere = sum(1 for r in t if r.get('quadratic_step', {}).get('on_sphere'))
    # linear search returning to best (rejecting new point)
    ls_best = sum(1 for r in t if r.get('linear_search', {}).get('method') == 'none-best')
    ls_t = [r.get('linear_search', {}).get('t') for r in t if r.get('linear_search')]
    # energy convergence: step at which E within 0.1 kcal/mol of final min
    e_settle = None
    for i, e in enumerate(energies):
        if (e - Emin) * HARTREE_KCAL < 0.1:
            e_settle = i + 1
            break
    # gradient rms trajectory
    grms = [crit(r, 'Gradient RMS')['value'] for r in t]
    # fletcher
    fl = [r.get('trust_update', {}).get('fletcher') for r in t]
    fl = [x for x in fl if x is not None]
    return {
        'name': name, 'family': family(name), 'atoms': REF[name]['atoms'],
        'n': n, 'conv': conv,
        'bottleneck': max(last_unsat, key=last_unsat.get),
        'last_unsat': last_unsat,
        'e_settle': e_settle, 'e_settle_frac': (e_settle / n) if e_settle else None,
        'tail_frac': round(1 - (e_settle / n), 2) if e_settle else None,
        'trust_min': min(trusts) if trusts else None,
        'trust_final': trusts[-1] if trusts else None,
        'neg_max': max(negs) if negs else None,
        'neg_nonzero_steps': sum(1 for x in negs if x > 0),
        'lowev_min': min(lowev) if lowev else None,
        'on_sphere': on_sphere, 'ls_best': ls_best,
        'fl_median': float(np.median(fl)) if fl else None,
        'grms_final': grms[-1],
    }

names = sorted(REF)
rows = [analyze(n) for n in names]

# summary by step count
rows_conv = [r for r in rows if r['conv']]
rows_noconv = [r for r in rows if not r['conv']]
print("=== NON-CONVERGED (ceiling/error) ===")
for r in sorted(rows_noconv, key=lambda x: -x['n']):
    print(f"{r['name']:20s} {r['family']:14s} atoms={r['atoms']:4d} n={r['n']:3d} "
          f"bottleneck={r['bottleneck']:16s} grms_final={r['grms_final']:.2e} "
          f"neg_max={r['neg_max']} trust_min={r['trust_min']}")

print("\n=== LONG CONVERGERS (n>=30) ===")
for r in sorted(rows_conv, key=lambda x: -x['n']):
    if r['n'] < 30: continue
    print(f"{r['name']:20s} {r['family']:12s} n={r['n']:3d} bott={r['bottleneck']:16s} "
          f"tail_frac={r['tail_frac']} e_settle={r['e_settle']}({r['e_settle_frac'] and round(r['e_settle_frac'],2)}) "
          f"onsph={r['on_sphere']:3d} lsbest={r['ls_best']:2d} neg_steps={r['neg_nonzero_steps']:3d} "
          f"trustmin={r['trust_min'] and round(r['trust_min'],3)} fl_med={r['fl_median'] and round(r['fl_median'],2)}")

# bottleneck tally across all converged
print("\n=== BOTTLENECK criterion tally (converged only) ===")
bt = defaultdict(int)
for r in rows_conv:
    bt[r['bottleneck']] += 1
for k, v in sorted(bt.items(), key=lambda x: -x[1]):
    print(f"  {k:18s}: {v}")

# bottleneck among long convergers
print("\n=== BOTTLENECK among long convergers (n>=30) ===")
bt = defaultdict(int)
for r in rows_conv:
    if r['n'] >= 30:
        bt[r['bottleneck']] += 1
for k, v in sorted(bt.items(), key=lambda x: -x[1]):
    print(f"  {k:18s}: {v}")

# tail fraction stats for long convergers
tails = [r['tail_frac'] for r in rows_conv if r['n'] >= 30 and r['tail_frac'] is not None]
print(f"\n=== TAIL fraction (n>=30): median={np.median(tails):.2f} mean={np.mean(tails):.2f} "
      f"min={min(tails):.2f} max={max(tails):.2f} ===")
print("(fraction of steps spent AFTER energy is within 0.1 kcal/mol of the final minimum)")

# per-family average steps vs atoms
print("\n=== per-family step scaling ===")
famrows = defaultdict(list)
for r in rows:
    famrows[r['family']].append(r)
for fam, rs in sorted(famrows.items()):
    rs = sorted(rs, key=lambda x: x['atoms'])
    steps = [f"{x['n']}{'*' if not x['conv'] else ''}" for x in rs]
    print(f"  {fam:16s}: {' '.join(steps)}")
