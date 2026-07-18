#!/usr/bin/env python3
import json, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

TR = os.environ.get('OLIGO_TRACES','./traces')
OUT = os.environ.get('OLIGO_FIG','.')
os.makedirs(OUT, exist_ok=True)
REF = json.load(open('/home/user/pyberny/src/berny/benchmarks/oligomers/reference.json'))
from berny.berny import BernyParams
P = BernyParams()

def load(n): return json.load(open(f'{TR}/xtb-oligomers-{n}.trace.json'))
def crit(r, nm):
    for c in r['convergence']['criteria']:
        if c['name'] == nm: return c

FLEX = {'glycine','alanine','serine','nylon'}
RIGID = {'acenes','polyynes','thiophene','ppe'}
ALI = {'ethylene','propylene','ethylene_glycol'}

# ---------- Fig 1: steps vs atoms ----------
fig, ax = plt.subplots(figsize=(7,5))
colors = {'flexible chains (peptides, nylon)':'#d1495b','rigid π (acenes, thiophene, PPE, polyynes)':'#2e86ab','aliphatic (PE, PP, PEG)':'#5a5a5a'}
def bucket(fam):
    if fam in FLEX: return 'flexible chains (peptides, nylon)'
    if fam in RIGID: return 'rigid π (acenes, thiophene, PPE, polyynes)'
    return 'aliphatic (PE, PP, PEG)'
seen=set()
for name in sorted(REF):
    n = len(load(name)); conv = load(name)[-1]['converged']
    atoms = REF[name]['atoms']; b = bucket(REF[name]['family'])
    lab = b if b not in seen else None; seen.add(b)
    mk = 'o' if conv else 'x'
    ax.scatter(atoms, n, c=colors[b], marker=mk, s=45, label=lab, alpha=0.8,
               edgecolors='none' if conv else colors[b], linewidths=1.5)
ax.axhline(100, ls='--', c='k', lw=0.8, alpha=0.6)
ax.text(6, 103, 'default maxsteps=100', fontsize=8)
ax.set_xlabel('number of atoms'); ax.set_ylabel('pyberny steps (GFN2-xTB)')
ax.set_title('Step count vs system size, oligomers benchmark\n(× = hit 100-step ceiling / errored)')
ax.legend(fontsize=8, loc='upper left'); ax.grid(alpha=0.2)
fig.tight_layout(); fig.savefig(f'{OUT}/steps_vs_atoms.png', dpi=130); plt.close(fig)

# ---------- Fig 2: convergence trajectory of polyglycine_n6 ----------
name='polyglycine_n6'; t=load(name)
steps=[r['step'] for r in t]
gmax=[crit(r,'Gradient maximum')['value'] for r in t]
grms=[crit(r,'Gradient RMS')['value'] for r in t]
smax=[(crit(r,'Step maximum') or {}).get('value',np.nan) for r in t]
srms=[(crit(r,'Step RMS') or {}).get('value',np.nan) for r in t]
fig, ax = plt.subplots(figsize=(7.5,5))
ax.semilogy(steps,gmax,'-',c='#2e86ab',label='gradient max')
ax.semilogy(steps,grms,'-',c='#6bb8d6',label='gradient RMS')
ax.semilogy(steps,smax,'-',c='#d1495b',label='step max (displacement)')
ax.semilogy(steps,srms,'-',c='#e69aa8',label='step RMS')
ax.axhline(P.gradientmax, ls=':', c='#2e86ab', lw=1); ax.axhline(P.gradientrms, ls=':', c='#6bb8d6', lw=1)
ax.axhline(P.stepmax, ls=':', c='#d1495b', lw=1); ax.axhline(P.steprms, ls=':', c='#e69aa8', lw=1)
# shade the wasted tail
ax.axvspan(55, 66, color='#ffd166', alpha=0.3)
ax.text(56, 3e-4, 'forces converged;\ndisplacement tail\n(11 steps, ΔE<0.03 kcal/mol)', fontsize=8)
ax.set_xlabel('step'); ax.set_ylabel('convergence metric (a.u., log)')
ax.set_title(f'{name}: forces converge at step ~55, displacement lags to step 66')
ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.2, which='both')
fig.tight_layout(); fig.savefig(f'{OUT}/trajectory_polyglycine_n6.png', dpi=130); plt.close(fig)

print("figures written:", os.listdir(OUT))
