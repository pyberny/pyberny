#!/usr/bin/env python3
"""Confirm the soft-mode tail mechanism: in the wasted tail (forces already
converged) is the predicted energy change negligible while the step still
exceeds stepmax, and is the step limited by the trust sphere / soft modes?"""
import json, os
import numpy as np
from berny.berny import BernyParams

TR = os.environ.get('OLIGO_TRACES','./traces')
P = BernyParams(); HK = 627.5

def load(n): return json.load(open(f'{TR}/xtb-oligomers-{n}.trace.json'))
def crit(r, nm):
    for c in r['convergence']['criteria']:
        if c['name'] == nm: return c

for name in ['polyserine_n5', 'polyglycine_n6', 'nylon6_n4']:
    t = load(name); n = len(t)
    print(f"\n===== {name} (n={n}) — last 12 steps =====")
    print(f"{'step':>4}{'gRMS':>10}{'gMax':>10}{'stepMax':>10}{'onSph':>6}{'trust':>8}{'lowEV':>9}{'|dE_pred|kcal':>13}{'conv?':>6}")
    for r in t[-12:]:
        g1 = crit(r,'Gradient RMS'); g2 = crit(r,'Gradient maximum')
        smax = crit(r,'Step maximum')
        q = r.get('quadratic_step', {})
        dEp = abs(q.get('predicted_energy_change', 0))*HK
        smax_v = smax['value'] if smax else float('nan')
        forces_ok = g1['value']<P.gradientrms and g2['value']<P.gradientmax
        print(f"{r['step']:>4}{g1['value']:>10.2e}{g2['value']:>10.2e}{smax_v:>10.2e}"
              f"{str(q.get('on_sphere')):>6}{q.get('trust_radius',0):>8.3f}"
              f"{q.get('lowest_eigenvalue',0):>9.4f}{dEp:>13.5f}"
              f"{'*'+('C' if r['converged'] else '') if forces_ok else '':>6}")

# Aggregate: across ALL wasted-tail steps (forces converged, not yet converged),
# what is the predicted energy change and step/threshold ratio?
print("\n===== Aggregate over all wasted-tail steps (forces already < threshold) =====")
dEs=[]; ratios=[]; lowevs=[]; onsph=0; tot=0
for name in json.load(open('/home/user/pyberny/src/berny/benchmarks/oligomers/reference.json')):
    t = load(name); n=len(t)
    # find contiguous force-ok tail
    ok_from=None
    for i in range(n):
        g1=crit(t[i],'Gradient RMS'); g2=crit(t[i],'Gradient maximum')
        if all(crit(t[j],'Gradient RMS')['value']<P.gradientrms and
               crit(t[j],'Gradient maximum')['value']<P.gradientmax for j in range(i,n)):
            ok_from=i; break
    if ok_from is None: continue
    for r in t[ok_from:n-1]:  # exclude the final converged step
        q=r.get('quadratic_step',{}); smax=crit(r,'Step maximum')
        if smax is None:  # on-sphere
            onsph+=1; tot+=1; continue
        tot+=1
        dEs.append(abs(q.get('predicted_energy_change',0))*HK)
        ratios.append(smax['value']/P.stepmax)
        lowevs.append(q.get('lowest_eigenvalue',0))
print(f"wasted-tail steps analysed: {tot} ({onsph} on-sphere)")
print(f"predicted |dE| in tail steps (kcal/mol): median={np.median(dEs):.5f} max={max(dEs):.4f}")
print(f"stepMax / threshold ratio in tail: median={np.median(ratios):.1f}x  (>1 means tripping the criterion)")
print(f"lowest Hessian eigenvalue in tail: median={np.median(lowevs):.4f} min={min(lowevs):.4f}")
print(f"(a soft mode lambda~0.005 turns a residual g~1e-4 into dq~g/lambda~0.02 a.u. >> stepmax {P.stepmax})")
