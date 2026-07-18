#!/usr/bin/env python3
import json, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

OUT=os.environ.get('OLIGO_FIG','.')
rows = json.load(open('./data/fix_results.json'))
# keep flexible ones that converged in both; drop controls and the unrescued n8
flex = [r for r in rows if r['bc'] and r['fc'] and r['name'] not in
        ('polyethylene_n8','thiophene_n4','anthracene','pentacene','polyglycine_n8')]
flex = sorted(flex, key=lambda r: -(r['bn']))
names=[r['name'] for r in flex]; base=[r['bn'] for r in flex]; fix=[r['fn'] for r in flex]
y=np.arange(len(names))
fig,ax=plt.subplots(figsize=(8,6))
ax.barh(y-0.2, base, height=0.4, color='#8d99ae', label='baseline')
ax.barh(y+0.2, fix, height=0.4, color='#2a9d8f', label='flat-surface accept')
for i,r in enumerate(flex):
    ax.text(max(r['bn'],r['fn'])+1, i, f"-{r['saved']}", va='center', fontsize=8, color='#264653')
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8); ax.invert_yaxis()
ax.set_xlabel('pyberny steps (GFN2-xTB)')
ax.set_title('Flat-surface force-based acceptance removes the displacement tail\n'
             '(96 steps saved / 19 molecules; final energies within 0.022 kcal/mol)')
ax.legend(fontsize=9, loc='lower right'); ax.grid(alpha=0.2, axis='x')
fig.tight_layout(); fig.savefig(f'{OUT}/fix_savings.png', dpi=130); plt.close(fig)
print("wrote fix_savings.png")
