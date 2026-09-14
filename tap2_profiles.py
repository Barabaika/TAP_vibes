"""Explicit separation of the published TAP2 specification and observed web behavior."""
import math
from tap2_metrics import compute, CHARGE

PAPER_THRESHOLDS={
 'L_tot':{'red_low':37.,'amber_low':42.,'amber_high':55.,'red_high':63.},
 'PSH':{'red_low':95.58,'amber_low':110.11,'amber_high':168.06,'red_high':201.59},
 'PPC':{'amber_high':1.32,'red_high':4.22},
 'PNC':{'amber_high':2.00,'red_high':4.42},
 'SFvCSP':{'red_low':-30.60,'amber_low':-6.00},
}

def paper_flags(metrics):
    result={}
    for k,v in metrics.items():
        if not math.isfinite(v): raise ValueError(f'Non-finite metric {k}')
        t=PAPER_THRESHOLDS[k]
        red=('red_low' in t and v<t['red_low']) or ('red_high' in t and v>t['red_high'])
        amber=('amber_low' in t and v<=t['amber_low']) or ('amber_high' in t and v>=t['amber_high'])
        result[k]='RED' if red else 'AMBER' if amber else 'GREEN'
    return result

def score(path,psa,profile='paper',expected=None):
    if profile not in ('paper','web-compatible'): raise ValueError(profile)
    r=compute(path,psa,vicinity_radius=4.5 if profile=='paper' else 4.0,
              neutralize_patch_charges=profile=='paper',expected=expected)
    if profile=='web-compatible':
        # Observed on 10 independent service jobs dated 2026-09-14, not an
        # assertion about the implementation of the separately licensed TAP.
        q={ch:sum(CHARGE.get(x['aa'],0.) for x in r['residues'] if x['exposed'] and x['chain']==ch)
           for ch in 'HL'}
        r['VH_surface_charge']=q['H']; r['VL_surface_charge']=q['L']
        r['metrics']['SFvCSP']=q['H']*q['L']
        r['protocol']['name']='TAP-web-observed-2026-09-14'
        r['protocol']['neutralize_sfv_charge']=False
        # Do not label thresholds from the inconsistent service homepage as
        # the actual service flags. Actual flags are retained in web responses.
    else:
        r['protocol']['name']='TAP2-paper-specification-2024'
        r['protocol']['neutralize_sfv_charge']=True
        r['paper_2024_flags']=paper_flags(r['metrics'])
        r['flag_thresholds_source']='https://www.nature.com/articles/s42003-023-05744-8/tables/1'
    return r
