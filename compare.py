import argparse, csv, json, os, re
from pathlib import Path
import numpy as np
from bs4 import BeautifulSoup
from tap2_profiles import score
from tap2_metrics import METRICS
from tap2_web import atomic_json

def actual_flags(html):
    text=BeautifulSoup(html,'html.parser').get_text(' ',strip=True).replace('\xa0',' ')
    labels=['Total IMGT CDR Length:', 'Patch CDR Surface Hydrophobicity score is:',
            'Patch CDR Positive Charge score is:', 'Patch CDR Negative Charge score is:', 'SFvCSP score is:']
    flags={}
    for k,label in zip(METRICS,labels):
        m=re.search(re.escape(label)+r'\s*[-\d.]+\s*\((GREEN|AMBER|RED) flag\)',text)
        if not m: raise ValueError(f'Missing actual web flag: {k}')
        flags[k]=m[1]
    return flags

def main(web,local,out,psa):
    out.mkdir(parents=True,exist_ok=True)
    refs=json.loads((web/'web_results.json').read_text())
    if len(refs)!=10 or len({r['pair_sha256'] for r in refs})!=10: raise ValueError('Expected exactly 10 unique pairs')
    rows=[]; detailed=[]
    for ref in refs:
        id=ref['id']; d=web/id
        label=ref['name'].strip() or f'HuAb348 row {ref["source_row"]} (unnamed)'
        flags=actual_flags((d/'result.html').read_text())
        expected={ch:ref[ch] for ch in 'HL'}
        same_web=score(d/'web_model.pdb',psa,'web-compatible',expected)
        same_paper=score(d/'web_model.pdb',psa,'paper',expected)
        local_web=score(local/id/'model.pdb',psa,'web-compatible',expected)
        local_paper=score(local/id/'model.pdb',psa,'paper',expected)
        for profile,value in [('same_pdb_web',same_web),('same_pdb_paper',same_paper),
                              ('local_fold_web',local_web),('local_fold_paper',local_paper)]:
            atomic_json(out/f'{id}_{profile}.json',value)
            for k in METRICS:
                rows.append(dict(id=id,name=label,profile=profile,metric=k,value=value['metrics'][k],
                    web_value=ref['metrics'][k],delta=value['metrics'][k]-ref['metrics'][k],
                    actual_web_flag=flags[k],url=ref['url']))
        detailed.append(dict(id=id,name=label,source=ref['source'],source_row=ref['source_row'],url=ref['url'],
            web=ref['metrics'],web_flags=flags,same_pdb=same_web['metrics'],paper_on_web_pdb=same_paper['metrics'],
            local=local_web['metrics'],paper_on_local_pdb=local_paper['metrics'],
            paper_flags_local=local_paper['paper_2024_flags']))
    summary={}
    for profile in ('same_pdb_web','same_pdb_paper','local_fold_web','local_fold_paper'):
        summary[profile]={}
        for k in METRICS:
            rr=[r for r in rows if r['profile']==profile and r['metric']==k]
            delta=np.array([r['delta'] for r in rr])
            summary[profile][k]=dict(n=len(rr),mae=float(abs(delta).mean()),max_abs=float(abs(delta).max()),
                 rmse=float(np.sqrt(np.mean(delta**2))))
    summary['same_pdb_tolerance']=1e-4
    summary['same_pdb_pass']=all(v['max_abs']<=1e-4 for v in summary['same_pdb_web'].values())
    summary['independent_fold_note']='Structural variation measured descriptively; no arbitrary pass/fail tolerance.'
    atomic_json(out/'summary.json',summary); atomic_json(out/'comparison.json',detailed)
    with (out/'comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(json.dumps(summary,indent=2),flush=True)
    if not summary['same_pdb_pass']: raise AssertionError('Same-PDB web compatibility regression failed')

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--web',type=Path,required=True); p.add_argument('--local',type=Path,required=True)
    p.add_argument('--psa',required=True); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); main(a.web,a.local,a.output,a.psa)
