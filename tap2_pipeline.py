"""Paired sequences -> ABodyBuilder2 refined Fv -> all five TAP2 metrics."""
import argparse, csv, hashlib, importlib.metadata, json, os, platform, random, time
from pathlib import Path
from tap2_web import atomic_json, read_pairs
from tap2_metrics import compute, METRICS
from tap2_profiles import score

def fold_pairs(pairs,out,weights,psa,device='cuda',threads=4,profile='paper'):
    import torch, numpy as np
    from ImmuneBuilder import ABodyBuilder2
    from tap2_compat import load_refine
    refine, refinement_compat = load_refine()
    from ImmuneBuilder.util import add_errors_as_bfactors
    if device=='cuda' and not torch.cuda.is_available(): raise RuntimeError('Requested CUDA is unavailable')
    if device=='cpu' and torch.cuda.is_available(): raise RuntimeError('Set CUDA_VISIBLE_DEVICES= for CPU mode')
    torch.set_num_threads(threads)
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    builder=ABodyBuilder2(weights_dir=str(weights),numbering_scheme='imgt')
    if profile not in ('paper','web-compatible'): raise ValueError(profile)
    provenance=dict(python=platform.python_version(),device=builder.device,profile=profile,
        refinement_compat=refinement_compat,
        packages={p:importlib.metadata.version(p) for p in ['ImmuneBuilder','torch','numpy','biopython','anarci','openmm','pdbfixer']},
        weights={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(weights).glob('antibody_model_*')},
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        metrics_sha256=hashlib.sha256(Path(__file__).with_name('tap2_metrics.py').read_bytes()).hexdigest(),
        profiles_sha256=hashlib.sha256(Path(__file__).with_name('tap2_profiles.py').read_bytes()).hexdigest(),
        psa_sha256=hashlib.sha256(Path(psa).read_bytes()).hexdigest())
    atomic_json(out/'provenance.json',provenance)
    results=[]
    for pair in pairs:
        d=out/pair['id']; d.mkdir(exist_ok=True)
        if (d/'result.json').exists():
            cached=json.loads((d/'result.json').read_text())
            if cached['pair_sha256']!=pair['pair_sha256']: raise ValueError('Cache sequence mismatch')
            if cached.get('provenance')!=provenance: raise ValueError('Cache protocol changed; use a fresh output directory')
            results.append(cached); continue
        atomic_json(d/'input.json',pair); start=time.monotonic()
        random.seed(2026); np.random.seed(2026); torch.manual_seed(2026)
        antibody=builder.predict({'H':pair['H'],'L':pair['L']})
        numbered={ch:[(list(n),aa) for n,aa in seq] for ch,seq in antibody.numbered_sequences.items()}
        atomic_json(d/'numbered_sequences.json',numbered)
        fv={ch:''.join(aa for n,aa in seq) for ch,seq in antibody.numbered_sequences.items()}
        if fv!={'H':pair['H'],'L':pair['L']}: raise ValueError('Input contains unnumbered flanks; provide Fv-only chains')
        np.save(d/'predicted_error_squared.npy',antibody.error_estimates.mean(0).cpu().numpy())
        final=d/'model.pdb'; success=False
        for rank,index in enumerate(antibody.ranking):
            unref=d/f'rank{rank}_unrefined.pdb'; antibody.save_single_unrefined(str(unref),index=index)
            for attempt in range(2):
                success=refine(str(unref),str(final),check_for_strained_bonds=True,n_threads=threads)
                if success: break
            if success: break
        if not success: raise RuntimeError(f'ABodyBuilder2 refinement failed: {pair["name"]}')
        add_errors_as_bfactors(str(final),antibody.error_estimates.mean(0).sqrt().cpu().numpy(),
                              header=['REMARK  MODELLED USING ABODYBUILDER2\n'])
        scores=score(final,psa,profile,expected=fv)
        result=dict(**pair,**scores,provenance=provenance,refined_rank=rank,elapsed_seconds=time.monotonic()-start)
        atomic_json(d/'result.json',result); results.append(result)
        atomic_json(out/'local_results.json',results)
        print('LOCAL_COMPLETE',pair['name'],scores['metrics'],f'{result["elapsed_seconds"]:.1f}s',flush=True)
    atomic_json(out/'local_results.json',results)
    write_metrics_csv(out/'metrics.csv',results)
    return results

def write_metrics_csv(path,rows):
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['id','name','pair_sha256',*METRICS]); w.writeheader()
        for r in rows: w.writerow({k:r[k] for k in ['id','name','pair_sha256']}|r['metrics'])

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    inp=ap.add_mutually_exclusive_group(required=True); inp.add_argument('--pairs-json'); inp.add_argument('--input')
    ap.add_argument('--kind',default='humanization'); ap.add_argument('--limit',type=int,default=10)
    ap.add_argument('--output',required=True); ap.add_argument('--weights',required=True); ap.add_argument('--psa',required=True)
    ap.add_argument('--device',choices=['cuda','cpu'],default='cuda'); ap.add_argument('--threads',type=int,default=4)
    ap.add_argument('--profile',choices=['paper','web-compatible'],default='paper')
    a=ap.parse_args()
    pairs=json.loads(Path(a.pairs_json).read_text()) if a.pairs_json else read_pairs(a.input,a.kind,a.limit)
    fold_pairs(pairs,a.output,a.weights,a.psa,a.device,a.threads,a.profile)

if __name__=='__main__':
    main()
