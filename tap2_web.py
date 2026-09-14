"""Small, resumable TAP2 web-service validation client (explicit opt-in, max 10)."""
import argparse, csv, hashlib, json, os, re, time
from pathlib import Path
from urllib.parse import urljoin, urlparse
import requests
from bs4 import BeautifulSoup

URL = 'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabpred/tap'
METRICS = ['L_tot','PSH','PPC','PNC','SFvCSP']
LABELS = ['Total CDR Length','CDR Vicinity PSH Score (Kyte & Doolittle)',
          'CDR Vicinity PPC Score','CDR Vicinity PNC Score','SFvCSP Score']

def atomic_json(path, obj):
    path=Path(path); tmp=path.with_suffix('.tmp')
    tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8'); tmp.replace(path)

def read_pairs(path, kind='humanized', limit=10):
    if limit<0: raise ValueError('limit must be non-negative (0 means all)')
    pairs=[]; seen=set()
    with open(path) as f:
        for index,r in enumerate(csv.DictReader(f)):
            row_kind=r.get('Specific',r.get('type',''))
            if kind and row_kind != kind: continue
            h=r.get('hseq',r.get('h_seq',r.get('H',''))).strip().upper()
            l=r.get('lseq',r.get('l_seq',r.get('L',''))).strip().upper()
            if not h or not l or set(h+l)-set('ACDEFGHIKLMNPQRSTVWY'):
                raise ValueError(f'Invalid paired sequence at row {index+2}')
            key=hashlib.sha256((h+'|'+l).encode()).hexdigest()
            if key in seen: continue
            seen.add(key)
            pairs.append(dict(id=f'{len(pairs):03d}_{key[:12]}',name=r.get('name',str(index)),
                              H=h,L=l,pair_sha256=key,source=str(path),
                              source_row=index+2,source_type=row_kind))
            if limit and len(pairs)>=limit: break
    if not pairs: raise ValueError('No matching pairs')
    return pairs

def parse_result(html):
    soup=BeautifulSoup(html,'html.parser'); values={}; flags={}
    for tr in soup.select('tr'):
        cells=tr.find_all(['th','td'])
        if len(cells)<2: continue
        label=cells[0].get_text(' ',strip=True)
        if label in LABELS:
            metric=METRICS[LABELS.index(label)]
            text=cells[1].get_text(' ',strip=True)
            values[metric]=float(text)
            flags[metric]=str(cells[1].get('style',''))
    links=[urljoin(URL,a['href']) for a in soup.find_all('a',href=True)
           if '/tap_results/' in a['href'] and a['href'].endswith('/model')]
    return values,flags,links

def run(pairs,out,timeout=1800):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    if len(pairs)>10: raise ValueError('Web validation is restricted to at most 10 antibodies')
    session=requests.Session()
    session.headers['User-Agent']='TAP2-validation/1.0 (10 paired antibodies; sequential requests)'
    home=session.get(URL,timeout=60); home.raise_for_status()
    (out/'service_home.html').write_text(home.text)
    results=[]
    for p in pairs:
        d=out/p['id']; d.mkdir(exist_ok=True)
        if (d/'input.json').exists():
            previous=json.loads((d/'input.json').read_text(encoding='utf-8'))
            if previous['pair_sha256']!=p['pair_sha256']:
                raise ValueError('Web cache sequence mismatch; use a new output directory')
        atomic_json(d/'input.json',p)
        if (d/'result.json').exists():
            cached=json.loads((d/'result.json').read_text(encoding='utf-8'))
            if cached['pair_sha256']!=p['pair_sha256']:
                raise ValueError('Web result sequence mismatch')
            results.append(cached); continue
        if (d/'submission.json').exists():
            sub=json.loads((d/'submission.json').read_text())
        else:
            # Never automatically retry an uncertain POST: it can create duplicate jobs.
            if (d/'submission_started.json').exists() or (d/'submission_uncertain.json').exists():
                raise RuntimeError('Uncertain prior POST; inspect and recover its job URL before retry')
            atomic_json(d/'submission_started.json',dict(time=time.time(),pair=p['pair_sha256']))
            try:
                r=session.post(URL,data={'hchain':p['H'],'lchain':p['L']},timeout=90)
                r.raise_for_status()
                (d/'submission.html').write_text(r.text)
                if '/tap_results/' not in r.url:
                    raise RuntimeError(f'No TAP job URL returned: {r.url}')
                sub=dict(url=r.url,time=time.time())
                atomic_json(d/'submission.json',sub)
            except Exception as e:
                atomic_json(d/'submission_uncertain.json',dict(error=str(e)))
                raise
        print('WEB',p['name'],sub['url'],flush=True)
        started=time.monotonic()
        while True:
            r=session.get(sub['url'],timeout=90); r.raise_for_status()
            (d/'result.html').write_text(r.text)
            metrics,flags,links=parse_result(r.text)
            if len(metrics)==5 and links:
                model=session.get(links[0],timeout=90); model.raise_for_status()
                if b'ATOM ' not in model.content: raise RuntimeError('Downloaded model is not PDB')
                (d/'web_model.pdb').write_bytes(model.content)
                result=dict(**p,url=sub['url'],metrics=metrics,web_flag_styles=flags,
                            pdb_sha256=hashlib.sha256(model.content).hexdigest(),retrieved_at=time.time())
                atomic_json(d/'result.json',result); results.append(result)
                print('WEB_COMPLETE',p['name'],metrics,flush=True); break
            soup=BeautifulSoup(r.text,'html.parser')
            message=soup.get_text(' ',strip=True)
            if re.search(r'Job (?:failed|error)|TAP failed',message,re.I):
                raise RuntimeError(f'TAP job failed: {sub["url"]}')
            if time.monotonic()-started>timeout: raise TimeoutError(sub['url'])
            time.sleep(20)
        atomic_json(out/'web_results.json',results)
    atomic_json(out/'web_results.json',results)
    return results

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input',required=True)
    ap.add_argument('--kind',default='humanized'); ap.add_argument('--limit',type=int,default=10)
    ap.add_argument('--output',required=True); ap.add_argument('--submit-web',action='store_true')
    args=ap.parse_args()
    if not 1<=args.limit<=10: ap.error('Validation limit must be 1..10')
    out=Path(args.output)
    out.mkdir(parents=True,exist_ok=True)
    pairs=read_pairs(args.input,args.kind,args.limit); atomic_json(out/'pairs.json',pairs)
    if args.submit_web: run(pairs,out/'web')

if __name__=='__main__':
    main()
