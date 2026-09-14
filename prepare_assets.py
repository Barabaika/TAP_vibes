"""Fetch pinned PSA and optionally original ABodyBuilder2 weights; verify SHA256."""
import argparse
import hashlib
import json
from pathlib import Path
import requests

COMMIT='29dcac72f1380e8538e8870f45a699d3c6156162'
PSA_SHA='9f38d64c4e533ab5da8c7d608312065e8dc0fb7b41c7dd27216ce40f6fc220e3'
WEIGHTS={
    'antibody_model_1':'ecefff604457a105fe405dcca44a76b046609ff2ac3c355329974659f6881e09',
    'antibody_model_2':'91dd9f5ebca8cc68f292d3a40ccdbc8a86974e7cebc4853e19565ecb8a8d508b',
    'antibody_model_3':'33085a728cfb1059d9e817b785a21f1dda74eff8f083a0c2a033ae88bdbdcd45',
    'antibody_model_4':'e2d666a22b1f0a7219baec778e67f36e448ad0cd8eed8dae044af20eec163a04',
}

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def download(url,path,sha=None):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        if sha and digest(path)!=sha: raise ValueError(f'Checksum mismatch: {path}')
        return
    temporary=path.with_suffix('.part')
    with requests.get(url,stream=True,timeout=(30,120)) as response:
        response.raise_for_status()
        with temporary.open('wb') as f:
            for chunk in response.iter_content(1024*1024): f.write(chunk)
    if sha and digest(temporary)!=sha: raise ValueError(f'Download checksum mismatch: {url}')
    temporary.replace(path)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('assets'))
    parser.add_argument('--psa-only',action='store_true',help='Do not download the 674 MB folding weights')
    args=parser.parse_args(); root=args.output
    base=f'https://raw.githubusercontent.com/oxpig/TNP/{COMMIT}'
    download(base+'/bin/psa',root/'psa',PSA_SHA)
    download(base+'/LICENCE',root/'TNP-LICENCE')
    (root/'psa').chmod(0o755)
    manifest={'psa_commit':COMMIT,'psa_sha256':PSA_SHA,'weights':{}}
    if not args.psa_only:
        for name,sha in WEIGHTS.items():
            url=f'https://zenodo.org/records/7258553/files/{name}?download=1'
            path=root/'weights'/name
            download(url,path,sha)
            manifest['weights'][name]={'url':url,'sha256':sha,'bytes':path.stat().st_size}
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(manifest,indent=2))

if __name__=='__main__': main()
