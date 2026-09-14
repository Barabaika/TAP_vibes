"""Score an existing IMGT-numbered paired Fv PDB without structure prediction."""
import argparse,json
from pathlib import Path
from tap2_profiles import score

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('pdb'); p.add_argument('--psa',required=True)
    p.add_argument('--profile',choices=['paper','web-compatible'],default='paper')
    p.add_argument('--output',required=True)
    a=p.parse_args(); r=score(a.pdb,a.psa,a.profile)
    Path(a.output).parent.mkdir(parents=True,exist_ok=True)
    Path(a.output).write_text(json.dumps(r,indent=2),encoding='utf-8')
    print(json.dumps(r['metrics'],indent=2))

if __name__=='__main__':
    main()
