"""TAP2 descriptors on paired, IMGT-numbered Fv structures.

Raybould et al. PNAS 2019, doi:10.1073/pnas.1810576116;
Raybould et al. Commun Biol 2024, doi:10.1038/s42003-023-05744-8.
SASA uses the PSA binary redistributed by OPIG TNP (BSD-3-Clause).
"""
from pathlib import Path
import hashlib, json, math, subprocess, tempfile
import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from scipy.spatial.distance import cdist

METRICS=('L_tot','PSH','PPC','PNC','SFvCSP')
CDRS=((27,38),(56,65),(105,117))
KD=dict(zip('ACDEFGHIKLMNPQRSTVWY',[1.8,2.5,-3.5,-3.5,2.8,-0.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-0.8,-0.7,4.2,-0.9,-1.3]))
CHARGE={'K':1.,'R':1.,'H':.1,'D':-1.,'E':-1.}
DONOR={'K':('NZ',),'R':('NH1','NH2')}
ACCEPTOR={'D':('OD1','OD2'),'E':('OE1','OE2')}

def in_cdr(number,padding=0):
    return any(a-padding<=number<=b+padding for a,b in CDRS)

def load_fv(path):
    structure=PDBParser(QUIET=True).get_structure('fv',str(path))
    if len(structure)!=1: raise ValueError('Expected exactly one structural model')
    model=structure[0]
    if not all(ch in model for ch in 'HL'): raise ValueError('PDB must contain H and L chains')
    residues=[]
    for ch in 'HL':
        for r in model[ch]:
            if r.id[0]!=' ': continue
            aa=seq1(r.resname)
            if aa not in KD: raise ValueError(f'Unsupported residue: {r.resname}')
            atoms={a.name:np.asarray(a.coord,dtype=np.float64) for a in r
                   if a.element not in ('H','D') and not a.name.lstrip('0123456789').startswith('H')}
            if not {'N','CA','C','O'}<=atoms.keys(): raise ValueError(f'Missing backbone atoms: {ch}{r.id}')
            residues.append(dict(chain=ch,number=r.id[1],insertion=r.id[2].strip(),aa=aa,
                                 name=r.resname,atoms=atoms))
    return residues

def check_imgt(residues, expected=None):
    from anarci import anarci
    for ch in 'HL':
        rr=[r for r in residues if r['chain']==ch]
        seq=''.join(r['aa'] for r in rr)
        if expected and seq!=expected[ch]: raise ValueError(f'{ch} structure sequence differs from expected Fv')
        numbered,details,_=anarci([(ch,seq)],scheme='imgt')
        if not numbered[0] or len(numbered[0])!=1: raise ValueError(f'Expected one {ch} variable domain')
        typ=details[0][0]['chain_type']
        if (ch=='H' and typ!='H') or (ch=='L' and typ not in ('K','L')): raise ValueError('Chain identity mismatch')
        nums=[(n[0],n[1].strip(),aa) for n,aa in numbered[0][0][0] if aa!='-']
        actual=[(r['number'],r['insertion'],r['aa']) for r in rr]
        if nums!=actual: raise ValueError(f'{ch} PDB is not correctly IMGT numbered')
        for a,b in CDRS:
            if not any(a<=r['number']<=b for r in rr): raise ValueError('Missing CDR')

def sasa(residues,psa):
    """PSA reports relative SIDECHAIN area vs Ala-X-Ala, not total residue RSA.

    Use unique sequential IDs only in a temporary SASA copy, preserving chain
    identity, coordinates and inter-chain occlusion. Return values by row order.
    """
    with tempfile.TemporaryDirectory(prefix='tap2_sasa_') as td:
        p=Path(td)/'fv.pdb'; lines=[]; serial=0
        for i,r in enumerate(residues,1):
            for name,xyz in r['atoms'].items():
                serial+=1; x,y,z=xyz; element=name[0]
                atom_field=f' {name:<3s}' if len(name)<4 else name
                lines.append(f'ATOM  {serial:5d} {atom_field} {r["name"]:3s} {r["chain"]}{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2s}\n')
        p.write_text(''.join(lines)+'END\n')
        result=subprocess.run([str(psa),'-t',str(p)],capture_output=True,text=True,check=True)
        parsed={}
        for line in result.stdout.splitlines():
            if line.startswith('ACCESS'):
                i=int(line[6:11])-1
                if i in parsed: raise ValueError('Duplicate PSA residue')
                parsed[i]=(float(line[55:61]),float(line[61:67]))
        if set(parsed)!=set(range(len(residues))):
            raise ValueError(f'PSA did not return every residue: {len(parsed)}/{len(residues)}. {result.stderr[:500]}')
        values=np.array([parsed[i] for i in range(len(residues))])
        if not np.isfinite(values).all(): raise ValueError('Non-finite PSA values')
        return values,result.stdout

def compute(path,psa,vicinity_radius=4.5,verify_numbering=True,expected=None,neutralize_patch_charges=True):
    residues=load_fv(path)
    if verify_numbering: check_imgt(residues,expected)
    areas,raw=sasa(residues,psa)
    exposed=areas[:,1]>=7.5
    cdr=np.array([in_cdr(r['number']) for r in residues])
    anchor=np.array([in_cdr(r['number'],2) for r in residues])
    n=len(residues); distances=np.full((n,n),np.inf)
    for i in range(n):
        for j in range(i+1,n):
            distances[i,j]=distances[j,i]=cdist(list(residues[i]['atoms'].values()),list(residues[j]['atoms'].values())).min()
    if np.any(distances<0.1): raise ValueError('Severe inter-residue atom overlap')
    seed=exposed & anchor
    vicinity=exposed & (seed | (distances[:,seed]<vicinity_radius).any(axis=1))
    q=np.array([CHARGE.get(r['aa'],0.) for r in residues]); original_q=q.copy()
    hyd=np.array([1+(KD[r['aa']]+4.5)/9. for r in residues])
    bridges=[]
    for i,r in enumerate(residues):
        if not exposed[i] or r['aa'] not in DONOR: continue
        for j,s in enumerate(residues):
            if not exposed[j] or s['aa'] not in ACCEPTOR: continue
            donor=[r['atoms'][a] for a in DONOR[r['aa']] if a in r['atoms']]
            acceptor=[s['atoms'][a] for a in ACCEPTOR[s['aa']] if a in s['atoms']]
            if not donor or not acceptor: raise ValueError('Missing salt-bridge side-chain atoms')
            d=float(cdist(donor,acceptor).min())
            if d<=3.2:
                q[[i,j]]=0.; hyd[[i,j]]=1+(KD['G']+4.5)/9.
                bridges.append(dict(i=i,j=j,distance=d))
    pairmask=vicinity[:,None]&vicinity[None,:]&(distances<7.5)
    # TAP convention: ordered pairs i != j (both i,j and j,i).
    weight=np.where(pairmask,1./distances**2,0.)
    vh=float(q[exposed & np.array([r['chain']=='H' for r in residues])].sum())
    vl=float(q[exposed & np.array([r['chain']=='L' for r in residues])].sum())
    if not neutralize_patch_charges: q=original_q
    metrics=dict(L_tot=int(cdr.sum()),PSH=float(hyd@weight@hyd),
                 PPC=float(np.maximum(q,0)@weight@np.maximum(q,0)),
                 PNC=float(np.minimum(q,0)@weight@np.minimum(q,0)),SFvCSP=vh*vl)
    annotations=[]
    for i,r in enumerate(residues):
        annotations.append({k:v for k,v in r.items() if k!='atoms'}|dict(
            sidechain_sasa=float(areas[i,0]),relative_sidechain_sasa=float(areas[i,1]),
            exposed=bool(exposed[i]),cdr=bool(cdr[i]),cdr_vicinity=bool(vicinity[i]),
            charge=float(q[i]),hydrophobicity=float(hyd[i])))
    return dict(metrics=metrics,VH_surface_charge=vh,VL_surface_charge=vl,
                cdr_lengths={f'{ch}{k+1}':sum(r['chain']==ch and a<=r['number']<=b for r in residues)
                             for ch in 'HL' for k,(a,b) in enumerate(CDRS)},
                salt_bridges=bridges,residues=annotations,psa_output=raw,
                protocol=dict(name='TAP2-paper-2024',vicinity_radius_A=vicinity_radius,
                    exposure_percent=7.5,patch_cutoff_A=7.5,salt_bridge_cutoff_A=3.2,
                    pair_sum='ordered',hydrophobicity='Kyte-Doolittle normalized [1,2]',
                    neutralize_patch_charges=neutralize_patch_charges,
                    psa_sha256=hashlib.sha256(Path(psa).read_bytes()).hexdigest()),
                pdb_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest())
