"""Real TAP responses, not expectations calculated by the implementation under test."""
import hashlib
import json
from pathlib import Path
import pytest
from compare import actual_flags
from tap2_profiles import score
from tap2_web import parse_result

FIXTURES=Path(__file__).parent/'fixtures'
REFERENCES=json.loads((FIXTURES/'web/web_results.json').read_text())

def test_exactly_ten_real_paired_references():
    assert len(REFERENCES)==10
    assert len({r['pair_sha256'] for r in REFERENCES})==10
    assert all(r['source_type']=='humanized' for r in REFERENCES)

@pytest.mark.parametrize('reference',REFERENCES,ids=[r['name'] or r['id'] for r in REFERENCES])
def test_reference_integrity_and_original_html(reference):
    directory=FIXTURES/'web'/reference['id']
    assert hashlib.sha256((reference['H']+'|'+reference['L']).encode()).hexdigest()==reference['pair_sha256']
    assert hashlib.sha256((directory/'web_model.pdb').read_bytes()).hexdigest()==reference['pdb_sha256']
    html=(directory/'result.html').read_text()
    metrics,_,links=parse_result(html)
    assert metrics==reference['metrics']
    assert len(actual_flags(html))==5
    assert links

@pytest.mark.structure
@pytest.mark.parametrize('reference',REFERENCES,ids=[r['name'] or r['id'] for r in REFERENCES])
def test_all_five_metrics_against_archived_web_service(reference,psa):
    pdb=FIXTURES/'web'/reference['id']/'web_model.pdb'
    result=score(pdb,psa,'web-compatible',expected={ch:reference[ch] for ch in 'HL'})
    assert set(result['metrics'])==set(reference['metrics'])
    for name,value in reference['metrics'].items():
        assert result['metrics'][name]==pytest.approx(value,abs=1e-4,rel=0), (reference['id'],name)
