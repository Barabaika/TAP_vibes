"""Manual network test: one of the SAME ten antibodies; never runs in normal CI."""
import json
from pathlib import Path
import pytest
from tap2_web import run
from tap2_profiles import score

@pytest.mark.live
def test_new_service_request_and_recompute_all_five_metrics(tmp_path,psa):
    pairs=json.loads((Path(__file__).parent/'fixtures/pairs.json').read_text())
    pair=pairs[0]
    results=run([pair],tmp_path/'live',timeout=1800)
    reference=results[0]
    local=score(tmp_path/'live'/pair['id']/'web_model.pdb',psa,'web-compatible',
                expected={ch:pair[ch] for ch in 'HL'})
    for name,value in reference['metrics'].items():
        assert local['metrics'][name]==pytest.approx(value,abs=1e-4,rel=0),name
    print('Live TAP reference:',reference['url'])
