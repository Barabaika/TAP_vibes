import os
from pathlib import Path
import shutil
import pytest

def pytest_addoption(parser):
    parser.addoption('--run-structure',action='store_true',help='Recompute all five metrics on 10 archived service PDBs')
    parser.addoption('--run-live',action='store_true',help='Submit the first archived antibody to the public TAP service')
    parser.addoption('--psa',default=os.environ.get('TAP_PSA','assets/psa'))

def pytest_collection_modifyitems(config,items):
    for item in items:
        for marker,option in [('structure','--run-structure'),('live','--run-live')]:
            if marker in item.keywords and not config.getoption(option):
                item.add_marker(pytest.mark.skip(reason=f'Explicit opt-in required: {option}'))

@pytest.fixture(scope='session')
def psa(request):
    path=Path(request.config.getoption('--psa')).resolve()
    if not path.is_file(): pytest.fail(f'Missing PSA: {path}; run tap-assets --psa-only')
    if not shutil.which('hmmscan'): pytest.fail('HMMER hmmscan must be on PATH')
    return path
