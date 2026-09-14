"""Narrow, in-memory workaround for ImmuneBuilder 1.2 CPU refinement typo."""
import hashlib
import importlib
from pathlib import Path

BAD = "{'Threads', str(n_threads)}"
GOOD = "{'Threads': str(n_threads)}"

def load_refine():
    module = importlib.import_module('ImmuneBuilder.refine')
    source = Path(module.__file__).read_text()
    info = {'upstream_sha256': hashlib.sha256(source.encode()).hexdigest(),
            'cpu_threads_set_to_dict': BAD in source}
    if BAD not in source:
        return module.refine, info
    if source.count(BAD) != 1:
        raise RuntimeError('Unexpected ImmuneBuilder CPU typo count; inspect dependency')
    # Execute a private module namespace; installed ImmuneBuilder remains unchanged.
    namespace = {'__name__': 'tap2_immune_refine', '__file__': module.__file__,
                 '__package__': module.__package__}
    exec(compile(source.replace(BAD, GOOD), module.__file__ + ':tap2-cpu-fix', 'exec'), namespace)
    return namespace['refine'], info
