import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import json
from pathlib import Path

for root in ("ref", "ref15"):
    path = Path(_config_path('${PROJECT_ROOT}/work')) / root / 'annot/cache/fm_all.stats.json'
    if path.exists():
        obj = json.loads(path.read_text())
        print(root, 'top_level_fields', sorted(obj))
        for key in ('cols', 'columns', 'feature_names', 'annotations'):
            value = obj.get(key)
            if isinstance(value, list) and all(isinstance(s, str) for s in value):
                print(root, key, value)
