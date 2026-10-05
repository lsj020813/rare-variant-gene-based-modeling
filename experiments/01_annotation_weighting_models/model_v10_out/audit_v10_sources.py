import ast
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parent
files = ['l1_train_v8.py', 'l1_train_v9.py', 'l1_train_v10.py', 'v10_smoke.py']
report = {}
for name in files:
    path = root / name if name in ('l1_train_v10.py', 'v10_smoke.py') else root.parent / name.replace('l1_train_', 'model_').replace('.py', '_out') / name
    source = path.read_text()
    ast.parse(source)
    report[name] = dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        syntax='PASS', lines=len(source.splitlines()))
target = root / 'V10_SOURCE_AUDIT.json'
target.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
