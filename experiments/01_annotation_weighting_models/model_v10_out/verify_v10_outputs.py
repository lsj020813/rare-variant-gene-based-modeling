import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SUITE = ROOT / '.smoke/suite_03'
APPLICATION = ROOT / '.smoke/apply_cli'
CLI = ROOT / '.smoke/cli_train'

def sha(path, compressed=False):
    result = hashlib.sha256()
    opener = gzip.open if compressed else open
    with opener(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()

def read(path):
    return json.loads(path.read_text())

marker = read(SUITE / 'SMOKE_SUITE_PASS')
result = read(SUITE / 'RESULT.json')
checks = result['smoke_checks']
assert all(c['pass_'] for c in checks.values())
assert marker['code_sha256'] == sha(ROOT / 'l1_train_v10.py')
assert marker['result_sha256'] == sha(SUITE / 'RESULT.json')
assert read(SUITE / 'L1_DONE')['result_sha256'] == marker['result_sha256']
assert read(APPLICATION / 'L1_DONE')['result_sha256'] == sha(APPLICATION / 'RESULT.json')
assert read(CLI / 'L1_DONE')['result_sha256'] == sha(CLI / 'RESULT.json')
assert read(CLI / 'RESULT.json')['verdict'].get('unsealed') == ['alpha']
identical_export = sha(SUITE / 'phi_scores_1-5.tsv.gz', True) == sha(APPLICATION / 'phi_scores_1-5.tsv.gz', True)
assert identical_export
fit_status = []
for path in SUITE.rglob('*.selection.json'):
    fit_status.extend(c['optimizer']['success'] for c in read(path)['candidates'].values())
for path in SUITE.rglob('RESULT.json'):
    record = read(path)
    if 'final_fit' in record:
        fit_status.extend(record['final_fit'][a]['success'] for a in ('C2', 'phi'))
assert fit_status and all(fit_status)
legacy = read(ROOT / 'legacy_source_hashes.json')
assert all(sha(ROOT / item['path']) == item['sha256'] for item in legacy.values())
report = dict(status='PASS', synthetic=True, smoke_checks=len(checks),
              oracle_correlation=marker['oracle_correlation'], null_z_mean=marker['null_z_mean'],
              null_z_details=checks['null_shuffle_z_centering']['details'],
              converged_fits=len(fit_status), all_fits_converged=all(fit_status),
              standalone_export_identical=identical_export, snapshot_sources_match_manifest=True,
              normal_cli_and_validated_completed_run='PASS',
              threads=result['guard']['threads'], memory_cap_bytes=result['guard']['memory_cap_bytes'],
              learning_curve_points=len(result['learning_curves']['rows']),
              phi_parameters=len(next(iter(read(SUITE / 'phi.json')['coefficients'].values()))),
              D3_modes={m: len(result['simulation'][m]) for m in ('null', 'oracle', 'power')},
              result_path=str(SUITE / 'RESULT.json'), result_sha256=marker['result_sha256'])
(ROOT / 'V10_VERIFICATION.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
