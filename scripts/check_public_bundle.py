"""Validate tracked release assets, portable sources, and public-file hygiene."""
import ast
import csv
import gzip
import hashlib
import io
from pathlib import Path
import re
import subprocess
import zipfile


def main():
    root = Path(__file__).resolve().parents[1]
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    with (root / 'study/configuration_manifest.tsv').open() as f:
        assets = list(csv.DictReader(f, delimiter='\t'))
    for row in assets:
        assert sha(root / row['repository_path']) == row['sha256'], row['logical_path']
    with (root / 'study/native_code_manifest.tsv').open() as f:
        native = list(csv.DictReader(f, delimiter='\t'))
    for row in native:
        p = root / row['path']
        assert sha(p) == row['portable_sha256'], row['path']
        if p.suffix == '.py':
            ast.parse(p.read_text())
    paths = subprocess.check_output(['git', 'ls-files', '-z'], cwd=root).decode().split('\0')
    personal_path = re.compile(rb'/(?:home|Users|Volumes)/[^\s\x00\"\'<>]+')
    secret = re.compile(rb'(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{50,}|-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----)')
    count = 0
    for name in paths:
        if not name:
            continue
        p = root / name
        assert p.suffix.lower() not in ('.doc', '.docx', '.pem', '.key', '.h5ad', '.bam', '.sra'), name
        assert not any(x in p.parts for x in ('__pycache__', '.venv', '.pytest_cache')), name
        assert p.stat().st_size < 10 * 1024 * 1024, name
        raw = p.read_bytes(); texts = [raw]
        if p.suffix == '.gz':
            texts = [gzip.decompress(raw)]
        if p.suffix == '.xlsx':
            with zipfile.ZipFile(io.BytesIO(raw)) as z:
                texts = [z.read(n) for n in z.namelist() if n.endswith(('.xml', '.rels'))]
        assert not any(personal_path.search(s) or secret.search(s) for s in texts), name
        count += 1
    forbidden_current=[p for p in (root/'study').rglob('*') if p.is_file() and p.suffix in ('.gz','.npz','.xlsx','.h5ad')]
    assert not forbidden_current, 'Current study release contains code/configuration only'
    print(f'PASS: {count} tracked files; {len(assets)} configuration aliases; {len(native)} native scripts; no new matrices/complete result bundles')


if __name__ == '__main__':
    main()
