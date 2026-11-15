"""Fetch an explicitly named public GEO file; no implicit large-data download."""
import argparse
import hashlib
from pathlib import Path
import re
from urllib.parse import quote
from urllib.request import Request, urlopen


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--accession', required=True)
    p.add_argument('--category', choices=['matrix', 'suppl', 'soft'], required=True)
    p.add_argument('--filename', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--sha256')
    a = p.parse_args()
    if not re.fullmatch(r'GSE\d+', a.accession):
        p.error('Expected a GEO series accession, e.g. GSE246622')
    if Path(a.filename).name != a.filename or a.filename in ('.', '..'):
        p.error('Filename must be a single deposited basename')
    family = a.accession[:-3] + 'nnn'
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{family}/{a.accession}/{a.category}/{quote(a.filename)}'
    out = Path(a.output_dir); out.mkdir(parents=True, exist_ok=True)
    target = out / a.filename
    if target.exists():
        raise FileExistsError('Existing input is protected: ' + str(target))
    temporary = target.with_name(target.name + '.partial')
    digest = hashlib.sha256()
    try:
        with urlopen(Request(url, headers={'User-Agent': 'CRRM388/0.2.0'}), timeout=120) as response, temporary.open('xb') as f:
            while block := response.read(1024 * 1024):
                digest.update(block); f.write(block)
        if a.sha256 and digest.hexdigest() != a.sha256.lower():
            raise ValueError('Downloaded file differs from the supplied source SHA256')
        temporary.rename(target)
    finally:
        temporary.unlink(missing_ok=True)
    print(digest.hexdigest(), target)


if __name__ == '__main__':
    main()
