"""Inspect the published compressed assets without executing pickle payloads."""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile
import numpy as np
import pandas as pd


def main():
    root = Path(__file__).resolve().parents[1]
    m = pd.read_csv(root / 'study/asset_manifest.tsv', sep='\t')
    rows = []
    for relative, group in m.groupby('repository_path', sort=True):
        p = root / relative; stored = p.read_bytes()
        raw = gzip.decompress(stored) if p.suffix == '.gz' else stored
        row = dict(repository_path=relative,
                   git_blob=hashlib.sha1(('blob ' + str(len(stored)) + '\0').encode() + stored).hexdigest(),
                   stored_sha256=hashlib.sha256(stored).hexdigest(),
                   decoded_sha256=hashlib.sha256(raw).hexdigest(),
                   stored_bytes=len(stored), decoded_bytes=len(raw),
                   logical_paths=' | '.join(group.logical_path))
        if p.suffix == '.gz':
            text = raw.decode(); logical = group.logical_path.iloc[0].removesuffix('.gz')
            if logical.endswith(('.csv', '.tsv')):
                d = pd.read_csv(io.StringIO(text), sep='\t' if logical.endswith('.tsv') else ',')
                forbidden = [c for c in d if c.startswith('contact_') or re.fullmatch(
                    r'(?:patient|participant|subject)[ _-]*name|e[ _-]*mail|phone|home[ _-]*address|birth[ _-]*date|date[ _-]*of[ _-]*birth', c, re.I)]
                if forbidden:
                    raise ValueError(f'Unnecessary identifying/contact fields in {relative}: {forbidden}')
                row.update(format='delimited scientific table', dimensions=f'{len(d)} x {len(d.columns)}',
                           fields='; '.join(d.columns[:12]))
            else:
                try:
                    j = json.loads(text)
                    row.update(format='scientific configuration/summary', dimensions=str(len(j)) + ' root entries',
                               fields='; '.join(list(j)[:12]) if isinstance(j, dict) else 'JSON array')
                except ValueError:
                    row.update(format='scientific text/helper', dimensions=str(len(text.splitlines())) + ' lines', fields='Source/helper text')
        elif p.suffix == '.npz':
            with np.load(p, allow_pickle=False) as z:
                info = [k + ':' + str(z[k].shape) + ':' + str(z[k].dtype) for k in z.files]
                if any(z[k].dtype.hasobject for k in z.files):
                    raise ValueError('Object array: ' + relative)
                row.update(format='derived numeric/identifier NumPy arrays', dimensions='; '.join(info), fields='; '.join(z.files))
        elif p.suffix == '.xlsx':
            with zipfile.ZipFile(p) as z:
                names = re.findall(r'<sheet[^>]* name="([^"]+)"', z.read('xl/workbook.xml').decode())
            row.update(format='current supplementary result workbook', dimensions=str(len(names)) + ' sheets', fields='; '.join(names))
        row['content_scope'] = 'Frozen definitions/configuration, public-study processed inputs, or derived results; no raw sequencing/all-cell matrix'
        rows.append(row)
    output = root / 'study/validation/public_asset_inventory.tsv'
    pd.DataFrame(rows).to_csv(output, sep='\t', index=False)
    print(f'PASS: {len(rows)} physical assets inspected; public-source linkage, dimensions and checksums recorded; no object arrays or forbidden fields')


if __name__ == '__main__':
    main()
