"""Download original design records; expression inputs remain local and unchanged."""
from pathlib import Path
from urllib.request import Request, urlopen
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib, json

OUT = (Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result6_analysis')
URLS = {
    'GSE179353_family.soft.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179353/soft/GSE179353_family.soft.gz',
    'GSE291797_family.soft.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291797/soft/GSE291797_family.soft.gz',
    'SRP326693_ENA_samples.tsv': 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP326693&result=read_run&fields=run_accession,experiment_accession,sample_accession,secondary_sample_accession,experiment_title,sample_title,library_name&format=tsv',
    'SRP326693_SRA.xml': 'https://www.ebi.ac.uk/ena/browser/api/xml/SRP326693?download=false&includeLinks=true',
    'BRD4_article.html': 'https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2021.728661/full',
    'BRD4_PubMed.xml': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=34765643&retmode=xml',
}

def fetch(item):
    name, url = item
    dest = OUT / 'data' / name
    try:
        if dest.exists() and dest.stat().st_size:
            data = dest.read_bytes()
        else:
            with urlopen(Request(url, headers={'User-Agent': 'Mozilla/5.0'}), timeout=40) as response:
                data = response.read()
            dest.write_bytes(data)
        return dict(file=name, url=url, status='OK', bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
    except Exception as exc:
        return dict(file=name, url=url, status='ERROR', reason=str(exc))

if __name__ == '__main__':
    with ThreadPoolExecutor(5) as pool:
        records = list(pool.map(fetch, URLS.items()))
    (OUT / 'audit' / 'metadata_fetch.json').write_text(json.dumps(dict(completed_at=datetime.now(timezone.utc).isoformat(), records=records), indent=2) + '\n')
    print(json.dumps(records, indent=2))
