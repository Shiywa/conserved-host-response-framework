#!/usr/bin/env python3
import csv,gzip,hashlib,json,re,shutil,urllib.request,concurrent.futures
from pathlib import Path
R=Path(str(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'phase06_external_validation_figure1'))
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
allrows=[];chosen=[];raw={}
for g in ['GSE270463','GSE309353']:
 text=gzip.decompress((R/'evidence/public_metadata'/f'{g}_family.soft.gz').read_bytes()).decode()
 samples=[];current=None
 for line in text.splitlines():
  if line.startswith('^SAMPLE = '):
   current={'gsm':line.split(' = ',1)[1]};samples.append(current)
  elif current is not None and line.startswith('!Sample_') and ' = ' in line:
   k,v=line.split(' = ',1);current.setdefault(k,[]).append(v)
 raw[g]=samples
 for s in samples:
  title=s['!Sample_title'][0]
  if g=='GSE270463':
   include=title.startswith('Donor ');donor=re.search(r'Donor\s+(\d+)',title);rep='1';condition='mock' if 'Mock' in title else 'RSV' if include else 'other';time='72'
  else:
   include=bool(re.search(r',\s*(mock|wt),',title,re.I));donor=re.search(r'Donor\s+(\d+)',title);m=re.search(r'rep\s*(\d+)',title,re.I);rep=m.group(1) if m else '';condition='mock' if ', mock,' in title.lower() else 'RSV' if ', wt,' in title.lower() else 'mutant';time='20'
  urls=[v for k,values in s.items() if k.startswith('!Sample_supplementary_file') for v in values]
  row=dict(dataset=g,gsm=s['gsm'],title=title,donor_id=donor.group(1) if donor else '',condition=condition,culture_replicate=rep,time_h=time,
           included_primary=include,exclusion_reason='' if include else 'CUTandRUN assay' if g=='GSE270463' else 'NS1-mutant arm',
           biosample=';'.join(v for v in s.get('!Sample_relation',[]) if 'BioSample' in v),supplementary_url=';'.join(urls),local_file='')
  if include:
   assert len(urls)==1,(g,title,urls)
   url=urls[0].replace('ftp://','https://');filename=url.rsplit('/',1)[-1];dest=R/'data/downloads'/g/filename;dest.parent.mkdir(exist_ok=True)
   row['local_file']=str(dest);chosen.append((url,dest))
  allrows.append(row)
save(R/'evidence/public_metadata/parsed_samples.json',raw)
with (R/'tables/01_external_sample_manifest.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
for g,n in [('GSE270463',6),('GSE309353',12)]:assert sum(x['dataset']==g and x['included_primary'] for x in allrows)==n
def fetch(job):
 url,p=job
 try:
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'CRRM-public-data-retrieval/1.0'}),timeout=45) as h:data=h.read();status=h.status
  p.write_bytes(data)
  if p.suffix=='.gz':gzip.decompress(data)
  d=dict(url=url,file=str(p),bytes=len(data),sha256=sha(p),status=status,gzip_integrity=True)
 except Exception as e:d=dict(url=url,file=str(p),error=repr(e))
 save(p.with_name(p.name+'.receipt.json'),d);return d
results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
 for out in ex.map(fetch,chosen):results.append(out);print(json.dumps(out,ensure_ascii=False))
save(R/'logs/count_download_summary.json',results)
shutil.copy2(Path(__file__),R/'scripts/01_manifest_and_download_counts.py')
assert all('error' not in x for x in results)
