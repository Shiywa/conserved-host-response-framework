from pathlib import Path
import requests, hashlib, json, concurrent.futures, datetime

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'phase10_functional_comparison_20261112')
SOURCES={
 'Wang2026_main.pdf':'https://spiral.imperial.ac.uk/bitstreams/6875af4b-5b83-4d49-a216-78c9323df6c7/download',
 'Wang2026_publisher.html':'https://www.journalofinfection.com/article/S0163-4453(26)00021-6/fulltext',
 'Wang2026_mmc1.docx':'https://ars.els-cdn.com/content/image/1-s2.0-S0163445326000216-mmc1.docx',
 'Wang2026_mmc1.pdf':'https://ars.els-cdn.com/content/image/1-s2.0-S0163445326000216-mmc1.pdf',
 'Wang2026_mmc1.xlsx':'https://ars.els-cdn.com/content/image/1-s2.0-S0163445326000216-mmc1.xlsx',
 'hallmark_2026_1.gmt':'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt',
 'Reactome_pathways.gmt.zip':'https://reactome.org/download/current/ReactomePathways.gmt.zip',
 'Reactome_release.json':'https://reactome.org/ContentService/data/database/version',
}
def get(item):
 name,url=item;p=OUT/'data/public'/name
 if p.exists():
  b=p.read_bytes();return dict(name=name,url=url,status='cached',bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
 try:
  r=requests.get(url,timeout=45)
  meta=dict(name=name,url=url,http_status=r.status_code,content_type=r.headers.get('Content-Type',''),bytes=len(r.content),downloaded_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
  valid=r.status_code==200 and (name.endswith('.html') or (name.endswith('.pdf') and r.content.startswith(b'%PDF')) or (name.endswith(('.docx','.xlsx','.zip')) and r.content.startswith(b'PK')) or (name.endswith('.gmt') and b'HALLMARK_' in r.content[:1000]) or name.endswith('.json'))
  if valid:p.write_bytes(r.content);meta['sha256']=hashlib.sha256(r.content).hexdigest();meta['status']='saved'
  else:meta['status']='not_retrieved';meta['response_start']=r.text[:180]
  return meta
 except Exception as e:return dict(name=name,url=url,status='error',error=str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:results=list(pool.map(get,SOURCES.items()))
(OUT/'evidence/public_source_downloads.json').write_text(json.dumps(results,ensure_ascii=False,indent=2))
for r in results:print(json.dumps(r,ensure_ascii=False))
