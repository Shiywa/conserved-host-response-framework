from pathlib import Path
import csv,json,hashlib
OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'phase10_functional_comparison_20261112')
cfg=json.loads((OUT/'parameters/functional_comparison_plan_v1.json').read_text())
rows=list(csv.DictReader(open(cfg['inputs']['universe']['path'],encoding='utf-8-sig')));un={x['display_symbol']:x for x in rows}
members=list(csv.DictReader(open(cfg['inputs']['priority388']['path'],encoding='utf-8-sig')));ids={x['human_ensembl'] for x in members}
def write(name,rows):
 with (OUT/'tables'/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
reported=[
 ('RTP4','Blood shared','Complete four-gene blood intersection','p2, blood Results'),
 ('BST1','Blood shared','Complete four-gene blood intersection','p2, blood Results'),
 ('STAT1','Blood shared','Complete four-gene blood intersection','p2, blood Results'),
 ('TMEM40','Blood shared','Complete four-gene blood intersection','p2, blood Results'),
 ('USP18','Adult cross-tissue','Named gene among adult blood/nasal/lung overlap','p3, Fig4B-C; ISGF3 is a complex and is not assigned to an invented gene'),
 ('IRF9','Adult cross-tissue','Named gene among adult blood/nasal/lung overlap','p3, Fig4B-C; ISGF3 is a complex and is not assigned to an invented gene'),
 ('CEBPB','IL-17 / S100','Named selected IL17-associated gene, not a full author signature','p5, Fig4J discussion'),
 ('CXCL10','IL-17 / S100','Named selected IL17-associated gene, not a full author signature','p5, Fig4J discussion'),
 ('IL1B','IL-17 / S100','Named selected IL17-associated gene, not a full author signature','p5, Fig4J discussion'),
 ('S100A8','IL-17 / S100','Named selected IL17/S100 effector, not a full author signature','p5, Fig4J and experimental validation'),
 ('S100A9','IL-17 / S100','Named selected IL17/S100 effector, not a full author signature','p5, Fig4J and experimental validation'),
]
out=[]
for gene,group,scope,source in reported:
 r=un.get(gene);assert r is not None,gene
 x=dict(display_symbol=gene,Wang2026_context=group,comparison_scope=scope,published_direction='up',paper_location=source,source_DOI='10.1016/j.jinf.2026.106696',human_ensembl=r['human_ensembl'],CRRM_member=r['human_ensembl'] in ids)
 for fam in ['A549','mouse','cotton_rat']:
  for key in ['mapped_and_tested','log2FC','FDR']:x[fam+'_'+key]=r[fam+'_'+key]
  try:x[fam+'_qualifying_up']=r[fam+'_mapped_and_tested'].lower() in ['true','t','1'] and float(r[fam+'_FDR'])<.05 and float(r[fam+'_log2FC'])>.58
  except ValueError:x[fam+'_qualifying_up']=False
 out.append(x)
write('10_Wang2026_reported_gene_comparison.csv',out)
design=[
 ('Scientific objective','Compare children, adult challenge and mice; identify disease-associated therapeutic targets','Define a fixed cross-model RSV host-response program for later application'),
 ('Discovery sources','Human adult and child blood/nasal; adult lung; mouse blood/nasal/lung','A549 GSE231788; mouse lung GSE197868; cotton-rat lung GSE297903'),
 ('Gene evidence','limma; P<0.05 and FDR<0.1; no nonzero effect-size cutoff; gene-symbol matching','Native edgeR QL, FDR<0.05 and abs(log2FC)>0.58 in every one of3 families; unique orthologue projection'),
 ('Direction','DE includes up and down; shared-pathway and intersection analyses emphasize upregulated genes','Both directions explicitly retained:280up and108down; signed score'),
 ('Central reported biology','IFN/cytokines,TLR,NFkB; primary-infection IL17/S100A8/A9; drug/intervention evidence','Function results quantified independently from frozen388; no intervention evidence attributed to this gene set'),
 ('Reusable scoring object','No single complete frozen cross-system scoring list or scoring formula identified in published main text and downloaded supplement','Fixed388 members and direction multipliers; equal-weight signed z score'),
 ('Human primary airway checking','Clinical tissue comparisons and adult challenge are constituent parts of study','Three separate primary-human studies after membership freeze; 20,72,96h are different study/model conditions'),
 ('Public gene overlap scope','Four named blood-intersection genes; selected named IL17/S100 and adult cross-tissue genes are directly auditable','Exact scoped overlap and discovery effect evidence; no claim of full40/49/30-gene-set overlap'),
]
write('11_Wang2026_study_design_comparison.csv',[dict(dimension=a,Wang2026=b,CRRM388=c) for a,b,c in design])
cohorts=[
 ('Human adult','Blood','Day0 N9; Day7 N9','9 symptomatic volunteers per methods','GSE17156','Public accession in article','Figure1 legend instead gives baseline11, peak9; exact author subset not reconstructed'),
 ('Human adult','Nasal scrape','Day0 N9; Day7 N9','9 donors','GSE155237 (Day0 only as specified)','Day7 described as unpublished in Table1','Do not assume public baseline series contains author Day7 comparison'),
 ('Human adult','Lung biopsy','Day0 N5; Day7 N5','5 donors','TBC in availability statement','Unpublished per Table1','No study accession identified in article'),
 ('Human child','Blood','Control45; RSV42','Female subset in Results','E-MTAB-12793 in Methods','Availability statement still TBC','Declared sample counts copied from Table1; no sample-level reanalysis'),
 ('Human child','Nasal mucosa','Control22; RSV24','Separate donors, natural infection','phs003053.v1.p1','dbGaP resource cited','Restricted individual-level data not downloaded'),
 ('Mouse','Blood','Control5; RSV5','6-7week females; d7','38382577 printed as ArrayExpress locator','Not a valid ArrayExpress-style accession','Figure1 legend gives4+4; source discrepancy retained'),
 ('Mouse','Nasal mucosa','Control5; RSV5','6-7week females; d7','Unpublished/current study in Table1','No standalone accession identified','No assumed public count matrix'),
 ('Mouse','Lung','Control5; RSV5','6-7week females; d7','38382577 printed as ArrayExpress locator','Not a valid ArrayExpress-style accession','Article cites Wang2024 Mucosal Immunology for provenance'),
]
write('12_Wang2026_cohort_and_overlap_audit.csv',[dict(species_age=a,tissue=b,Table1_declared_samples=c,design=d,reported_locator=e,public_scope=f,note=g,overlap_with_current_three_discovery_accessions=False,overlap_with_current_three_external_accessions=False) for a,b,c,d,e,f,g in cohorts])
evidence=dict(paper_DOI='10.1016/j.jinf.2026.106696',main_pdf='data/public/Wang2026_main.pdf',supplement='data/public/Wang2026_mmc1.pdf',supplement_contents=['FigureS1 six-way orthologue overlap','FigureS2 flow cytometry gating'],complete_shared_gene_lists_not_found=[dict(n=40,scope='three-system nasal up intersection'),dict(n=49,scope='mouse blood/nasal/lung up intersection'),dict(n=30,scope='mouse/adult lung up intersection')],complete_scoring_signature_found=False,comparison_panel_genes=len(out),scoped_complete_blood_set=dict(n_author_genes=4,n_CRRM_members=sum(x['CRRM_member'] for x in out if x['Wang2026_context']=='Blood shared')),representative_genes_not_full_signature=True,ISGF3_handling='Author reports an ISGF3 complex; not treated as a single unambiguous human gene or silently replaced by IRF9/STAT1/STAT2',sample_sizes='Table1 values and differing Figure1 legend values both retained; no invented exact analysis sample manifest',no_full_signature_performance_benchmark=True)
(OUT/'Wang2026_comparison_summary.json').write_text(json.dumps(evidence,ensure_ascii=False,indent=2))
print(json.dumps(evidence,ensure_ascii=False,indent=2));print([(x['display_symbol'],x['CRRM_member']) for x in out])
