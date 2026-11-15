"""Portable entry points for the current paper's curated native workflows."""
from __future__ import annotations
from pathlib import Path
import gzip,hashlib,json,os,shutil,subprocess,sys
import pandas as pd

CLINICAL="clinical_otherwisehealthy_262_r1"
PLOTS=[("phase16_identifier12pt_20261112","01_build_figures.py"),
       ("phase18_Result2_external_response_20261112","05_plot_figures.py"),
       (f"{CLINICAL}/main","07_plot_clinical_figures.py"),
       (f"{CLINICAL}/incremental","07_plot_exploratory_figure.py"),
       ("Result4_analysis","20_plot_cellular_focus_r1.py"),
       ("Result5_analysis","15_plot_crossvirus_focus_r2.py"),
       (f"{CLINICAL}/host","17_plot_component_focus_r3.py")]

def study_directory(explicit=None):
    candidates=[Path(explicit)] if explicit else []
    if os.environ.get("CRRM_STUDY_DIR"):candidates.append(Path(os.environ["CRRM_STUDY_DIR"]))
    candidates += [p/"study" for p in [Path(__file__).resolve().parents[2],*Path.cwd().resolve().parents,Path.cwd()]]
    for p in candidates:
        if (p/"native_code_manifest.tsv").is_file():return p.resolve()
    raise FileNotFoundError("Locate the source-release study directory using --study-dir or CRRM_STUDY_DIR")

def input_study_directory(explicit=None):
    """Locate separately held processed inputs; these are not in the code release."""
    value=explicit or os.environ.get("CRRM_INPUT_BUNDLE")
    if not value:raise FileNotFoundError("Processed inputs are held separately. Supply --input-bundle or CRRM_INPUT_BUNDLE for data-dependent replay; see study/WORKFLOWS.md")
    p=Path(value).resolve()
    for candidate in [p/"study",p]:
        if (candidate/"asset_manifest.tsv").is_file():return candidate
    raise FileNotFoundError("The local input bundle must contain study/asset_manifest.tsv and its separately held files")

def study_input_path(relative,*,input_bundle=None):
    return input_study_directory(input_bundle)/relative

def prepare(workspace, *, study_dir=None, data_root=None, resource_root=None,input_bundle=None):
    bundle=study_directory(study_dir);root=Path(workspace).resolve();root.mkdir(parents=True,exist_ok=True)
    data=Path(data_root).resolve() if data_root else root/"external_data"
    resource=Path(resource_root).resolve() if resource_root else root/"external_resource"
    marker=root/"PORTABLE_STUDY.json"
    if marker.exists():raise FileExistsError("Use a new workspace; an executed study workspace is never silently replaced")
    supplied=input_bundle or os.environ.get("CRRM_INPUT_BUNDLE")
    source_bundle=input_study_directory(supplied) if supplied else None
    manifest=pd.read_csv(source_bundle/"asset_manifest.tsv",sep="\t") if source_bundle else pd.DataFrame(columns=["logical_path","source_sha256"])
    for row in manifest.itertuples():
        source=source_bundle.parent/row.repository_path;raw=source.read_bytes()
        if hashlib.sha256(raw).hexdigest()!=row.stored_sha256:raise ValueError("Asset checksum mismatch: "+row.logical_path)
        content=gzip.decompress(raw) if row.compression=="gzip" else raw
        if hashlib.sha256(content).hexdigest()!=row.sha256:raise ValueError("Decoded checksum mismatch: "+row.logical_path)
        if row.compression=="gzip":
            text=content.decode().replace("@workspace",str(root)).replace("@data",str(data)).replace("@resource",str(resource))
            content=text.encode()
        dest=root/row.logical_path;dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.name.endswith('.gz'):content=gzip.compress(content,mtime=0)
        dest.write_bytes(content)
    configuration=pd.read_csv(bundle/"configuration_manifest.tsv",sep="\t")
    for row in configuration.itertuples():
        dest=root/row.logical_path
        if dest.exists():continue
        source=bundle.parent/row.repository_path
        if hashlib.sha256(source.read_bytes()).hexdigest()!=row.sha256:raise ValueError("Configuration checksum mismatch: "+row.logical_path)
        content=source.read_text().replace("@workspace",str(root)).replace("@data",str(data)).replace("@resource",str(resource))
        dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(content)
    selection=root/CLINICAL/"source_data/GSE246622_included_262.tsv"
    if not selection.exists():
        selection.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(bundle/"selection/GSE246622_included_262.tsv",selection)
    for source in (bundle/"native").rglob('*'):
        if source.is_file():
            dest=root/source.relative_to(bundle/"native");dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
    # Normalize the fingerprints of explicitly packaged transport artifacts.
    # Unbundled external raw inputs keep their original checksums.
    hashes={}
    for r in manifest.itertuples():hashes[r.source_sha256]=hashlib.sha256((root/r.logical_path).read_bytes()).hexdigest()
    for r in configuration.itertuples():hashes.setdefault(r.source_sha256,hashlib.sha256((root/r.logical_path).read_bytes()).hexdigest())
    native=pd.read_csv(bundle/'native_code_manifest.tsv',sep='\t')
    for r in native.itertuples():hashes[r.source_sha256]=r.portable_sha256
    for _ in range(3):
        for p in root.rglob('parameters/*.json'):
            old=p.read_bytes()
            def transform(x):
                if isinstance(x,dict):return {k:transform(v) for k,v in x.items()}
                if isinstance(x,list):return [transform(v) for v in x]
                return hashes.get(x,x) if isinstance(x,str) else x
            try:value=transform(json.loads(old))
            except (ValueError,UnicodeError):continue
            p.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
            oldhash=hashlib.sha256(old).hexdigest();newhash=hashlib.sha256(p.read_bytes()).hexdigest()
            for k,v in list(hashes.items()):
                if v==oldhash:hashes[k]=newhash
    for stage in set(manifest.logical_path.astype(str).str.split('/').str[0])|set(configuration.logical_path.str.split('/').str[0])|{CLINICAL}:
        for sub in ["data","figures","tables","models","parameters","logs","audit","reports"]:(root/stage/sub).mkdir(parents=True,exist_ok=True)
    for part in ["main","incremental","host"]:
        for sub in ["data/processed","figures","tables","models","parameters","logs","audit","reports"]:(root/CLINICAL/part/sub).mkdir(parents=True,exist_ok=True)
    # New transport receipts contain unchanged scientific settings and the supplied
    # input hashes. They do not claim a new prospective analysis declaration.
    p=root/CLINICAL/"parameters/population_revision_plan.json"
    value=json.loads(p.read_text());value["portable_transport"]="Paths normalized for the published study; original scientific settings retained"
    p.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
    (p.parent/"population_revision_receipt.json").write_text(json.dumps({"plan_sha256":hashlib.sha256(p.read_bytes()).hexdigest(),"new_analysis_started":False})+'\n')
    marker.write_text(json.dumps(dict(bundle=str(bundle),data_root=str(data),resource_root=str(resource),
                                    population="otherwise-healthy262",processed_inputs_loaded=bool(source_bundle),input_manifest_sha256=hashlib.sha256((source_bundle/"asset_manifest.tsv").read_bytes()).hexdigest() if source_bundle else None),indent=2)+'\n')
    return root

def run_native(workspace,stage,script,*,rscript="Rscript",arguments=()):
    root=Path(workspace).resolve();meta=json.loads((root/"PORTABLE_STUDY.json").read_text())
    p=root/stage/"scripts"/script
    if not p.is_file():raise FileNotFoundError(p)
    env=os.environ.copy();env.update(CRRM_WORKSPACE=str(root),CRRM_DATA_ROOT=meta['data_root'],CRRM_RESOURCE_ROOT=meta['resource_root'],
      OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MPLCONFIGDIR=str(root/"matplotlib_cache"))
    command=[rscript if p.suffix.lower()=='.r' else sys.executable,str(p),*map(str,arguments)]
    subprocess.run(command,env=env,cwd=root,check=True)

def clinical(workspace,*,regenerate_splines=True,rscript="Rscript"):
    if regenerate_splines:run_native(workspace,CLINICAL,"01_generate_spline_designs.R",rscript=rscript)
    for name in ["02_rebuild_clinical.py","04_rebuild_incremental.py","05_update_clinical_bridge.py"]:
        run_native(workspace,CLINICAL,name,rscript=rscript)
    # Clinical bridge affects the common result asset consumed by later scripts.
    root=Path(workspace)
    for name in ["clinical_gene_alignment.csv","clinical_gene_alignment_points.csv"]:
        shutil.copy2(root/CLINICAL/"host/results"/name,root/"Result6_analysis/results"/name)

def figures(workspace,*,rscript="Rscript"):
    for stage,script in PLOTS:run_native(workspace,stage,script,rscript=rscript)
    return [p for p in Path(workspace).rglob('figures/*.pdf')]

def export_tables(workspace,output,*,study_dir=None):
    root=Path(workspace);out=Path(output);out.mkdir(parents=True,exist_ok=True)
    sources=sorted((root/'tables').glob('TableS*.xlsx'))
    if not sources:raise FileNotFoundError("Generate or supply the local supplementary workbooks first; they are not part of the code release")
    rows=[]
    for source in sources:
        dest=out/source.name;shutil.copy2(source,dest)
        rows.append(dict(table=dest.name,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),role='Traceable export of an existing local table;statistics are not recomputed'))
    s=study_directory(study_dir)/'tables'
    for p in s.iterdir():
        if p.is_file():shutil.copy2(p,out/p.name)
    (out/'EXPORT_MANIFEST.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows
