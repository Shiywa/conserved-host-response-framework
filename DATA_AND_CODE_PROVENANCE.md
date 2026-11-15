# Code and input provenance

The public release contains code, frozen definitions, readable parameters,
data/sample manifests, software-environment records and validation summaries.
Current expression matrices, large cellular objects, complete result tables
and supplementary workbooks are not distributed in the code repository.

`study/native_code_manifest.tsv` records original/portable script checksums and
adaptations. Runtime roots replace local paths; author-document backup and
baseline bookkeeping are omitted. Public GEO submitter contact fields are
omitted from two source-table exports; all other cell strings and scientific
calculations are unchanged. The clinical raw-preparation selection uses the
current 262-profile ledger. The retired clinical driver is retained only as
statistical functions actually reused by the current code.

`study/configuration_manifest.tsv` maps stage configurations; fixed gene lists
contain no sample expression. `study/input_requirements.tsv` lists paths and
checksums from the recorded local processed-input replay. It is an input
contract, not a claim that those files are distributed or hosted. Readers
obtain public inputs from the accession/source records and run the corresponding
preparation scripts. Author-side regression can use a separately held local
bundle through CRRM_INPUT_BUNDLE. No such bundle is uploaded with this release.

Discovery sources: GSE231788, GSE197868 and GSE297903. Primary-human evaluation:
GSE231666. GSE155237/GSE166161 form one external challenge family. Clinical
sources have distinct designs: GSE246622 otherwise-healthy selection, GSE69606
severity/recovery, and GSE105450 hospitalization context. The cellular sources
include GSE283746 (infant PBMCs) and GSE334417 (DC cultures). Cross-virus whole
blood uses GSE38900/GPL6884 and GSE42026. Perturbation sources include GSE179353
and GSE291797. GEO records use the stable accession viewer at
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=ACCESSION.

Public author CyTOF/clinical files are identified by
[RESCEU source commit c16315dbd28fcab727c70a10f238dbfbb2a0ce6e](https://github.com/JanssenDiscoveryOmics/RESCEU_ImmuneProfiling/tree/c16315dbd28fcab727c70a10f238dbfbb2a0ce6e).
Their original source Git blob hashes are fcd0cdef404a83ea39f652c3a645f3f32a56209c
(datTraits.txt) and 3876600dcb61a64ff258075098f06efac9050a50
(datTraits_counts.txt). Local source bytes matched both during audit; the
clinical records and large author objects remain outside this code release.

Raw and processed source data retain provider terms. MIT applies to software,
with the original copyright unchanged. No Word draft, private review-access
link, credential or font is part of the code release. The immutable v0.1.0
history retains its original publicly distributed example inputs.
