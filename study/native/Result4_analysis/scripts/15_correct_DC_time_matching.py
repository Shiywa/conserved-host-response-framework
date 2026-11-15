from pathlib import Path
from datetime import datetime,timezone
import json,subprocess,shutil
import numpy as np,pandas as pd
OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result4_analysis');PY=__import__('sys').executable
receipt=OUT/'audit/DC_time_matching_correction.json'
if not receipt.exists():
    old=pd.read_csv(OUT/'results/primary_cellular_effects.csv')
    old.to_csv(OUT/'audit/primary_effects_before_DC_time_matching.csv',index=False)
    shutil.copy2(OUT/'results/primary_unit_effects.csv',OUT/'audit/primary_units_before_DC_time_matching.csv')
    receipt.write_text(json.dumps(dict(status='CORRECTION_DECLARED_AFTER_PRIMARY_RESULTS',declared_at=datetime.now(timezone.utc).isoformat(),trigger='Final timepoint audit found Untreated0h in CD1c_DC, whereas RSV is observed only at4/8/18h; marginal within-donor means mix unmatched baseline/time profiles',correction='For both paired sources, subtract within the same donor AND recorded time before equally averaging matched-time contrasts within each donor. Do not include0h as an RSV control; retain>=30cells and all eight original comparisons/families.',membership_changed=False,original_frozen_plan_preserved=True,reference_scaling_changed=False),indent=2)+'\n')
for script in ['04_score_and_test_profiles.py','05_cellular_interpretation.py','07_measurement_robustness.py','08_finalize_cellular_evidence.py']:
    subprocess.run([PY,str(OUT/'scripts'/script)],check=True)
subprocess.run(['/usr/local/bin/Rscript',str(OUT/'scripts/06_independent_statistics.R')],check=True)
old=pd.read_csv(OUT/'audit/primary_effects_before_DC_time_matching.csv');new=pd.read_csv(OUT/'results/primary_cellular_effects.csv')
cols=['comparison_id','score','case_units','control_units','estimate','lower','upper','P','resampling_P']
pd.testing.assert_frame_equal(old[~old.comparison_id.str.startswith('DC_')][cols].reset_index(drop=True),new[~new.comparison_id.str.startswith('DC_')][cols].reset_index(drop=True))
d=json.loads(receipt.read_text());d.update(status='PASS',completed_at=datetime.now(timezone.utc).isoformat(),non_DC_estimates_intervals_and_P_unchanged=True,independent_R_replay_passed=True,corrected_total_results=new[new.comparison_id.str.startswith('DC_')&new.score.eq('CRRM388')][['comparison_id','estimate','lower','upper','P','adjusted_P','resampling_P']].to_dict('records'))
receipt.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d,indent=2))
