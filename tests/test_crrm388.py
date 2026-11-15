"""Scientific migration checks; separate from the v0.1.0 compatibility tests."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import numpy as np
import pandas as pd
from scipy import stats
from conserved_host_response.reference import FrozenReference,score,CoveragePolicy
from conserved_host_response.crrm388 import reconstruct,verify_definition
from conserved_host_response.statistics import adjust,welch,hc3,exact_signflip
from conserved_host_response.study import study_directory,prepare,clinical,CLINICAL,input_study_directory
from conserved_host_response.replay import verify

class SignedReferenceTests(unittest.TestCase):
    def setUp(self):
        self.genes=pd.DataFrame({'frozen_direction':['up','down']},index=['UP','DOWN'])
        self.x=pd.DataFrame({'a':[0.,2.],'b':[2.,0.]},index=self.genes.index)

    def test_down_members_are_reversed(self):
        ref=FrozenReference.fit(self.x,['a','b']);v,c,z=score(self.x,ref,self.genes)
        self.assertAlmostEqual(v.loc['a','CRRM388'],-1/np.sqrt(2))
        self.assertAlmostEqual(z.loc['DOWN','a'],-1/np.sqrt(2))
        self.assertTrue(np.allclose(v.CRRM388,v['Up contribution']+v['Down contribution']))

    def test_projecting_new_samples_never_refits_reference(self):
        ref=FrozenReference.fit(self.x,['a','b']);before=ref.fingerprint()
        extended=self.x.assign(recovery=[800.,-300.]);projected=ref.project(extended)
        self.assertEqual(before,ref.fingerprint())
        np.testing.assert_array_equal(projected[['a','b']],ref.project(self.x))
        with TemporaryDirectory() as d:
            p=Path(d)/'reference.json';ref.save(p);loaded=FrozenReference.load(p)
            self.assertEqual(ref.fingerprint(),loaded.fingerprint())

    def test_missing_direction_fails_eligibility(self):
        # Overall coverage is 80%: this must fail specifically because the
        # down direction is unmeasured, not because of the overall 70% gate.
        ids=[f'U{i}' for i in range(8)]+['D0','D1']
        genes=pd.DataFrame({'frozen_direction':['up']*8+['down']*2},index=ids)
        x=pd.DataFrame({'a':np.arange(8.),'b':np.arange(8.)+2},index=ids[:8])
        ref=FrozenReference.fit(x,['a','b']);v,c,_=score(x,ref,genes)
        self.assertTrue(v.CRRM388.isna().all())
        self.assertEqual(c.status.iloc[0],'UNESTIMABLE_COVERAGE')

    def test_invariant_members_are_not_zero_imputed(self):
        x=self.x.copy();x.loc['DOWN']=1.;ref=FrozenReference.fit(x,['a','b'])
        self.assertFalse(ref.parameters.loc['DOWN','eligible'])
        self.assertEqual(ref.project(x).index.tolist(),['UP'])

    def test_initial_human_sd_and_count_support_are_explicit(self):
        x=self.x*1e-9
        initial=FrozenReference.fit(x,['a','b'],sd_threshold=0.,nonzero_support=pd.Series({'UP':True,'DOWN':False}))
        standard=FrozenReference.fit(x,['a','b'])
        self.assertTrue(initial.parameters.loc['UP','eligible']);self.assertFalse(initial.parameters.loc['DOWN','eligible'])
        self.assertFalse(standard.parameters.eligible.any())

    def test_overlapping_annotations_do_not_duplicate_core_genes(self):
        ids=[f'G{i}' for i in range(12)];genes=pd.DataFrame({'frozen_direction':['up']*6+['down']*6},index=ids)
        x=pd.DataFrame(np.arange(48).reshape(12,4)**1.17,index=ids,columns=list('abcd'))
        ref=FrozenReference.fit(x,list('abcd'));v,c,z=score(x,ref,genes,{'A':ids[:8],'B':ids[4:]})
        np.testing.assert_allclose(v.CRRM388,z.mean(0))
        np.testing.assert_allclose(v.CRRM388,v['Up contribution']+v['Down contribution'])
        self.assertEqual(c.scored.iloc[0],12)

    def test_unit_equal_reference_uses_population_sd(self):
        x=pd.DataFrame({'a1':[0.],'a2':[2.],'a3':[4.],'b1':[8.]},index=['G'])
        units=pd.Series({'a1':'A','a2':'A','a3':'A','b1':'B'})
        ref=FrozenReference.fit(x,x.columns.tolist(),units=units)
        self.assertAlmostEqual(ref.parameters.loc['G','mean'],5.)
        self.assertAlmostEqual(ref.parameters.loc['G','sd'],np.sqrt(31/3))
        self.assertAlmostEqual(sum(ref.metadata['profile_weights'][s] for s in ['a1','a2','a3']),.5)

    def test_parent_parameters_are_inherited_on_omission(self):
        ref=FrozenReference.fit(self.x,['a','b']);before=ref.fingerprint()
        score(self.x,ref,self.genes,omitted=['UP']);self.assertEqual(before,ref.fingerprint())

    def test_duplicate_members_are_rejected(self):
        ref=FrozenReference.fit(self.x,['a','b'])
        with self.assertRaises(ValueError):score(self.x,ref,pd.concat([self.genes,self.genes.iloc[:1]]))

    def test_subset_minimum_is_enforced_even_with_complete_coverage(self):
        ids=[f'G{i}' for i in range(8)];genes=pd.DataFrame({'frozen_direction':['up']*4+['down']*4},index=ids)
        x=pd.DataFrame(np.arange(24).reshape(8,3),index=ids,columns=list('abc'))
        ref=FrozenReference.fit(x,list('abc'));v,c,_=score(x,ref,genes,{'four':ids[:4]})
        self.assertTrue(v.four.isna().all());self.assertFalse(c.set_index('set').loc['four','eligible'])

    def test_original_adjustment_denominator_keeps_unavailable_slots(self):
        q=adjust([.01,np.nan,.04],'BH',family_size=4)
        np.testing.assert_allclose(q[[0,2]],[.04,.08]);self.assertTrue(np.isnan(q[1]))
        np.testing.assert_allclose(adjust([.01,.04],'holm',family_size=4),[.04,.12])

    def test_current_welch_is_not_the_legacy_effect_size(self):
        a=np.array([1.,3.,5.,8.]);b=np.array([-1.,0.,2.])
        r=welch(a,b);tt=stats.ttest_ind(a,b,equal_var=False)
        self.assertAlmostEqual(r['effect'],a.mean()-b.mean());self.assertAlmostEqual(r['p_value'],tt.pvalue)
        self.assertEqual(exact_signflip([1.,2.,3.]),.25)


class FrozenDefinitionTests(unittest.TestCase):
    def test_current_frozen_definition_and_annotation_contract(self):
        bundle=study_directory()
        result=verify_definition(pd.read_csv(bundle/'frozen/crrm388_genes.tsv',sep='\t'),pd.read_csv(bundle/'frozen/crrm388_subsets.tsv',sep='\t'))
        self.assertEqual((result['members'],result['up'],result['down']),(388,280,108))
        self.assertEqual(result['nominal_subset_sizes']['IFN-list overlap'],79)
        self.assertEqual(result['nominal_subset_sizes']['Outside IFN lists'],309)
        self.assertEqual(result['nominal_subset_sizes']['Neutrophil degranulation'],25)


class CurrentStudyRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bundle=study_directory()
        try:cls.inputs=input_study_directory()
        except FileNotFoundError as e:raise unittest.SkipTest(str(e))
        cls.temp=TemporaryDirectory();cls.root=Path(cls.temp.name)/'study_run'
        prepare(cls.root,study_dir=cls.bundle)
        # The independent saved R basis is used here. CLI's default separately
        # regenerates it with R; the actual R run is recorded in study/environment.
        clinical(cls.root,regenerate_splines=False)
        cls.report=verify(cls.bundle,cls.root)

    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    def test_actual_discovery_and_primary_human_regressions(self):
        self.assertEqual(self.report['discovery_selected'],388)
        self.assertEqual(self.report['discovery_background'],8607)
        self.assertLess(self.report['primary_human_score_max_error'],1e-12)

    def test_real_otherwisehealthy_population_and_direction_coverage(self):
        self.assertEqual(self.report['whole_blood_profiles'],262)
        self.assertEqual((self.report['whole_blood_scored'],self.report['up'],self.report['down']),(369,267,102))
        ledger=pd.read_csv(self.inputs/'examples/GSE246622_included_262.tsv',sep='\t')
        scores=pd.read_csv(self.root/CLINICAL/'main/tables/07_sample_scores.csv')
        wb=scores[scores.dataset.eq('GSE246622')]
        self.assertEqual(set(wb.sample_id),set(ledger.laboratory_id))
        self.assertEqual(wb.sample_id.nunique(),262)
        self.assertEqual(wb[wb.group.isin(['Severe','Mild'])].group.value_counts().to_dict(),{'Mild':79,'Severe':36})

    def test_saved_clinical_effect_and_separate_P_families(self):
        r=self.report['whole_blood_total']
        self.assertAlmostEqual(r['effect'],.095171119063625,places=12)
        self.assertAlmostEqual(r['p_value'],.0233212655353988,places=12)
        self.assertAlmostEqual(r['primary_Holm_P'],.0466425310707977,places=12)
        self.assertEqual(r['resampling_P'],.0393);self.assertEqual(r['resampling_Holm_P'],.0786)

    def test_HC3_kernel_matches_independently_saved_clinical_model(self):
        base=self.root/CLINICAL/'main';d=pd.read_csv(base/'models/GSE246622_adjusted_design.csv')
        s=pd.read_csv(base/'tables/07_sample_scores.csv');s=s[s.dataset.eq('GSE246622')].set_index('sample_id')
        X=d.iloc[:,1:].to_numpy();contrast=np.zeros(X.shape[1]);contrast[d.columns[1:].tolist().index('high')]=1
        r=hc3(s.loc[d.sample_id,'score_Total'].to_numpy(),X,contrast)
        saved=pd.read_csv(self.inputs/'reference_results/17_covariate_adjusted_effects.csv');saved=saved[saved.dataset.eq('GSE246622')&saved.model.eq('adjusted')].iloc[0]
        for k in ['effect','CI95_low','CI95_high','p_value']:self.assertAlmostEqual(r[k],saved[k],places=10)

    def test_all_regenerated_numeric_tables_match_saved_results(self):
        self.assertGreaterEqual(self.report['clinical_tables_compared'],50)
        self.assertLess(self.report['clinical_max_numeric_error'],1e-10)

    def test_leave_one_reference_and_unchanged_recovery(self):
        self.assertEqual((self.report['leave_one_evaluation_n'],self.report['leave_one_training_n'],self.report['leave_one_reference_n']),(115,114,261))
        self.assertEqual(self.report['GSE69606_recovery_pairs'],15)
        prediction=pd.read_csv(self.root/CLINICAL/'incremental/tables/19_leave_one_unit_out_predictions.csv')
        metadata=pd.read_csv(self.root/CLINICAL/'main/tables/07_sample_scores.csv')
        ids=set(metadata.loc[metadata.dataset.eq('GSE246622')&metadata.group.isin(['Severe','Mild']),'sample_id'])
        for view in ['GSE246622','GSE246622_RMA']:
            d=prediction[prediction.dataset.eq(view)]
            self.assertEqual(set(d.sample_id),ids)
            self.assertTrue(d.gene_reference_excludes_test.all())
            self.assertTrue(d.model_scaling_fitted_on_training.all())
            for model in ['Baseline','Augmented']:
                self.assertEqual(d.loc[d.model.eq(model),'sample_id'].nunique(),115)

    def test_S7E_preserves_culture_label_permutations(self):
        self.assertAlmostEqual(self.report['S7E_whole_blood_rho'],.4460869565217391,places=12)
        self.assertLess(self.report['S7E_max_numeric_error'],1e-12)

if __name__=='__main__':unittest.main()
