import csv, math, tempfile, unittest
from pathlib import Path
from unittest.mock import patch, Mock
import numpy as np
from tap2_metrics import compute, in_cdr
from tap2_profiles import score, paper_flags
from tap2_web import read_pairs, parse_result, run

def residue(ch,n,aa,atom,xyz):
    return dict(chain=ch,number=n,insertion='',aa=aa,name={'K':'LYS','E':'GLU','D':'ASP'}[aa],
                atoms={atom:np.array(xyz,dtype=float)})

class MetricsTests(unittest.TestCase):
    def test_web_cache_cannot_silently_reuse_other_sequences(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td); (p/'pair').mkdir()
            (p/'pair/input.json').write_text('{"pair_sha256":"old"}')
            session=Mock(); session.headers={}; session.get.return_value.text='<html></html>'
            with patch('tap2_web.requests.Session',return_value=session):
                with self.assertRaisesRegex(ValueError,'cache sequence mismatch'):
                    run([dict(id='pair',name='test',H='AAA',L='CCC',pair_sha256='new')],p)
            session.post.assert_not_called()

    def test_cached_web_results_rebuild_aggregate(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td); (p/'pair').mkdir()
            (p/'pair/result.json').write_text('{"pair_sha256":"test","metrics":{}}')
            session=Mock(); session.headers={}; session.get.return_value.text='<html></html>'
            with patch('tap2_web.requests.Session',return_value=session):
                result=run([dict(id='pair',name='test',H='AAA',L='CCC',pair_sha256='test')],p)
            self.assertTrue((p/'web_results.json').is_file())
            self.assertEqual(len(result),1)
            session.post.assert_not_called()
    def test_cpu_refinement_properties_are_dictionary_without_changing_dependency(self):
        from tap2_compat import load_refine
        source="def refine(n_threads):\n    return {'Threads', str(n_threads)}\n"
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'refine.py'; p.write_text(source)
            module=Mock(__file__=str(p),__package__='ImmuneBuilder')
            with patch('tap2_compat.importlib.import_module',return_value=module):
                refine,info=load_refine()
            self.assertEqual(refine(4),{'Threads':'4'})
            self.assertTrue(info['cpu_threads_set_to_dict'])
            self.assertEqual(p.read_text(),source)
    def test_salt_bridge_changes_paper_charges_but_web_keeps_raw_charge(self):
        rr=[residue('H',27,'K','NZ',(0,0,0)),residue('L',27,'E','OE1',(3,0,0)),
            residue('L',28,'D','OD1',(3,5,0))]
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'fixture'; p.write_bytes(b'fixture')
            with patch('tap2_metrics.load_fv',return_value=rr),patch('tap2_metrics.check_imgt'), \
                 patch('tap2_metrics.sasa',return_value=(np.array([[10.,100.]]*3),'')):
                paper=score(p,p,'paper'); web=score(p,p,'web-compatible')
        self.assertEqual(paper['metrics']['L_tot'],3)
        self.assertEqual(len(paper['salt_bridges']),1)
        self.assertEqual(paper['metrics']['PNC'],0.)
        self.assertAlmostEqual(web['metrics']['PNC'],2/25)
        self.assertEqual(paper['metrics']['SFvCSP'],0.)
        self.assertEqual(web['metrics']['SFvCSP'],-2.)
        self.assertAlmostEqual(paper['residues'][0]['hydrophobicity'],1+4.1/9)
    def test_imgt_insertions_and_anchor_are_separate(self):
        self.assertTrue(in_cdr(111)); self.assertFalse(in_cdr(26)); self.assertTrue(in_cdr(26,2))
        self.assertFalse(in_cdr(120,2))
    def test_paper_threshold_boundaries(self):
        self.assertEqual(paper_flags({'L_tot':37})['L_tot'],'AMBER')
        self.assertEqual(paper_flags({'L_tot':36})['L_tot'],'RED')
        self.assertEqual(paper_flags({'L_tot':43})['L_tot'],'GREEN')
        self.assertEqual(paper_flags({'PPC':4.22})['PPC'],'AMBER')
        self.assertEqual(paper_flags({'PPC':4.2201})['PPC'],'RED')
        with self.assertRaises(ValueError): paper_flags({'PNC':float('nan')})
    def test_csv_pairing_filter_and_duplicate_handling(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'input.csv'
            p.write_text('Specific,name,hseq,lseq\nmouse,original,AAA,CCC\nhumanization,a,DDD,EEE\nhumanization,duplicate,DDD,EEE\nhumanization,b,FFF,GGG\n')
            r=read_pairs(p,'humanization',10)
            self.assertEqual([(x['H'],x['L']) for x in r],[('DDD','EEE'),('FFF','GGG')])
            self.assertEqual([x['source_row'] for x in r],[3,5])
    def test_invalid_amino_acid_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'input.csv'; p.write_text('type,name,h_seq,l_seq\nhumanized,x,AA*,CCC\n')
            with self.assertRaises(ValueError): read_pairs(p)
    def test_partial_service_response_not_accepted(self):
        values,_,links=parse_result('<table><tr><td>Total CDR Length</td><td>43</td></tr></table>')
        self.assertEqual(values,{'L_tot':43.}); self.assertFalse(links)
    def test_crashed_post_is_not_sent_twice(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td); (p/'pair').mkdir(); (p/'pair/submission_started.json').write_text('{}')
            session=Mock(); session.headers={}; session.get.return_value.text='<html></html>'
            with patch('tap2_web.requests.Session',return_value=session):
                with self.assertRaisesRegex(RuntimeError,'Uncertain prior POST'):
                    run([dict(id='pair',name='test',H='AAA',L='CCC',pair_sha256='test')],p)
            session.post.assert_not_called()

if __name__=='__main__': unittest.main()
