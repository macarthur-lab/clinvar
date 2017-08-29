import unittest
from group_by_allele import group_alleles, group_by_allele
from pprint import pprint
from StringIO import StringIO
from parse_clinvar_xml import HEADER

class TestGroupByAllele(unittest.TestCase):

    def setUp(self):
        self.header = HEADER
        self.r1=['1','55518316','C','T','55518316','55518316','+','Variant','36672','RCV000030351','SCV000323050;SCV000053018;SCV000358232','45333','PCSK9','NM_174936.3:c.658-7C>T','','NM_174936.3:c.658-7C>T:intron variant','Conflicting interpretations of pathogenicity','benign;benign;likely benign','0','0','0','1','2','criteria provided, conflicting interpretations','criteria provided, single submitter;criteria provided, single submitter;criteria provided, single submitter','14/06/2016','Laboratory Corporation of America,;Cardiovascular Research Group,Instituto Nacional de Saude Doutor Ricardo Jorge;Illumina Clinical Services Laboratory,Illumina','Laboratory Corporation of America,;Cardiovascular Research Group,Instituto Nacional de Saude Doutor Ricardo Jorge;Illumina Clinical Services Laboratory,Illumina','Familial hypercholesterolemia;Familial Hypercholesterolemia;Familial Hypercholesterolemia','12730697;15177124;17094996;21600525;22364837;23725921;23788249;24404629;24418289;24636176;25053660;25356965;25404096;25741868;27854360','Autosomal dominant inheritance','','1-9 / 1 000 000','','germline;not provided','GeneReviews:NBK174884;Genetic Alliance:Familial+Hypercholesterolemia/2746;Genetic Testing Registry (GTR):GTR000203962;Genetic Testing Registry (GTR):GTR000260633;Genetic Testing Registry (GTR):GTR000320944;Genetic Testing Registry (GTR):GTR000321567;Genetic Testing Registry (GTR):GTR000500234;Genetic Testing Registry (GTR):GTR000500235;Genetic Testing Registry (GTR):GTR000500236;Genetic Testing Registry (GTR):GTR000512385;Genetic Testing Registry (GTR):GTR000520922;Genetic Testing Registry (GTR):GTR000521881;Genetic Testing Registry (GTR):GTR000523352;Genetic Testing Registry (GTR):GTR000523354;Genetic Testing Registry (GTR):GTR000525924;Genetic Testing Registry (GTR):GTR000528325;Genetic Testing Registry (GTR):GTR000528380;Genetic Testing Registry (GTR):GTR000528697;Genetic Testing Registry (GTR):GTR000528701;Genetic Testing Registry (GTR):GTR000552211;MedGen:C0020445;OMIM:143890;OMIM:144400;OMIM:600946.0028;Orphanet:391665;SNOMED CT:397915002;SNOMED CT:398036000','2011-08-18;2016-03-01;2016-06-14']
        self.r2=['1','55518316','C','T','55518316','55518316','+','Variant','36672','RCV000252912','SCV000519557;SCV000316578','45333','PCSK9','NM_174936.3:c.658-7C>T','','NM_174936.3:c.658-7C>T:intron variant','Benign','benign;benign','0','0','0','0','2','criteria provided, multiple submitters, no conflicts','criteria provided, single submitter;criteria provided, single submitter','28/09/2016','PreventionGenetics,PreventionGenetics;GeneDx','PreventionGenetics,PreventionGenetics;GeneDx','not specified;NOT SPECIFIED;not specified','12730697;17094996;23788249;25741868;27854360','','','','','germline','MedGen:CN169374','0000-00-00;2016-09-28']
        self.r3=['1','55518316','C','T','55518316','55518316','+','Variant','36672','RCV000300438','SCV000358231','45333','PCSK9','NM_174936.3:c.658-7C>T','','NM_174936.3:c.658-7C>T:intron variant','Likely benign','likely benign','0','0','0','1','0','criteria provided, single submitter','criteria provided, single submitter','14/06/2016','Illumina Clinical Services Laboratory,Illumina','Illumina Clinical Services Laboratory,Illumina','Familial hypobetalipoproteinemia;Familial Hypobetalipoproteinemia','12730697;17094996;23788249;27854360','','','','','germline','Genetics Home Reference:familial-hypobetalipoproteinemia;MedGen:C1862596;SNOMED CT:60193003','14/06/2016']

    def test_group_alleles(self):
        data1 = dict(zip(self.header, self.r1))
        data2 = dict(zip(self.header, self.r2))
        data3 = dict(zip(self.header, self.r3))

        combined_data = group_alleles(data1, data2)
        combined_data2 = group_alleles(combined_data, data3)

        self.assertEqual(combined_data2['rcv'], "RCV000030351;RCV000252912;RCV000300438")
        self.assertEqual(combined_data2['scv'],"SCV000323050;SCV000053018;SCV000358232;SCV000519557;SCV000316578;SCV000358231")
        self.assertEqual(combined_data2['benign'], "4")
        self.assertEqual(combined_data2['likely_benign'],"2")

    def test_group_by_allele(self):
        input_rows = ["\t".join(row)+"\n" for row in [self.header, self.r1, self.r2, self.r3]]

        outfile = StringIO()
        group_by_allele(iter(input_rows), outfile)

        for i, output_row in enumerate(outfile.getvalue().split('\n')):
            if i == 0:
                output_header = output_row.split('\t')
                self.assertEqual(output_header, self.header)

            elif i == 1:
                combined_data = dict(zip(output_header, output_row.split('\t')))
                self.assertEqual(combined_data['rcv'], "RCV000030351;RCV000252912;RCV000300438")
                self.assertEqual(combined_data['scv'],"SCV000323050;SCV000053018;SCV000358232;SCV000519557;SCV000316578;SCV000358231")
                self.assertEqual(combined_data['benign'], "4")
                self.assertEqual(combined_data['likely_benign'], "2")

            elif i == 2:
                self.assertEqual(output_row, '')


if __name__ == '__main__':
    unittest.main()