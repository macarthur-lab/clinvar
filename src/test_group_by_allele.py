import unittest
from group_by_allele import group_alleles, group_by_allele
from StringIO import StringIO
from parse_clinvar_xml import HEADER


class TestGroupByAllele(unittest.TestCase):

    def setUp(self):
        self.header = HEADER

        self.r1 = ['chrom1', 'pos1', 'ref1', 'alt1']
        for h in self.header[4:]:
            self.r1.append(h + "1a;" + h + "1b")

        # should be joined with r1, but data is different
        self.r2 = ['chrom1', 'pos1', 'ref1', 'alt1']
        for h in self.header[4:]:
            self.r2.append(h + "1c;" + h + "1d")

        # should not be joined with r1/r2 because alt is different
        self.r3 = ['chrom1', 'pos1', 'ref1', 'alt2']
        for h in self.header[4:]:
            self.r3.append(h + "2a;" + h + "2b")

        # TODO pathogenic, benign, conflicted, gold_stars should be re-evaluated?

    def test_group_alleles(self):
        data1 = dict(zip(self.header, self.r1))
        data2 = dict(zip(self.header, self.r1))
        data3 = dict(zip(self.header, self.r2))
        data4 = dict(zip(self.header, self.r3))

        combined_data = group_alleles(data1, data2)
        combined_data2 = group_alleles(combined_data, data3)

        with self.assertRaises(ValueError):
            group_alleles(data1, data4)  # different alts!

        for h in self.header[4:]:
            if h.endswith('_ordered'):
                # no deduping should have occurred
                self.assertEqual(
                    combined_data2[h],
                    "{h}1a;{h}1b;{h}1a;{h}1b;{h}1c;{h}1d".format(h=h))
            else:
                # deduping *should* have occurred
                self.assertEqual(
                    combined_data2[h],
                    "{h}1a;{h}1b;{h}1c;{h}1d".format(h=h))

    def test_group_by_allele(self):
        row_data = [self.header, self.r1, self.r1, self.r2, self.r3]
        input_rows = ["\t".join(row) + "\n" for row in row_data]

        outfile = StringIO()
        group_by_allele(iter(input_rows), outfile)

        for i, output_row in enumerate(outfile.getvalue().split('\n')):
            if i == 0:
                output_header = output_row.split('\t')
                self.assertEqual(output_header, self.header)

            elif i == 1:
                combined_data = dict(zip(output_header, output_row.split('\t')))
                for h in self.header[4:]:
                    if h.endswith('_ordered'):
                        # no deduping should have occurred
                        self.assertEqual(
                            combined_data[h],
                            "{h}1a;{h}1b;{h}1a;{h}1b;{h}1c;{h}1d".format(h=h))
                    else:
                        # deduping *should* have occurred
                        self.assertEqual(
                            combined_data[h],
                            "{h}1a;{h}1b;{h}1c;{h}1d".format(h=h))

            elif i == 2:
                # this row should not have been merged with previous ones
                self.assertEqual(self.r3, output_row.split('\t'))

            elif i == 3:
                self.assertEqual(output_row, '')


if __name__ == '__main__':
    unittest.main()
