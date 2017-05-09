import unittest
import logging
from mock import patch, MagicMock, mock
from tempfile import NamedTemporaryFile
import os
import filecmp
from cohort_matcher import checkConfig, main, parseArguments, readSamples, vcfToIntervals

class TestCohortMatcher(unittest.TestCase):
    @patch('os.path.isdir')
    @patch('os.path.exists')
    def test_checkConfig(self, mock_exists, mock_isdir):
        # Set up test case
        config = MagicMock(log_level=logging.INFO, reference='', reference2=None, vcf2=None, chromosome_map=None)
        # Set up supporting mocks
        # Test
        retVal = checkConfig(config)
        # Check results
        self.assertTrue(retVal)

    def test_compareSamples(self):
        self.skipTest("not yet implemented")

    def test_genotypeSamples(self):
        self.skipTest("not yet implemented")

    @patch('cohort_matcher.parseArguments')
    @patch('cohort_matcher.checkConfig')
    @patch('cohort_matcher.readSamples')
    @patch('cohort_matcher.vcfToIntervals')
    @patch('cohort_matcher.genotypeSamples')
    @patch('cohort_matcher.compareSamples')
    def test_main(self, mock_compareSamples, mock_genotypeSamples,
                  mock_vcfToIntervals, mock_readSamples,
                  mock_checkConfig, mock_parseArguments):
        # Set up test case
        argv = []
        # Set up supporting mocks
        # Test
        with patch('logging.basicConfig') as mock_basicConfig:
            retval = main(argv)
        # Check results
        self.assertEqual(retval, 0)

    @patch('os.path.isfile')
    def test_readSamples(self, mock_isfile):
        # Set up test parameters
        sampleSheetFile = MagicMock(name="sampleSheetFile")
        # Set up supporting mocks
        mock_isfile.return_value = True
        sampleSheetText = ['sample1\ts3://some/path/to/sample1.bam\n',
                           'sample2\ts3://some/path/to/sample2.bam\n']
        mock_open = mock.mock_open(read_data=''.join(sampleSheetText))
        mock_open.return_value.__iter__ = lambda self: iter(self.readline, '')
        # Test
        with patch('__builtin__.open', mock_open):
            actual_sampleSheet = readSamples(sampleSheetFile)
        # Check results
        self.assertEqual(len(actual_sampleSheet), 2)
        self.assertEqual(actual_sampleSheet[0]['name'], 'sample1')
        self.assertEqual(actual_sampleSheet[1]['name'], 'sample2')
        self.assertEqual(actual_sampleSheet[0]['bam'], 's3://some/path/to/sample1.bam')
        self.assertEqual(actual_sampleSheet[1]['bam'], 's3://some/path/to/sample2.bam')
        
    def test_parseArguments(self):
        # Set up test case
        args = ['-S1', 'set1', '-S2', 'set2', '-V', 'somevcf', '-R', 'someref']
        # Set up supporting mocks
        # Test
        actual_args = parseArguments(args)
        # Check results
        self.assertEqual(actual_args.set1, 'set1')
        self.assertEqual(actual_args.set2, 'set2')
        self.assertEqual(actual_args.vcf, 'somevcf')
        self.assertEqual(actual_args.reference, 'someref')

    def test_vcfToIntervals(self):
        # Set up test parameters
        vcfFile = "hg19.exome.highAF.1511.vcf"
        bedFile = NamedTemporaryFile(delete=False)
        # Set up supporting mocks
        # Test
        vcfToIntervals(vcfFile, bedFile.name)
        # Check results
        self.assertTrue(filecmp.cmp(bedFile.name, "test_data/hg19.exome.highAF.1511.bed"))
        # Clean up
        os.remove(bedFile.name)