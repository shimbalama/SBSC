"""
Unit tests for SBSC.output.filt module
"""
import importlib.resources as pkg_resources
import json
import tempfile

from copy import deepcopy
from unittest import TestCase

import numpy as np
import pandas as pd

from tests import resources
from tests.test_utils import filediff
from SBSC.output.filt import distance, mergepos, filt

# this data has 953 SNV records, of which 57 pass default filters.
# Of the 57 passing, there are 4 expected MNV alleles: 
# 121191 TG>AA,AC
# 121454 TC>AA
# 121692 GA>TT
RAW = json.loads(pkg_resources.read_text(resources,
                                         'chr22_KI270734v1_random.json'))

class Expando(object):
    pass

# mock args with defalt values
ARGS = Expando()
ARGS.MNVs = True
ARGS.P_value = 0.001
ARGS.min_depth_tumour = 12
ARGS.min_depth_normal = 10
ARGS.min_mean_base_quality = 9

class TestMNVs(TestCase):
    """
    Unit tests for MNV functions
    """
    def test_distance_identity(self):
        r1 = {'chrom': 'chr1', 'pos': 1}
        r2 = {'chrom': 'chr1', 'pos': 1}
        self.assertEqual(distance(r1, r2), 0)

    def test_distance_contiguous(self):
        r1 = {'chrom': 'chr1', 'pos': 1}
        r2 = {'chrom': 'chr1', 'pos': 2}
        self.assertEqual(distance(r1, r2), 1)

    def test_distance_empty(self):
        r1 = {}
        r2 = {'chrom': 'chr1', 'pos': 1}
        self.assertEqual(distance(r1, r2), np.inf)

    def test_distance_diff_chrom(self):
        r1 = {'chrom': 'chr1', 'pos': 1}
        r2 = {'chrom': 'chr2', 'pos': 1}
        self.assertEqual(distance(r1, r2), np.inf)

    def test_mergepos_monoallelic_snv(self):
        self.assertEqual(mergepos([['A']]), ['A'])

    def test_mergepos_polyallelic_snv(self):
        self.assertEqual(mergepos([['A', 'C']]), ['A', 'C'])

    def test_mergepos_monoallelic_mnv(self):
        self.assertEqual(mergepos([['A'], ['C']]), ['AC'])

    def test_mergepos_polyallelic_mnv(self):
        self.assertEqual(mergepos([['A', 'C'], ['A'], ['C', 'G']]),
                         ['AAC', 'AAG', 'CAC', 'CAG'])

    def test_filt_mnvs(self):
        with tempfile.TemporaryDirectory() as tempdir:
            args = deepcopy(ARGS)
            args.output = f'{tempdir}/output'
            filt(RAW, args)
            filtered = f'{tempdir}/output_001_9_12_10.tsv'
            with pkg_resources.path(
                    resources,
                    'chr22_KI270734v1_random_001_9_12_10.tsv') as expected:
                diff = filediff(filtered, expected)
                if diff:
                    self.fail(diff)

    def test_filt_no_mnvs(self):
        with tempfile.TemporaryDirectory() as tempdir:
            args = deepcopy(ARGS)
            args.MNVs = False
            args.output = f'{tempdir}/output'
            filt(RAW, args)
            filtered = f'{tempdir}/output_001_9_12_10.tsv'
            with pkg_resources.path(
                    resources,
                    'chr22_KI270734v1_random_no-MNVs_001_9_12_10.tsv') as expected:
                diff = filediff(filtered, expected)
                if diff:
                    self.fail(diff)
