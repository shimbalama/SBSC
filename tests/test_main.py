"""
Integration tests
"""

import filecmp
import gzip
import importlib.resources as pkg_resources
import json
import shlex
import subprocess
import tempfile
from pathlib import Path
from unittest import TestCase

from tests import resources

from jsondiff import diff


class TestMain(TestCase):
    """
    Integration tests
    """

    def test_help(self):
        """
        Test that the tool can be called by the documented name
        """
        result = subprocess.run(shlex.split('python3 bin/SBSCall.py --help'),
                                capture_output=True, text=True, check=True)
        self.assertTrue(result.stdout.startswith('usage: SBSCall.py'))

    def test_call(self):
        """
        Test that the call subcommand produces the expected outputs
        """

        # resource names
        refgz = 'chr17_RNF157.fa.gz'
        npu_ = 'HCC1937_n_RNF157.pileup.gz'
        tpu_ = 'HCC1937_t_RNF157.pileup.gz'
        filt_ = 'output_001_9_12_10.tsv'
        unfilt_ = 'output.json'

        with tempfile.TemporaryDirectory() as tempdir:

            # set up inputs as real files on the filesystem
            ref = f'{tempdir}/{Path(refgz).stem}'
            with open(ref, 'wb') as fd:
                byts = pkg_resources.read_binary(resources, refgz)
                fd.write(gzip.decompress(byts))
            npu = f'{tempdir}/{npu_}'
            tpu = f'{tempdir}/{tpu_}'
            for ext in ('', '.tbi'):
                for rsrc, fil in ((npu_, npu), (tpu_, tpu)):
                    with open(fil + ext, 'wb') as fd:
                        byts = pkg_resources.read_binary(resources, rsrc + ext)
                        fd.write(byts)

            subprocess.run(
                shlex.split(
                    f'python3 bin/SBSCall.py -o {tempdir}/output -x chr17 '
                    f'call -t 4 -r {ref} -n {npu} -c {tpu}'),
                check=True)

            with open(f'{tempdir}/output.json') as fd:
                unfiltered = json.load(fd)
            expected = json.loads(pkg_resources.read_text(resources, unfilt_))
            self.assertFalse(diff(unfiltered, expected))

            filtered = f'{tempdir}/output_001_9_12_10.tsv'
            with pkg_resources.path(resources, filt_) as expected:
                if not filecmp.cmp(filtered, expected):
                    msg = [f'{filtered} differs from expected:']
                    found = False
                    with open(filtered) as f1, open(expected) as f2:
                        for n, (line1, line2) in enumerate(zip(f1, f2), 1):
                            if line1 != line2:
                                found = True
                                msg.append(f'first difference at line {n}:')
                                msg.append(f' - {line2}')
                                msg.append(f' + {line1}')
                        if not found:
                            f1.seek(0)
                            f2.seek(0)
                            ll1, ll2 = len(f1.readlines()), len(f2.readlines())
                            assert ll1 != ll2
                            msg.append(
                                f'first {n} lines match but {filtered} has '
                                f'{ll1} lines and {expected} has {ll2}')
                    self.fail('\n'.join(msg))

    def test_filt(self):
        """
        Test that the filt subcommand produces the expected outputs;
        specifically that running it on the unfiltered outputs of the call
        subcommand generates the same filtered outputs
        """

        # resource names
        filt_ = 'output_001_9_12_10.tsv'
        unfilt_ = 'output.json'

        with tempfile.TemporaryDirectory() as tempdir:

            # set up input as real file on the filesystem
            unfilt = f'{tempdir}/{unfilt_}'
            with open(unfilt, 'wb') as fd:
                byts = pkg_resources.read_binary(resources, unfilt_)
                fd.write(byts)

            subprocess.run(
                shlex.split(
                    f'python3 bin/SBSCall.py -o {tempdir}/output '
                    f'filt {unfilt}'),
                check=True)

            filtered = f'{tempdir}/output_001_9_12_10.tsv'
            with pkg_resources.path(resources, filt_) as expected:
                if not filecmp.cmp(filtered, expected):
                    msg = [f'{filtered} differs from expected:']
                    found = False
                    with open(filtered) as f1, open(expected) as f2:
                        for n, (line1, line2) in enumerate(zip(f1, f2), 1):
                            if line1 != line2:
                                found = True
                                msg.append(f'first difference at line {n}:')
                                msg.append(f' - {line2}')
                                msg.append(f' + {line1}')
                        if not found:
                            f1.seek(0)
                            f2.seek(0)
                            ll1, ll2 = len(f1.readlines()), len(f2.readlines())
                            assert ll1 != ll2
                            msg.append(
                                f'first {n} lines match but {filtered} has '
                                f'{ll1} lines and {expected} has {ll2}')
                    self.fail('\n'.join(msg))
