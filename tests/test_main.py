"""
Integration tests
"""

import filecmp
import gzip
import importlib.resources as pkg_resources
import shlex
import subprocess
import tempfile
from pathlib import Path
from unittest import TestCase

from tests import resources


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

    def test_end_to_end(self):
        """
        Test that the tool produces the expected outputs
        """

        # resource names
        refgz = 'chr17_BRCA1.fa.gz'
        npu_ = 'HCC1937_n_BRCA1.pileup.gz'
        tpu_ = 'HCC1937_t_BRCA1.pileup.gz'
        out_ = 'out.tsv'

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
                    f'python3 bin/SBSCall.py -t 4 -x chr17 -r {ref} '
                    f'-n {npu} -c {tpu} -o {tempdir}/out'),
                check=True)

            outfile = f'{tempdir}/out.tsv'
            with pkg_resources.path(resources, out_) as expected:
                if not filecmp.cmp(outfile, expected):
                    msg = [f'{outfile} differs from expected:']
                    found = False
                    with open(outfile) as f1, open(expected) as f2:
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
                                f'first {n} lines match but {outfile} has '
                                f'{ll1} lines and {expected} has {ll2}')
                    self.fail('\n'.join(msg))
