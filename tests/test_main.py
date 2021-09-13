from unittest import TestCase
import shlex
import subprocess


class TestMain(TestCase):                                                       
    """                                                                         
    Integration tests                                                           
    """                                                                         
    def test_help(self):                                                      
        """
        Test that the script can be called by the documented name
        """
        result = subprocess.run(shlex.split('python3 bin/SBSCall.py --help'),
                                capture_output=True, text=True)
        self.assertTrue(result.stdout.startswith('usage: SBSCall.py'))
