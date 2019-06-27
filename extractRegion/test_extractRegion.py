#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import subprocess
import unittest

class MpTestCase(unittest.TestCase):
    """ tests for happie
    """

    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "tmp_test_output")
        self.multifasta = os.path.join(
            os.path.dirname(__file__), "multientry.fasta")

    def test_simple(self):
        cmd = "extractRegion NODE_1_length_778_cov_41.3402:1:10 -f {self.multifasta}".format(**locals())
        res =  subprocess.run([cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
        print(res.stdout.decode())
        assert res.stdout.decode().split("\n")[0] ==  ">NODE_1_length_778_cov_41.3402:1:10"
        assert res.stdout.decode().split("\n")[1] ==  "ATGGGAAAAG"

    def test_pipe(self):
        cmd = "echo 'NODE_1_length_778_cov_41.3402:1:10'  | extractRegion -f {self.multifasta}".format(**locals())
        res =  subprocess.run([cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
        print(res.stdout.decode())
        assert res.stdout.decode().split("\n")[0] ==  ">NODE_1_length_778_cov_41.3402:1:10"
        assert res.stdout.decode().split("\n")[1] ==  "ATGGGAAAAG"

    def test_pipe_w_name(self):
        cmd = "echo 'TES@NODE_1_length_778_cov_41.3402:1:10'  | extractRegion -f {self.multifasta}".format(**locals())
        res =  subprocess.run([cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
        print(res.stdout.decode())
        assert res.stdout.decode().split("\n")[0] ==  ">TES@NODE_1_length_778_cov_41.3402:1:10"
        assert res.stdout.decode().split("\n")[1] ==  "ATGGGAAAAG"

    def test_rc(self):
        cmd = "echo 'TES-RC@NODE_1_length_778_cov_41.3402:1:10'  | extractRegion -f {self.multifasta}".format(**locals())
        res =  subprocess.run([cmd],
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
        print(res.stdout.decode())
        assert res.stdout.decode().split("\n")[0] ==  ">TES-RC@NODE_1_length_778_cov_41.3402:1:10"
        assert res.stdout.decode().split("\n")[1] ==  "CTTTTCCCAT"
