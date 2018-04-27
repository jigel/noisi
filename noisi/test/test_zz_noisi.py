import os
import sys


def test_zz_noisi():


	if os.path.exists("test/testdata/testsrc/step_0/corr"):
		os.system("rm -rf test/testdata/testsrc/step_0/corr")

	if os.path.exists("test/testdata/testsrc/step_0/adjt"):
		os.system("rm -rf test/testdata/testsrc/step_0/adjt")

	if os.path.exists("test/testdata/testsrc/step_0/kern"):
		os.system("rm -rf test/testdata/testsrc/step_0/kern")

	if os.path.exists("test/testdata/testsrc/step_0/grad"):
		os.system("rm -rf test/testdata/testsrc/step_0/grad")

