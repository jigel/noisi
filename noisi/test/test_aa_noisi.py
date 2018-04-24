import os
import sys

os.system("noisi")
if not os.path.exists("test/testdata/testsrc"):
	sys.exit("\n\n\nGo to noisi/noisi directory to run tests.\n\n\n")

def test_aa_noisi():


	if not os.path.exists("test/testdata/testsrc/step_0/corr"):
		os.mkdir("test/testdata/testsrc/step_0/corr")

	if not os.path.exists("test/testdata/testsrc/step_0/adjt"):
		os.mkdir("test/testdata/testsrc/step_0/adjt")

	if not os.path.exists("test/testdata/testsrc/step_0/kern"):
		os.mkdir("test/testdata/testsrc/step_0/kern")

	if not os.path.exists("test/testdata/testsrc/step_0/grad"):
		os.mkdir("test/testdata/testsrc/step_0/grad")

