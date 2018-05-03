import os
import sys

os.system("noisi")
if not os.path.exists("test/testdata/testsrc"):
    sys.exit("\n\n\nGo to noisi/noisi directory to run tests.\n\n\n")

def test_aa_noisi():

    os.system('cp test/testdata/config_archived.json test/testdata/config.json')
    os.system('cp test/testdata/testsrc/source_config_archived.json \
test/testdata/testsrc/source_config.json')
    os.system('cp test/testdata/testsrc/measr_config_archived.json \
test/testdata/testsrc/measr_config.json')
    
    os.system('rm -rf test/testdata/testsrc/wavefield_processed/')
    
    os.system("rm -rf test/testdata/testsrc/step_0/corr")
    os.mkdir("test/testdata/testsrc/step_0/corr")

    os.system("rm -rf test/testdata/testsrc/step_0/adjt")
    os.mkdir("test/testdata/testsrc/step_0/adjt")

    os.system("rm -rf test/testdata/testsrc/step_0/kern")
    os.mkdir("test/testdata/testsrc/step_0/kern")

    os.system("rm -rf test/testdata/testsrc/step_0/grad")
    os.mkdir("test/testdata/testsrc/step_0/grad")


