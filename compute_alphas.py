#!/usr/bin/env python
from glob import glob
import re,sys
from eventlet.green import subprocess
from eventlet.greenpool import GreenPool

# This is a one-off script, no real point in
# a config file
infiles = glob('SNAP_ASCII/RUN_ANALYSIS001/gen*.clones')
outfiles = glob('SNAP_ASCII/RUN_ANALYSIS002/gen*.clones')
infiles.sort()
outfiles.sort()

# ifname = file with the "ingroups" (tested for polymorphism)
# ofname = file with the "outgroups" (diff between in/out => divergance)
have_printed_header = False
MonKeyTestPath = glob('analys*/src/MonKeyTest')[0]
def process_generation(ifname,ofname):
    global have_printed_header
    gen_num = int(re.findall('gen(\d+)',ifname)[0])
    gen_num2 = int(re.findall('gen(\d+)',ofname)[0])
    assert(gen_num==gen_num2)
    args = [MonKeyTestPath,'-1',ifname,'-2',ofname]
    mkout = subprocess.check_output(args)
    mklines = mkout.splitlines()
    if not have_printed_header:
        sys.stdout.write("generation\t"+mklines[0]+'\n')
        have_printed_header = True
    for line in mklines[1:]:
        if line=='': next
        sys.stdout.write("%i\t%s\n"%(gen_num,line))

gp = GreenPool(size=10)
for ifname,ofname in zip(infiles,outfiles):
    process_generation(ifname,ofname)
    # gp.spawn(process_generation,ifname,ofname)
gp.waitall()

