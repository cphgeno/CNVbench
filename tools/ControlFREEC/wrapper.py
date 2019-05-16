#!/usr/bin/env python3

import argparse
import subprocess

parser = argparse.ArgumentParser(description='Wrapper for running ControlFREEC. Check tool documentation for further help.')
parser.add_argument('--input', action='store')
parser.add_argument('--sex', action='store')
parser.add_argument('--outputDir', action='store')
parser.add_argument('--maxThreads', action='store')
parser.add_argument('--sambamba', action='store')
parser.add_argument('--SambambaThreads', action='store')
parser.add_argument('--coefficientOfVariation', action='store')
parser.add_argument('--minCNAlength', action='store')
parser.add_argument('--forceGCcontentNormalization', action='store')
parser.add_argument('--ploidy', action='store')
parser.add_argument('--chrFiles', action='store')
parser.add_argument('--chrLenFile', action='store')
args = parser.parse_args()


configfile = '{path}/config.txt'.format(path = args.outputDir)
with open(configfile, "w") as c:
	print("[general]\n", file = c)
	for arg, value in vars(args).items():
		if arg == "input":
			continue
		print(' = '.join([arg, value]), file = c)
	print("\n[sample]\n", file = c)
	print("mateFile = {}".format(args.input), file = c)
	print("inputFormat = {}".format(args.input.rsplit('.', 1)[1].upper()), file = c)
	print("mateOrientation = FR", file = c)

subprocess.check_call(["freec", "-conf", configfile])

