#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
import argparse
import matplotlib.pyplot as plt
import sys
import os

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", required = True, metavar = "input.sam", help = "Alignment on SAM format as input")
	parser.add_argument("-o", "--output", required = True, metavar = "output.png", help = "Read length distribution picture as output")
	args = parser.parse_args()

	if os.path.exists(args.input):
		inPut = args.input
	else:
		print "No such file : " + args.input
		sys.exit(1)

	graphOutputName = args.output

	reads_lgts = []
	entree = sys.stdin

	for line in entree:
		reads_lgts.append(int(line.strip()))

	L = max(reads_lgts) - min(reads_lgts)

	if L == 0:
		L = 1

	plt.hist(reads_lgts, histtype='bar', facecolor='b', bins = L)

	title = 'Reads length distribution from ' + inPut
	Xlabel = 'Length' + '\n' + 'Bins = ' + str(L)

	plt.title(title)
	plt.xlabel(Xlabel)
	plt.ylabel('Number of reads')
	plt.grid(True)
	plt.xlim(min(reads_lgts), max(reads_lgts) + 5)
	plt.savefig(graphOutputName)
