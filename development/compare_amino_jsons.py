#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 11:46:37 2021

@author: colin
"""

import json
import argparse
import os

def is_file(string):
	if os.path.isfile(string):
		return string
	else:
		raise ValueError("Not a file")

parser = argparse.ArgumentParser(description="compare two amino alphabet clusterings")

parser.add_argument('-f', '--files', type=is_file, nargs=2, default=None,
                    help='pass two .json file names seperated by a space')

args = parser.parse_args()

f1 = args.files[0]
f2 = args.files[1]

file1 = open(f1)
file2 = open(f2)

dict1 = json.load(file1)
dict2 = json.load(file2)

print('Amino Alphabet Differences for F1='+f1+', F2='+f2)
print('Amino\tF1\tF2')

for key, value in dict1.items():
    if dict2[key] is not value:
        print(key+'\t'+value+'\t'+dict2[key])

file1.close()
file2.close()