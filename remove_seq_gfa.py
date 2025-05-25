#!/usr/bin/env python3
import sys
for line in open(sys.argv[1]):
    if line.startswith("S"):
        arr = line.strip().split("\t")
        lseq = len (arr[2])
        new_arr = [arr[0],arr[1],"*",f"LN:i:{lseq}"]
        print("\t".join(new_arr))
    else: 
        print(line.strip())

