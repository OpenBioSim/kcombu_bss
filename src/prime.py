#!/usr/bin/env python

import sys
import os
import math

if (len(sys.argv)<2):
  print "prime.py [N]"
  sys.exit(1)

N = int(sys.argv[1])
nprime = 0
for i in range (2,N):
  j = 2
  div = 0 
  while ((div==0) and (j<i)): 
    if ((i%j)==0):
      div = 1
    else:
      j += 1
  if (div==0):
    sys.stdout.write("%4d,"%(i))
    nprime += 1
    if ((nprime%20)==0):
      sys.stdout.write("\n")
sys.stdout.write("\n")
