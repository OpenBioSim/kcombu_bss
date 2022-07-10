#!/usr/bin/env python

import sys
import os
import math

LastModDate = "Jan 8, 2016"

if (len(sys.argv)<2):
  print "decompose_prime.py [N]"
  print " coded by T.Kawabata. LastModDate:%s"%(LastModDate)
  sys.exit(1)

N = int(sys.argv[1])
nprime = 0

prime = {}

for i in range (2,N+1):
  prime[i] = 1 

for i in range (2,N+1):
  for j in range(2,(N+1)/i + 1):
    k = i *j
    prime[k] = 0

#for i in range (2,N+1):
#  if (prime[i]==1):
#    print "PRIME %d"%(i)

for i in range (2,N+1):
  if (prime[i]==1) and ((N%i)==0):
    shou = N/i
    Npow = 1
    while (shou%i==0):
      shou = shou/i
      Npow += 1
    sys.stdout.write("%4d ^ %2d\n"%(i,Npow))
