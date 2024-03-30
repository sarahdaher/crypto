#!/usr/bin/python3

# Usage: python3 subst_mono.py clef c/d phrase
# Returns the result without additional text


import sys
cle = sys.argv[1]
texte = sys.argv[3]
res = ''


if sys.argv[2] == 'c':
  
    for c in texte:
        res += cle[ord(c) - ord('A')]
    print(res)

if sys.argv[2] == 'd':
    
    for c in texte:
        res += chr(cle.index(c) + ord('A'))
    print(res)