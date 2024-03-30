#!/usr/bin/python3


# Usage: python3 cesar.py clef c/d phrase
# Returns the result without additional text

import sys
cle = sys.argv[1]
texte = sys.argv[3]
res = ''
decalage = ord(str(cle))

if sys.argv[2] == 'c':
  
    for c in texte:
        res += chr((decalage + ord(c))%(ord('Z') - ord('A')+1) + ord('A'))
    print(res)

if sys.argv[2] == 'd':
    

    for c in texte:
        res += chr((-decalage + ord(c))%(ord('Z') - ord('A')+1) + ord('A'))
    print(res)