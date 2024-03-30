#!/usr/bin/python3

# Usage: python3 frequence.py fichier_texte
import sys

alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
Occurences = {}
length = 0


f = open(sys.argv[1])
texte = f.read()

for c in texte:
	if c in Occurences:
		Occurences[c] += 1.0
	else:
		Occurences[c] = 1.0
	length += 1

# Print the frequences
for c in alphabet:
	if c in Occurences:
		print (c, Occurences[c] / length)
	else:
		print (c, 0.0)



