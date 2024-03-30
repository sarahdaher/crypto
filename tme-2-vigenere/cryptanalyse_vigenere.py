# Sorbonne Université 3I024 2023-2024
# TME 2 : Cryptanalyse du chiffre de Vigenere
#
# Etudiant.e 1 : CHEMALI MAISSA 28722554
# Etudiant.e 2 : DAHER SARAH 21100791

import sys, getopt, string, math

# Alphabet français
alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

# Fréquence moyenne des lettres en français
# À modifier
freq_FR = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

Occurences = {}
length = 0

f = open("germinal_nettoye")
texte = f.read()

for c in texte:
    if c in Occurences:
        Occurences[c] += 1.0
    else:
        Occurences[c] = 1.0
    length += 1

#stocker les frequences dans freq_FR
for c in alphabet:
    if c in Occurences:
        freq_FR[ord(c) - ord('A')] = Occurences[c]/length
    else:
        freq_FR[ord(c) - ord('A')] = 0.0


# Chiffrement César
def chiffre_cesar(txt, key):
    """
    Prend un texte (txt) nettoyé et un nombre (key), et renvoie le texte dont chaque lettre est décalée à droite de key caractères
    ce qui correspond a son chiffrement de césar de clef key
    """
    res = ''
    for c in txt:
        res += chr((ord(c) + key + ord('A'))%26 + ord('A'))
    return res

# Déchiffrement César
def dechiffre_cesar(txt, key):
    """
    Prend un texte (txt) nettoyé et un nombre (key), et renvoie le texte dont chaque lettre est décalée à gauche de key caractères
    ce qui correspond a son déchiffrement de césar de clef key
    """
    res = ''
    for c in txt:
        res += chr((ord(c) - key - ord('A'))%26 + ord('A'))
    return res

# Chiffrement Vigenere
def chiffre_vigenere(txt, key):
    """
    Prend un texte (txt) nettoyé et une liste de nombres (key) et renvoie le texte correspondant a son chiffrement de vigenere de clef key
    """
    res = ''
    n = len(key)
    i = 0
    for c in txt:
        res += chr((ord(c) + key[i] + ord('A'))%26 + ord('A'))
        i= (i+1)%n
    return res

# Déchiffrement Vigenere
def dechiffre_vigenere(txt, key):
    """
    Prend un texte (txt) nettoyé et une liste de nombres (key) et renvoie le texte correspondant a son dechiffrement de vigenere de clef key
    """
    res = ''
    n = len(key)
    i = 0
    for c in txt:
        res += chr((ord(c) - key[i] - ord('A'))%26 + ord('A'))
        i= (i+1)%n
    return res

# Analyse de fréquences
def freq(txt):
    """
    Prend un texte nettoye (txt) et renvoie une liste donnant le nombre d'occurences de chaque lettres de l'alphabet (en position i de la liste, on aura le nombre d'occurences de la (i+1)eme lettre de l'alphabet)
    """
    hist=[0.0]*len(alphabet)
    for c in txt :
        hist[(ord(c) - ord('A'))%26] +=1
    return hist

# Renvoie l'indice dans l'alphabet
# de la lettre la plus fréquente d'un texte
def lettre_freq_max(txt):
    """
    Prend un texte nettoye et renvoie l'indice dans l'alphabet de sa lettre la plus frequente
    """
    occ = freq(txt)
    imax = 0
    for i in range(len(occ)) :
        if occ[i] >occ[imax] :
            imax = i
    return imax

# indice de coïncidence
def indice_coincidence(hist):
    """
    Prend un tableau de nombres (hist) contenant les occurences des caractères d'un texte et renvoie son indice de coincidence
    """
    res = 0.0
    cpt = 0
    for i in range(26):
        res += hist[i]*(hist[i] - 1)
        cpt += hist[i]
    if cpt == 0 or cpt == 1:
        return 0.0
    return res /(cpt*(cpt - 1))

# Recherche la longueur de la clé
def longueur_clef(cipher):
    """
    Prend un texte chiffre (cipher) par Vigenere, et renvoie la taille probable de sa clef de chiffrement/dechiffrement
    """

    for k in range(1,21) :   # parcours des tailles de clé possibles
        ick = 0              # stockage de la somme des indices de coincidence de chaque colonne pour une longueur de clé k
        for i in range(k) :  # traitement de chacune des k colonnes
            coli = ''
            for j in range(len(cipher)) :  
                if j%k == i:              # on recupere chaque caractere du cipher dont l'indice est de la forme i+j*k
                    coli += cipher[j]
           
            ick += indice_coincidence(freq(coli)) # on rajoute l'indice de coincidence de la colonne traitée dans la somme ick
        if (ick/k) >0.06 : # si la moyenne est > 0.06, c’est qu’on a trouvé la bonne taille de clef
            return k
    return 0
   
# Renvoie le tableau des décalages probables étant
# donné la longueur de la clé
# en utilisant la lettre la plus fréquente
# de chaque colonne
def clef_par_decalages(cipher, key_length):
    """
   Prend un texte chiffre (cipher) par Vigenere et la longueur de sa clé (key_length), et renvoie la liste probable des decalages de chaque colonne
    """
    decalages=[0]*key_length
    for i in range(key_length) :  # traitement de chacune des key_length colonnes
            coli = ''
            for j in range(len(cipher)) :  
                if j%key_length == i:              # on recupere chaque caractere du cipher dont l'indice est de la forme i+j*key_length
                    coli += cipher[j]
            ioccMax = lettre_freq_max(coli)
            decalages[i] = (ioccMax - 4)%26 # E=4
    return decalages

# Cryptanalyse V1 avec décalages par frequence max
def cryptanalyse_v1(cipher):
    """
    Prend un texte chiffre par vigenere (cipher) et renvoie un texte correspond a une tentative de dechiffrement avec décalages par frequence max
    """
    key_length = longueur_clef(cipher)
    key = clef_par_decalages(cipher, key_length)
   
    return dechiffre_vigenere(cipher, key)

# Question 9

"""
18 textes (sur 100) sont correctement cryptanalysés avec cette fonction.
Ceci peut s'expliquer par le fait qu'on a choisi la première lettre d'occurence maximale comme étant le chiffré de E (dans lettre_freq_max). Il se peut qu'on ait plusieurs lettres ayant ce même nombre d'occurences maximales, et alors cela peut induire une erreur dans le choix du décalage. Egalement, il se peut que le texte ne soit pas completement representatif de la langue (et donc E ne serait
en fait pas l'element d'occurences maximales du clair correspondant exemple : la disparition de Perec).
Dans la fonction clef_par_decalages, on appelle lettre_freq_max sur chacune des colonnes (donc key_length fois) ce qui augmente le risque d'erreur.
Cela peut également s'expliquer par le fait que lorsqu'on devine la taille de la clef, on prend la premiere taille dont l'IC est > 0.06 (longeur_clef). Mais il se peut que ce ne soit qu'une coincidence et que la réelle taille soit plus grande.
"""
################################################################


### Les fonctions suivantes sont utiles uniquement
### pour la cryptanalyse V2.

# Indice de coincidence mutuelle avec décalage
def indice_coincidence_mutuelle(h1,h2,d):
    """
    Prend deux tableaux de fréquence (h1 d'un texte 1 et h2 d'un texte 2) et un décalage d entier et calcule l'indice de coincidence mutuelle de ces 2 textes où le texte 2 est décalé de d positions
    """
    res = 0.0 #ICM
    n1 = 0 #len(texte1)
    n2 = 0 #len(texte2)

    for i in range(26):
        res += h1[i]*h2[(i+d)%26] # numérateur de la formule de l'ICM, texte2 décalé de d positions
       
        #MAJ de la taille
        n1 += h1[i]
        n2 += h2[i]

    return res / (n1*n2) # formule de l'ICM

# Renvoie le tableau des décalages probables étant
# donné la longueur de la clé
# en comparant l'indice de décalage mutuel par rapport
# à la première colonne
def tableau_decalages_ICM(cipher, key_length):
    """
Prend un texte chiffré (cipher) et la longueur supposée de sa clef et renvoie le tableau de décalages de chaque colonne par rapport à la colonne 0, décalage[i] maximise l'ICM entre la colonne 0 et la colonne i
    """
    decalages=[0]*key_length
    col0 = ''
    for i in range(key_length) :  # traitement de chacune des key_length colonnes
            coli = ''

            # on construit et stocke la colonne 0 dans col0
            if i == 0:
                for j in range(len(cipher)) :  
                    if j%key_length == i:              # on recupere chaque caractere du cipher dont l'indice est de la forme i+j*key_length
                        col0 += cipher[j]
       
                decalages[0] = 0
                freqCol0 = freq(col0) #on stocke les fréquences des caractères de la colonne 0
           
            #on parcourt toutes les autres colonnes
            else :
                for j in range(len(cipher)) :  #construction de la colonne i
                    if j%key_length == i:              # on recupere chaque caractere du cipher dont l'indice est de la forme i+j*key_length
                        coli += cipher[j]

                dmax = 0
                ICMmax = 0
                for d in range(26): #calcul du décalage d maximisant l'ICM entre la colonne 0 et la colonne i décalée de d
                    ICMactuel = indice_coincidence_mutuelle(freqCol0, freq(coli), d)
                    if ICMmax < ICMactuel:
                        dmax = d
                        ICMmax = ICMactuel
                decalages[i] = dmax #le décalage de la colonne i est dmax
    return decalages

# Cryptanalyse V2 avec décalages par ICM
def cryptanalyse_v2(cipher):
    """
     Prend un texte chiffre par vigenere (cipher) et renvoie un texte correspond a une tentative de dechiffrement avec décalages par ICM
    """
    key_length = longueur_clef(cipher) # longueur de la clef deduite avec l'indice de coincidence
    decalage = tableau_decalages_ICM(cipher, key_length) # tableau des decalages de chaque colonne par rapport a la premiere (en utilisant l'ICM)

    # construction du texte où chaque colonne a été décalée selon le tableau decalage afin d'obtenir un texte aligné avec la premiere colonne ; ce nouveau texte correspond ainsi a un texte chiffré par cesar
    texteCesar = ''
    i = 0
    for c in cipher:
        texteCesar += chr((ord(c) - decalage[i%key_length] - ord('A') ) %26 +ord('A') )
        i += 1
   
   
    res = dechiffre_cesar(texteCesar, (lettre_freq_max(texteCesar) - 4)%26) # on peut donc dechiffrer ce texte comme un texte de cesar : on decale toutes lettres avec le decalage correspondant a celui entre la lettre d'occurence maximale et 'E' (on l'identifie donc a 'E')
    return res

# Question 12

"""
43 textes (sur 100) sont correctement cryptanalysés avec cette fonction.
On observe les même problèmes que pour la fonction cryptanalyseV1 (lettre_freq_max et longeur_clef). En revanche, on a avec cette deuxième version (beaucoup) plus de textes correctement cryptanalysés.
Cette différence s'explique par le fait qu'on déchiffre avec César une seule fois et non plus key_length fois comme précédemment (donc un seul appel à lettre_freq_max, ce qui diminue le risque d'erreur lié à cette fonction).
De plus, le déchiffrement avec César est plus concluant lorsque le texte est représentatif de la langue (donc assez long). Ici, au lieu de l'appliquer à chaque colonne on l'applique sur un texte entier (donc key_length fois plus long).
"""

################################################################


### Les fonctions suivantes sont utiles uniquement
### pour la cryptanalyse V3.

# Prend deux listes de même taille et
# calcule la correlation lineaire de Pearson
def correlation(L1,L2):
    """
    Prend 2 listes L1 et L2 de même longueur et renvoie leur corrélation de Pearson
    len(L1) = len(L2)
    """

    #Calcul de l'espérance des variables L1 et L2
    moyL1 = 0
    moyL2 = 0
    n = len(L1)
    for i in range(n):
        moyL1 += L1[i]
        moyL2 += L2[i]
    moyL1 /= n
    moyL2 /= n

    numerateur = 0
    a = 0
    b = 0 # denominateur =  (a*b)**0.5

    for i in range(n):
        numerateur += (L1[i] - moyL1)*(L2[i] - moyL2)
        a += (L1[i] - moyL1)**2
        b += (L2[i] - moyL2) ** 2

    return numerateur / ((a*b)**0.5)

# Renvoie la meilleur clé possible par correlation
# étant donné une longueur de clé fixée
def clef_correlations(cipher, key_length):
    """
   Prend un texte chiffré (cipher) et la longueur potentielle de sa clef (key_length) et renvoie la moyenne sur les colonnes des corrélations obtenues avec pour chacune d'elles le décalage maximisant la corrélation de Pearson de la colonne associée avec la langue française ; elle retourne aussi la clef associée
    """
    key=[0]*key_length
    score = 0.0

    for i in range(key_length) :  # traitement de chacune des key_length colonnes
            coli = ''
            for j in range(len(cipher)) :  
                if j%key_length == i:              # on recupere chaque caractere du cipher dont l'indice est de la forme i+j*key_length
                    coli += cipher[j]

            dmax = 0
            corrMax = 0
            for d in range(26):  #calcul du décalage d maximisant la correlation de Pearson entre la langue française et la colonne i décalée de d
                coliDec = ''
                for c in coli: #construction de la colonne i décalée de d
                    coliDec += chr((ord(c) - d - ord('A') ) %26 +ord('A') )

                corrActuelle =  correlation(freq_FR, freq(coliDec))
               
                if corrMax < corrActuelle:
                    dmax = d
                    corrMax = corrActuelle

            score += corrMax
            key[i] = dmax #key[i] contient le décalage maximisant la corrélation

    score /= key_length #moyenne des corrélations maximales
    return (score, key)

# Cryptanalyse V3 avec correlations
def cryptanalyse_v3(cipher):
    """
    Prend un texte (cipher) chiffré par Vigenere et retourne son clair probable en utilisant la corrélation de Pearson pour déterminer la clef de chiffrement
    """
    scoreMax = 0
    decalageMax = []

    for t in range(1, 21): # parcours des tailles possibles de la clef
        score, decalage = clef_correlations(cipher, t) 

        # on garde la bonne clef dans decalageMax : c'est celle maximisant le score scoreMax
        if score>scoreMax:
            scoreMax = score
            decalageMax = decalage
    print(decalageMax)

    # on retourne le texte dechiffré avec cette clef 
    return dechiffre_vigenere(cipher, decalageMax) 

# Question 15

"""
94 textes (sur 100) sont correctement cryptanalysés. En effet, contrairement à cryptalyseV2 ou la cryptalyse est basée (après s'être ramené à un texte de César) sur le caractère le plus fréquent (E), la corrélation de pearson fait correspondre la distribution des freéquences de tous les caractères avec ceux de la langue française, ce qui explique le meilleur résultat obtenu.
Analysons les textes qui ne sont pas correctement cryptanalysés (il s'agit des textes 81, 86, 88, 89, 94, 96) :
    - On trouve dans chacun de ces cas une taille de clef correcte
    Texte 81 : longeur de la clef = 11 ; nombre d'erreurs obtenues : 1 / 11 ; la première colonne contenant l'erreur n'est pas représentative de la langue française (par exemple le A, le P et le U sont trop présents par rapport au E qui ne figure qu'une fois )
    Texte 86 : longeur de la clef = 14 ; nombre d'erreurs obtenues : 1 / 14 ; la dixième colonne contenant l'erreur n'est pas représentative de la langue française (par exemple elle contient X Y Z (3 lettres sur 15) qui ne sont pas récurentes dans la langue française)
    Texte 88 : longeur de la clef = 10 ; nombre d'erreurs obtenues : 1 / 10 ...
    Texte 89 : longeur de la clef = 9 ; nombre d'erreurs obtenues : 2 / 9
    Texte 94 : longeur de la clef = 18 ; nombre d'erreurs obtenues : 6 / 18 ; le texte est court et la clef est longue donc les colonnes sont trop courtes et pas représentatives de la langue française
    Texte 96 : longeur de la clef = 11 ; nombre d'erreurs obtenues : 2 / 11
    - Donc le probleme vient de la correspondance des colonnes avec la langue française : plus le texte est court et la clef est longue moins les colonnes du texte sont représentatives.

"""

################################################################
# NE PAS MODIFIER LES FONCTIONS SUIVANTES
# ELLES SONT UTILES POUR LES TEST D'EVALUATION
################################################################


# Lit un fichier et renvoie la chaine de caracteres
def read(fichier):
    f=open(fichier,"r")
    txt=(f.readlines())[0].rstrip('\n')
    f.close()
    return txt

# Execute la fonction cryptanalyse_vN où N est la version
def cryptanalyse(fichier, version):
    cipher = read(fichier)
    if version == 1:
        return cryptanalyse_v1(cipher)
    elif version == 2:
        return cryptanalyse_v2(cipher)
    elif version == 3:
        return cryptanalyse_v3(cipher)

def usage():
    print ("Usage: python3 cryptanalyse_vigenere.py -v <1,2,3> -f <FichierACryptanalyser>", file=sys.stderr)
    sys.exit(1)

def main(argv):
    size = -1
    version = 0
    fichier = ''
    try:
        opts, args = getopt.getopt(argv,"hv:f:")
    except getopt.GetoptError:
        usage()
    for opt, arg in opts:
        if opt == '-h':
            usage()
        elif opt in ("-v"):
            version = int(arg)
        elif opt in ("-f"):
            fichier = arg
    if fichier=='':
        usage()
    if not(version==1 or version==2 or version==3):
        usage()

    print("Cryptanalyse version "+str(version)+" du fichier "+fichier+" :")
    print(cryptanalyse(fichier, version))
   
if __name__ == "__main__":
   main(sys.argv[1:])
