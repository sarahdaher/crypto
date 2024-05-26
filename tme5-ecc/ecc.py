# Sorbonne Université LU3IN024 2021-2022
# TME 5 : Cryptographie à base de courbes elliptiques
#
# Etudiant.e 1 : NOM ET NUMERO D'ETUDIANT : Sarah Daher 21100791
# Etudiant.e 2 : NOM ET NUMERO D'ETUDIANT : Maïssa Chemali 28722554

from math import sqrt
#import matplotlib.pyplot as plt
from random import randint

# Fonctions utiles

def exp(a, N, p):
    """Renvoie a**N % p par exponentiation rapide."""
    def binaire(N):
        L = list()
        while (N > 0):
            L.append(N % 2)
            N = N // 2
        L.reverse()
        return L
    res = 1
    for Ni in binaire(N):
        res = (res * res) % p
        if (Ni == 1):
            res = (res * a) % p
    return res


def factor(n):
    """ Return the list of couples (p, a_p) where p is a prime divisor of n and
    a_p is the p-adic valuation of n. """
    def factor_gen(n):
        j = 2
        while n > 1:
            for i in range(j, int(sqrt(n)) + 1):
                if n % i == 0:
                    n //= i
                    j = i
                    yield i
                    break
            else:
                if n > 1:
                    yield n
                    break

    factors_with_multiplicity = list(factor_gen(n))
    factors_set = set(factors_with_multiplicity)

    return [(p, factors_with_multiplicity.count(p)) for p in factors_set]


def inv_mod(x, p):
    """Renvoie l'inverse de x modulo p."""
    return exp(x, p-2, p)


def racine_carree(a, p):
    """Renvoie une racine carrée de a mod p si p = 3 mod 4."""
    assert p % 4 == 3, "erreur: p != 3 mod 4"

    return exp(a, (p + 1) // 4, p)


# Fonctions demandées dans le TME

def est_elliptique(E):
    """
    Renvoie True si la courbe E est elliptique et False sinon.

    E : un triplet (p, a, b) représentant la courbe d'équation
    y^2 = x^3 + ax + b sur F_p, p > 3
    """
    p,a,b = E

    return ((-16*(4*exp(a,3,p) + 27*exp(b,2,p))) % p) != 0


def point_sur_courbe(P, E):
    """Renvoie True si le point P appartient à la courbe E et False sinon."""
    if P == () :
        return True
    x,y=P
    p,a,b = E
    return exp(y,2,p) == (exp(x,3,p) + a*x +b)%p


"""
Montrons que le symbole de Legendre est donné par a^((p-1)/2) mod p
- Cas 1 : a est un résidu quadratique, dans ce cas il existe b tel que a = b^2 mod p donc a^((p-1)/2) = b^(p-1) mod p 
Par application du petit théorème de Fermat on a que b^(p-1) = 1 mod p et comme 1 est bien égale au symbole de Legendre dans ce cas alors la formule est vérifiée.
- Cas 2 : p divise a donc a = 0 mod p alors a^((p-1)/2) = 0 mod p. Comme dans ce cas 0 est bien égal au symbole de Legendre alors la formule est vérifiée.
- Cas 3 : a n'est pas un résidu quadratique, donc il n'existe pas de b tel que a = b^2 mod p.
Considérons le polynome x^((p-1)/2) - 1, comme l'anneau (Z/pZ)* est integre, ce polynome admet (p-1)/2 racines qui sont les résidus quadratiques.
Or a^(p-1) = 1 mod p par Fermat mais a n'est pas un résidu quadratique donc a^((p-1)/2) != 1.
Ainsi, nous avons forcèment a^((p-1)/2) = -1 (pour vérifier les 2 conditions précédentes) et -1 est bien égal dans ce cas au symbole de Legendre.
"""
def symbole_legendre(a, p):
    """Renvoie le symbole de Legendre de a mod p."""
    
    return exp(a, (p-1)//2, p)


def cardinal(E):
    """Renvoie le cardinal du groupe de points de la courbe E."""
    p,a,b = E
    cpt = 1 # pt a l'infini
    for x in range(p):
        z = (exp(x,3,p) + a*x + b) % p  
        legendre = symbole_legendre(z, p)
        if legendre == 1: 
            cpt += 2     
        elif legendre == 0:
            cpt += 1 
    return cpt
"""
Montrons que lorsque p premier est congru a 3 mod 4 et que a est un carre, on a que a**((p+1)/4) est une racine carree de a.
Calculons ( a**((p+1)/4) )**2 mod p et montrons que cela vaut a.
p congru a 3 mod 4, donc il existe k positif tel que p = 4k+3. 
Alors : ( a**((p+1)/4) )**2 = ( a**((4k+4)/4) )**2 = ( a**(k+1) )**2 = a**(2k+2) = a * a**(2k+1)
Il suffit donc de montrer que a**(2k+1) mod p vaut 1.
On a p-1 = 2*(2k+1).
Donc a**(2k+1) = a**((p-1)/2) = 1 mod p par la remarque de la question 3, si a possede une racine.
On a bien montre ce que l'on voulait.
"""

def liste_points(E):
    """Renvoie la liste des points de la courbe elliptique E."""
    p, a, b = E

    assert p % 4 == 3, "erreur: p n'est pas congru à 3 mod 4."

    lstpoints = [()]
    dic = dict()

    
    for x in range(p) :
        res = (exp(x,3,p) + a*x +b)%p
        
        if symbole_legendre(res,p)==0: 
            lstpoints.append((x,0))

        elif symbole_legendre(res,p)==1 : 
            elem = exp(res,(p+1)//4,p)
            if elem == 0 :
                lstpoints.append((x,0))
            else :
                lstpoints.append((x,elem))
                lstpoints.append((x,p-elem))
    
    return lstpoints

"""
Théorème de Hasse : Soit K = F_p avec p > 3 un nombre premier,
p+1-2*sqrt(p) <= card(E) <= p+1+2*sqrt(p)
"""
def cardinaux_courbes(p):
    """
    Renvoie la distribution des cardinaux des courbes elliptiques définies sur F_p.

    Renvoie un dictionnaire D où D[i] contient le nombre de courbes elliptiques
    de cardinal i sur F_p.
    """
    D = {}
    pp = int(p+1-2*sqrt(p))
    pg = int(p+1+2*sqrt(p))

    for i in range(pp+1, pg+1): #initialisation
        D[i] = 0

    for a in range(p):
        for b in range(p):
            E = p,a,b #testes toutes les valeurs possibles de a et b
            if est_elliptique(E):
                card = cardinal(E)
                D[card] += 1

    return D


def dessine_graphe(p):
    """Dessine le graphe de répartition des cardinaux des courbes elliptiques définies sur F_p."""
    bound = int(2 * sqrt(p))
    C = [c for c in range(p + 1 - bound, p + 1 + bound + 1)]
    D = cardinaux_courbes(p)

    plt.bar(C, [D[c] for c in C], color='b')
    plt.show()


def moins(P, p):
    """Retourne l'opposé du point P mod p."""
    if P==():
        return ()
    x,y = P
    return x, (p-y)%p


def est_egal(P1, P2, p):
    """Teste l'égalité de deux points mod p."""

    if P1 == P2:
        return True
    
    elif P1 != () and P2 != ():
        x1,y1 = P1
        x2, y2 = P2
        return (x1%p == x2%p) and (y1%p == y2%p)
    
    return False


def est_zero(P):
    """Teste si un point est égal au point à l'infini."""

    return P == ()


def addition(P1, P2, E):
    """Renvoie P1 + P2 sur la courbe E."""
    p,a,b = E

    if est_zero(P2):
        return P1
    elif est_zero(P1):
        return P2
    
    elif est_egal(P1, moins(P2, p), p) or est_egal(P2, moins(P1, p), p):
        return()
    
    x1,y1 = P1
    x2,y2 = P2

    if not est_egal(P1, P2, p):
        lamda = ((y2-y1)*inv_mod((x2-x1), p))%p
    else:
        lamda = ((3*x1**2+a)*inv_mod(2*y1, p))%p

    xres = (lamda**2 - x1 - x2)%p
    yres = lamda*(x1 - xres) - y1

    if point_sur_courbe((xres,yres), E):
        return (xres,yres)
    else:
        return ()

def multiplication_scalaire(k, P, E):
    """Renvoie la multiplication scalaire k*P sur la courbe E."""
    
    def binaire(N):
        L = list()
        while (N > 0):
            L.append(N % 2)
            N = N // 2
        L.reverse()
        return L
    
    p, a, b = E

    res = ()
    for Ni in binaire(abs(k)):
        res = addition(res, res, E)
        if (Ni == 1):
            res = addition(res, P, E)

    if k>=0:
        return res
    return moins(res, p)
    

def ordre(N, factors_N, P, E):
    """Renvoie l'ordre du point P dans les points de la courbe E mod p. 
    N est le nombre de points de E sur Fp.
    factors_N est la factorisation de N en produit de facteurs premiers."""

    if P==() :
        return 1
    produits = [1]
    for q,vq in factors_N :
        nv_produits = []
        for prod in produits :
            elem=prod
            res = multiplication_scalaire(prod, P, E)
            if res==() :
                return prod
            nv_produits.append(prod)
            for k in range(1,vq+1) :
                elem = elem *q
                res = multiplication_scalaire(q, res, E)
                if res == () :
                    return elem
                nv_produits.append(elem)
        produits = nv_produits
   

    
"""
On sait par le théorème de Hasse que le cardinal de E est de l'ordre de p. On a p² couples possibles sachant que chacun peut etre tire plusieurs fois.
La probabilité de tirer chacun de ces couples suit une loi uniforme donc = 1/p², ainsi la probabilite de tirer un point sur la courbe est egal p/p² = 1/p.
Donc il faudrait tirer environ p points pour en obtenir un sur la courbe.
Donc la complexité de la fonction est en O(p) (car randint est en O(1)). 
Donc plus p est grand plus la fonction prend de temps.
"""

def point_aleatoire_naif(E):
    """Renvoie un point aléatoire (différent du point à l'infini) sur la courbe E."""
    p,_,_ = E
    x = randint(0,p-1)
    y= randint(0,p-1)
    while not point_sur_courbe((x,y),E) :
        x = randint(0,p-1)
        y= randint(0,p-1)
    return (x,y)

#print(point_aleatoire_naif((360040014289779780338359, 117235701958358085919867, 18575864837248358617992)))

""" 
Le test prend énormément de temps car est en O(360040014289779780338359)
"""


def point_aleatoire(E):
    """Renvoie un point aléatoire (différent du point à l'infini) sur la courbe E."""
    p, a, b = E
    x = randint(0, p-1)
    y_carre = (exp(x,3,p) + a*x + b) % p #complexité : O(log(p)² * log(3))
    y = racine_carree(y_carre, p)  #complexité : O(log(p)³)
    return (x, y)

#print(point_aleatoire((360040014289779780338359, 117235701958358085919867, 18575864837248358617992)))

""" La fonction point_aléatoire est en O(log(p)³) qui est de meilleure complexité que O(p). Cette fois le test est plus rapide et se termine"""


def point_ordre(E, N, factors_N, n):
    """Renvoie un point aléatoire d'ordre N sur la courbe E.
    Ne vérifie pas que n divise N."""
    P = point_aleatoire(E)
    while ordre(N,factors_N,P,E) !=n :
        P = point_aleatoire(E)
    return P

def keygen_DH(P, E, n):
    """Génère une clé publique et une clé privée pour un échange Diffie-Hellman.
    P est un point d'ordre n sur la courbe E.
    """
    sec = randint(1,n-1)
    pub = multiplication_scalaire(sec, P, E) 
    
    return (sec, pub)

def echange_DH(sec_A, pub_B, E):
    """Renvoie la clé commune à l'issue d'un échange Diffie-Hellman.
    sec_A est l'entier secret d'Alice et pub_b est l'entier public de Bob."""

    return multiplication_scalaire(sec_A, pub_B, E)

"""Pour un echange de clefs Diffie Hellman, nous devons avoir une clef telle que son ordre divise N (Lagrange) et soit le plus grand possible pour avoir un générateur qui est difficile à casser.
Regardons si on a un générateur d'ordre N :
"""
sec, pub = point_ordre((248301763022729027652019747568375012323, 1, 0), 248301763022729027652019747568375012324, [(2, 2), (62075440755682256913004936892093753081, 1)], 2**2 * 62075440755682256913004936892093753081)
"""On voit ainsi que ce point existe. Donc notre groupe est cyclique et P = (sec, pub) en est un generateur qui est optimal pour l'echange DH."""
