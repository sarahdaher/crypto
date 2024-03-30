from prime import is_probable_prime
from math import sqrt, ceil
import random


#Exercice 1
#Q1
def bezout(a, b):
    u0,v0, u1,v1, r0,r1 = 1,0,0,1, a, b
    while r1 != 0:
        qi = r0 // r1
        u0,u1 = u1, u0 - qi*u1
        v0,v1 = v1, v0 - qi*v1
        r0, r1 = r1, r0 - qi*r1
    return r0, u0, v0


#Q2
def inv_mod(a, n):
    pgcd, u, v = bezout(a,n)
    if pgcd == 1:
        return u%n
    else:
        return 9999 # pas d inverse


def invertibles(N):
    return [i for i in range(N) if inv_mod(i, N) != 9999]


#Q3
def phi(N):
    return len(invertibles(N)) #inv ds Z sont les premiers

#Exercice 2
#Q1
def exp1(a, n, p): #recursif
    if n == 0:
        return 1
    if n%2 == 0:
        ak2 = exp(a, n//2, p)
        return (ak2*ak2)%p
    else:
        ak2 = exp(a, n//2, p)
        return (a*ak2*ak2)%p


def exp(a, n, p):
    res = 1
    N = n
    A = a
    while N>0:
        if N%2 == 1: 
            res =(res*A)%p
        A = (A*A)%p 
        N = N//2
    return res

#Q2
def factor(n):
    N= n
    res = []
    p = 2
    while p<= sqrt(n):
        cpt = 0
        while N%p == 0:
            cpt+=1
            N = N/p 
        if cpt >0:
            res.append((p, cpt))
        p += 1
    if N>1: #N premier
        res.append((N, 1))
    return res


#Q3
def order(a, p, factors_p_minus1):
    orda = p - 1
    for (e, ve) in factors_p_minus1:
        orda = orda // e**ve
        A = exp(a, orda, p)
        for i in range(ve +1):
            if A == 1:
                break
            A = exp(A, e, p)
            orda *= e
    return orda


#Q4
def find_generator(p, factors_p_minus1):
    for g in range(p):
        if order(g,p,factors_p_minus1) == p-1:
            return g


#Q5
def generate_safe_prime(k):
    q = random.randint(2**(k-1), 2**k - 1)
    while not is_probable_prime(q) or not is_probable_prime(2*q+1): #7 -> 15 
        q = random.randint(2**(k-1), 2**k - 1)
    return 2*q + 1


#Q6
def bsgs(n, g, p):
    s = ceil((p)**0.5)+1

    BS = {exp(g, i, p): i for i in range(s)}
    g_inv= inv_mod(g,p)
    GS= {(n*exp(g_inv, s*i, p))%p : i for i in range(s)}
    
    for x in BS:
        if x in GS:
            return (BS[x] + s*GS[x])%p
    return None
