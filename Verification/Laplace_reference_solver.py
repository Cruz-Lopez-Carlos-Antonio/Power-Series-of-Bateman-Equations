import itertools
import numpy as np
import math
from decimal import Decimal, getcontext

getcontext().prec = 32


def partitions_restricted(m, n, L):
    # Generates all n-tuples of nonnegative integers summing to m
    values = list(range(m + 1))
    for k in itertools.product(values, repeat=n):
        a = np.array(k, dtype=int)
        if np.sum(a) == m:
            L.append(list(a))


def chi(i, j, Mu, Lambd, L):  # 'j' is kept for compatibility; not used internally
    c = Decimal(0)

    # Exclude by index (NOT by value)
    Aux_L = Mu[:i] + Mu[i+1:]
    Aux2  = Lambd[:i] + Lambd[i+1:]

    for u in L:
        a = Decimal(1)
        for k in range(len(Aux_L)):
            b_f = Decimal(math.comb(Aux_L[k] + u[k], Aux_L[k]))
            dif = Decimal(1) / ( (Decimal(Lambd[i]) - Decimal(Aux2[k])) ** Decimal(int(u[k])) )
            a = a * b_f * dif
        c = c + a

    return c


def GPS(X0, DC, t):
    term1 = Decimal(X0) / Decimal(DC[-1])

    # Unique values preserving first-appearance order (deterministic)
    Au1 = list(dict.fromkeys(DC))

    Mu = [DC.count(u) - 1 for u in Au1]

    term2 = Decimal(1)
    for z in range(len(Au1)):
        term2 = term2 * (Decimal(Au1[z]) ** Decimal(Mu[z] + 1))

    s1 = Decimal(0)

    for i in range(len(Au1)):
        p1 = Decimal.exp(-Decimal(Au1[i]) * Decimal(t))

        # Exclude by index (NOT by value)
        Au2  = Au1[:i] + Au1[i+1:]
        Mu_m = Mu[:i] + Mu[i+1:]

        p2 = Decimal(1)
        for j in range(len(Au2)):
            p2 = p2 * (Decimal(1) / ((Decimal(Au2[j]) - Decimal(Au1[i])) ** Decimal(Mu_m[j] + 1)))

        s2 = Decimal(0)
        for l in range(Mu[i] + 1):
            z1 = Decimal(t) ** Decimal(l) / Decimal(math.factorial(l))

            L = []
            partitions_restricted(Mu[i] - l, len(Au2), L)

            z2 = chi(i, Mu[i] - l, Mu, Au1, L)
            s2 = s2 + z1 * z2

        s1 = s1 + p1 * p2 * s2

    solution = term1 * term2 * s1
    return solution


Half_lifes = [2,2,3,3,3,4,4,4,4]
D = [Decimal(math.log(2)) / Decimal(z) for z in Half_lifes]

Time_vector = [1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100]

for k in Time_vector:
    solucion = GPS(6.023E23, D, k)
    print(float(solucion))
