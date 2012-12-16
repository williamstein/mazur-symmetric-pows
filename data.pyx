#########################################################
# Efficiently compute raw sum for a given elliptic curve
# (c) William Stein
#########################################################
from sage.all import prime_range

cdef extern from "math.h":
    double log(double)
    double sqrt(double)

def raw_data(E, B):
    """
    Return the list of pairs (p, D_E(p)), for all primes p < B.

    E -- an elliptic curve over QQ
    B -- positive integer
    """
    result = []
    aplist = E.aplist(B)
    primes = prime_range(B)
    assert len(aplist) == len(primes)
    cdef double X
    cdef int i, cnt = 0
    for i in range(len(aplist)):
        if aplist[i] < 0:
            cnt -= 1
        elif aplist[i] > 0:
            cnt += 1
        X = primes[i]
        result.append((primes[i], (log(X) / sqrt(X)) * cnt))
    return result
