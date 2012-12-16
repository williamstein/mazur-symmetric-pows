#########################################################
# Efficiently compute raw sum for a given elliptic curve
# (c) William Stein
#########################################################
from sage.all import prime_range, line
import math

cdef extern from "math.h":
    double log(double)
    double sqrt(double)

cdef double pi = math.pi

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

def draw_plot(E, B, vanishing_symmetric_powers=None):
    if vanishing_symmetric_powers is None:
        vanishing_symmetric_powers = []
    r = E.rank()
    mean = 2/pi - 16/(3*pi)*r
    for n, order in vanishing_symmetric_powers:
        assert n%2 != 0
        k = (n-1)/2
        mean += (4/pi) * (-1)**(k+1)*(1/(2*k+1) + 1/(2*k+3))*order
    d = raw_data(E, B)
    g = line(d)
    g += line([(0,mean), (d[-1][0],mean)], color='darkred')
    return g
