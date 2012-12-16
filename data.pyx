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
    Return the list of pairs (p, D_E(p)), for all primes p < B, and
    another list of pairs (p, running_average).

    E -- an elliptic curve over QQ
    B -- positive integer
    """
    result = []; avgs = []
    aplist = E.aplist(B)
    primes = prime_range(B)
    assert len(aplist) == len(primes)
    cdef double X, running_sum = 0, val = 0
    cdef int i, cnt = 0
    for i in range(len(aplist)):
        if aplist[i] < 0:
            cnt -= 1
        elif aplist[i] > 0:
            cnt += 1

        if i == 0:
            running_sum = 0
        else:
            running_sum += val*(primes[i] - primes[i-1] - 1)

        X   = primes[i]
        val = (log(X) / sqrt(X)) * cnt
        running_sum += val
        result.append((primes[i], val))
        avgs.append((primes[i], running_sum/X))
    return result, avgs

def the_mean(r, vanishing_symmetric_powers):
    """
    INPUT:

    """
    mean = 2/pi - 16/(3*pi)*r
    for n, order in vanishing_symmetric_powers:
        assert n%2 != 0
        k = (n-1)/2
        mean += (4/pi) * (-1)**(k+1)*(1/(2*k+1) + 1/(2*k+3))*order
    return mean

def draw_plot(E, B, vanishing_symmetric_powers=None):
    if vanishing_symmetric_powers is None:
        vanishing_symmetric_powers = []
    mean = the_mean(E.rank(), vanishing_symmetric_powers)
    d, running_average = raw_data(E, B)
    g = line(d)
    g += line([(0,mean), (d[-1][0],mean)], color='darkred')
    g += line(running_average, color='green')
    return g

