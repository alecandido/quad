#!/usr/bin/env python3
import quad
import math
import time
import cProfile
import numpy as np
from scipy.integrate import quad_vec
from scipy.integrate import quad as q


def bench_scipy():
    max = 100
    var = range(max)

    t1 = 0.0

    t2 = 0.0

    for a in var:
        f = lambda x: (math.cos(x), math.sin(x), math.cos(x), math.sin(x))
        epsabs = 1.0e-3
        epsrel = 0.0
        a = 0.0
        b = 10000.0

        #  start = time.time()
        #  res_array = quad.qag_array(f,a,b,epsabs,epsrel,6,10000,4)
        #  end = time.time()
        #  t1 += end - start
        #  print(t1)

        start = time.time()
        res_vec = quad.qag_vec(f, a, b, epsabs, epsrel, 6, 10000)
        end = time.time()
        t2 += end - start
        print(t2)

        if a > max - 10:
            print(res_array)
            print(res_vec)

    print(t1 / max)
    print(t2 / max)


# f = lambda x : np.array([math.cos(1/(x*x)),math.sin(1/(x*x))])
#
# start = time.time()
# y , err = quad_vec(f,a,b,epsabs=epsabs,epsrel=epsrel,limit=1000000)
# end = time.time()
# print(end - start)
# print(y)
# print(err)

# f = lambda x : math.cos(x)
# start = time.time()
# y , err = q(f,a,b,epsabs=epsabs,epsrel=epsrel,limit=1000000)
# f = lambda x : math.sin(x)
# y , err = q(f,a,b,epsabs=epsabs,epsrel=epsrel,limit=1000000)
# end = time.time()
# print(end - start)
# print(y)
# print(err)


# cProfile.run('quad.qag(f,a,b,epsabs,epsrel,6,8)','restats')
# cProfile.run("quad_vec(f,a,b,epsabs=epsabs,epsrel=epsrel,limit=1000000)",'restats2')
