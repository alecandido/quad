import quad
import math
import time
import timeit

import numpy as np


def test_parallel():
    max = 3
    var = range(max)
    a = 0.0
    b = 10000.0
    limit = 100000

    for k in var:
        f = lambda x: (math.cos(x), math.sin(x))

        start = time.time()
        res = quad.qag_par(f, a, b, limit=limit, more_info=True)
        end = time.time()
        print("par time is", end - start)
        print("par res : ", res.result)
        print("par err : ", res.abserr)
        print("par last : ", res.more_info[1])

        start = time.time()
        res = quad.qag_nopar(f, a, b, limit=limit, more_info=True)
        end = time.time()
        print("serial time is", end - start)
        print("serial res : ", res.result)
        print("serial err : ", res.abserr)
        print("serial last : ", res.more_info[1])


def test_parallel_2():
    number = 100
    a = 0.0
    b = 10000.0
    limit = 1000000

    f = lambda x: (
        math.cos(x),
        math.sin(x),
        math.cos(x),
        math.cos(x),
        math.sin(x),
        math.cos(x),
    )

    def bench(test_call):
        times = timeit.repeat(test_call, number=number)
        print(times)
        res = test_call()
        return min(times), np.array(res.result), res.abserr

    stime, sres, serr = bench(lambda: quad.qag(f, a, b, limit=limit))
    ptime, pres, perr = bench(lambda: quad.qag_par(f, a, b, limit=limit))

    print(f"\nratio: {ptime / stime * 100:.2f}%")
    print("distances:", (pres - sres) / (perr + serr))
