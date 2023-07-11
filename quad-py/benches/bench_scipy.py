import quad
import math
import time
import timeit

import numpy as np
from scipy.integrate import quad_vec
from scipy.integrate import quad as q


def f1(x, delay):
    time.sleep(delay)
    return math.cos(x)


def f2(x, delay):
    time.sleep(delay)
    return math.sin(x)


def test_vs_scipy_quad_vec():
    number = 10
    a = 0.0
    b = 1000.0
    limit = 1000000
    delay = 1.0e-9
    key = 2

    def f(x):
        return (
            f1(x, delay),
            f2(x, delay),
            f1(x, delay),
            f1(x, delay),
            f2(x, delay),
            f1(x, delay),
        )

    def g(x):
        return np.array(
            [
                f1(x, delay),
                f2(x, delay),
                f1(x, delay),
                f1(x, delay),
                f2(x, delay),
                f1(x, delay),
            ]
        )

    def bench(test_call, flag: bool = False):
        times = timeit.repeat(test_call, number=number)
        print(times)
        res = test_call()
        if flag:
            return min(times), np.array(res[0]), res[1]
        return min(times), np.array(res.result), res.abserr

    quad_time, quad_res, quad_err = bench(
        lambda: quad.qag_par(f, a, b, limit=limit, key=key)
    )
    scipy_time, scipy_res, scipy_err = bench(
        lambda: quad_vec(g, a, b, limit=limit), True
    )

    print(f"\nratio: {scipy_time / quad_time * 100:.2f}%")
    print("distances:", (scipy_res - quad_res) / (scipy_err + quad_err))


def test_vs_scipy_quad():
    number = 1
    a = 0.0
    b = 100.0
    limit = 1000000
    delay = 1.0e-7
    # key is used only in quad.qag_vec not in q
    key = 2

    def f(x):
        return (
            f1(x, delay),
            f2(x, delay),
            f1(x, delay),
            f1(x, delay),
            f2(x, delay),
            f1(x, delay),
        )

    def g1(x):
        return f1(x, delay)

    def g2(x):
        return f2(x, delay)

    def bench(test_call, flag: bool = False):
        times = timeit.repeat(test_call, number=number)
        print(times)
        res = test_call()
        if flag:
            return (
                min(times),
                np.array(
                    [res[0][0], res[1][0], res[2][0], res[3][0], res[4][0], res[5][0]]
                ),
                math.sqrt(res[0][1] * res[0][1] + res[1][1] * res[1][1]),
            )

        return min(times), np.array(res.result), res.abserr

    quad_time, quad_res, quad_err = bench(
        lambda: quad.qag_par(f, a, b, limit=limit, key=key)
    )
    scipy_time, scipy_res, scipy_err = bench(
        lambda: (
            q(g1, a, b, limit=limit),
            q(g2, a, b, limit=limit),
            q(g1, a, b, limit=limit),
            q(g1, a, b, limit=limit),
            q(g2, a, b, limit=limit),
            q(g1, a, b, limit=limit),
        ),
        True,
    )

    print(f"\nratio: {scipy_time / quad_time * 100:.2f}%")
    print("distances:", (scipy_res - quad_res) / (scipy_err + quad_err))


