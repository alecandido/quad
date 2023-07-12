import numpy as np
import quad
import time
import math
from scipy.integrate import quad_vec
from scipy.integrate import quad as q


def f1(x, delay):
    time.sleep(delay)
    return math.cos(x)


def f2(x, delay):
    time.sleep(delay)
    return math.sin(x)


def f_ln6(x, delay, i):
    if i == 0:
        return (
            f1(x, delay),
            f2(x, delay),
            f1(x, delay),
            f1(x, delay),
            f2(x, delay),
            f1(x, delay),
        )
    elif i == 1:
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
    else:
        return [[f1(x, delay), 4], [f2(x, delay), 2]]


def f_ln2(x, delay, i):
    if i == 0:
        return (
            f1(x, delay),
            f2(x, delay),
        )
    elif i == 1:
        return np.array(
            [
                f1(x, delay),
                f2(x, delay),
            ]
        )
    else:
        return [[f1(x, delay), 1], [f2(x, delay), 1]]


class QagBench:
    timeout = 10000.0
    warmup_time = 5.0
    repeat = (1, 10, 1200.0)
    params = ([0.0, 1.0e-7, 1.0e-4], [100.0, 1000.0, 10000.0], [f_ln6, f_ln2])
    param_names = ["delay", "b", "fun"]

    def setup(self, delay, b, fun):
        self.a = 0.0
        self.limit = 1000000
        self.key = 2
        self.epsabs = 1.0e-6
        self.epsrel = 0.0

    def time_qag(self, delay, b, fun):
        f = lambda x: fun(x, delay, 0)
        quad.qag(
            f,
            self.a,
            b,
            limit=self.limit,
            key=self.key,
            epsabs=self.epsabs,
            epsrel=self.epsrel,
        )

    def time_scipy_vec(self, delay, b, fun):
        f = lambda x: fun(x, delay, 1)
        quad_vec(f, self.a, b, limit=self.limit, epsabs=self.epsabs, epsrel=self.epsrel)

    def time_scipy_1d(self, delay, b, fun):
        f_comp1 = lambda x: fun(x, delay, 2)[0][0]
        f_comp2 = lambda x: fun(x, delay, 2)[1][0]
        k1 = fun(1.0, delay, 2)[0][1]
        k2 = fun(1.0, delay, 2)[1][1]
        for k in range(k1):
            q(
                f_comp1,
                self.a,
                b,
                limit=self.limit,
                epsabs=self.epsabs,
                epsrel=self.epsrel,
            )
        for k in range(k2):
            q(
                f_comp2,
                self.a,
                b,
                limit=self.limit,
                epsabs=self.epsabs,
                epsrel=self.epsrel,
            )
