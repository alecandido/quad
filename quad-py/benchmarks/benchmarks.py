import timeit
import numpy as np
import quad
import time
import math
from scipy.integrate import quad_vec


def f1(x, delay):
    time.sleep(delay)
    return math.cos(x)


def f2(x, delay):
    time.sleep(delay)
    return math.sin(x)


def bench(test_call, number, flag: bool = False):
    times = timeit.repeat(test_call, number=number)
    res = test_call()
    if flag:
        return min(times), np.array(res[0]), res[1]
    return min(times), np.array(res.result), res.abserr


class TimeSuite:
    def setup(self):
        self.number = 10
        self.a = 0.0
        self.b = 1000.0
        self.limit = 1000000
        self.delay = 0.0
        self.key = 2

        def f(x):
            return (
                f1(x, self.delay),
                f2(x, self.delay),
                f1(x, self.delay),
                f1(x, self.delay),
                f2(x, self.delay),
                f1(x, self.delay),
            )

        def g(x):
            return np.array(
                [
                    f1(x, self.delay),
                    f2(x, self.delay),
                    f1(x, self.delay),
                    f1(x, self.delay),
                    f2(x, self.delay),
                    f1(x, self.delay),
                ]
            )

        self.f = lambda x: f(x)
        self.g = lambda x: g(x)

    def time_qag_vec(self):
        bench(
            lambda: quad.qag_vec(
                self.f, self.a, self.b, limit=self.limit, key=self.key
            ),
            self.number,
        )

    def time_scipy_quad_vec(self):
        bench(
            lambda: quad_vec(self.g, self.a, self.b, limit=self.limit),
            self.number,
            True,
        )
