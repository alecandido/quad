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

def f3(x,delay):
    return (
    f1(x, delay),
    f2(x, delay),
    f1(x, delay),
    f1(x, delay),
    f2(x, delay),
    f1(x, delay),
)


def bench(test_call, number, flag: bool = False):
    times = timeit.repeat(test_call, number=number)
    res = test_call()
    if flag:
        return min(times), np.array(res[0]), res[1]
    return min(times), np.array(res.result), res.abserr


class TimeSuite:
    timeout = 1800.0

    def setup(self):
        self.number = 10
        self.a = 0.0
        self.b = 1000.0
        self.limit = 1000000
        self.delay = 0.0
        self.key = 2
        self.number_of_thread = 1
        self.h = lambda x,delay: (
            f1(x, delay),
            f2(x, delay),
            f1(x, delay),
            f1(x, delay),
            f2(x, delay),
            f1(x, delay),
        )

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

    def time_qag_par(self):
        bench(
            lambda: quad.qag_par(
                self.f,
                self.a,
                self.b,
                limit=self.limit,
                key=self.key,
                number_of_thread=self.number_of_thread,
            ),
            self.number,
        )

    def time_scipy_quad_vec(self):
        bench(
            lambda: quad_vec(self.g, self.a, self.b, limit=self.limit),
            self.number,
            True,
        )

class TimeParam:
    timeout = 10000.0
    params = [0.0,1.0e-9,1.0e-8,1.0e-7]
    param_names = ['delay']

    def setup(self,delay):
        self.number = 10
        self.a = 0.0
        self.b = 1000.0
        self.limit = 1000000
        self.key = 2
        self.number_of_thread = 1
        self.f = lambda x: f3(x,delay)

    def time_qag_vec(self,delay):
        bench(
            lambda: quad.qag_vec(
                self.f, self.a, self.b, limit=self.limit, key=self.key
            ),
            self.number,
        )