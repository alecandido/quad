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


def f3(x, delay):
    return (
        f1(x, delay),
        f2(x, delay),
        f1(x, delay),
        f1(x, delay),
        f2(x, delay),
        f1(x, delay),
    )

def f_scipy(x, delay):
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

class QagVecVsQagPar:
    timeout = 10000.0
    warmup_time = 3.0
    repeat = (1,10,1200.0)
    params = ([0.0, 1.0e-7, 1.0e-4], [100.0, 1000.0, 10000.0])
    param_names = ['delay', 'b']

    def setup(self, delay, b):
        self.a = 0.0
        self.limit = 1000000
        self.key = 2
        self.number_of_thread = 1
        self.f = lambda x: f3(x, delay)
        self.f_scipy = lambda x: f_scipy(x, delay)
        self.f1_scipy = lambda x: f1(x, delay)
        self.f2_scipy = lambda x: f2(x, delay)

    def time_qag_par(self, delay, b):
        quad.qag_par(
            self.f,
            self.a,
            b,
            limit=self.limit,
            key=self.key,
            number_of_thread=self.number_of_thread,
        )

    def time_scipy_vec(self, delay, b):
        quad_vec(self.f_scipy, self.a, b, limit=self.limit)

    def time_scipy_1d(self, delay, b):
        q(self.f1_scipy, self.a, b, limit=self.limit)
        q(self.f2_scipy, self.a, b, limit=self.limit)
        q(self.f1_scipy, self.a, b, limit=self.limit)
        q(self.f1_scipy, self.a, b, limit=self.limit)
        q(self.f2_scipy, self.a, b, limit=self.limit)
        q(self.f1_scipy, self.a, b, limit=self.limit)
