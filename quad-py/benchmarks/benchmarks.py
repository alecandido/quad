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


def f3(x, delay):
    return (
        f1(x, delay),
        f2(x, delay),
        f1(x, delay),
        f1(x, delay),
        f2(x, delay),
        f1(x, delay),
    )


class TimeParam:
    timeout = 10000.0
    warmup_time = 3.0
    repeat = (1,10,1200.0)
    params = ([0.0, 1.0e-7, 1.0e-4], [100.0, 1000.0, 10000.0])
    param_names = ['delay', 'b']

    def setup(self, delay, b):
        self.number = 10
        self.a = 0.0
        self.limit = 1000000
        self.key = 2
        self.number_of_thread = 1
        self.f = lambda x: f3(x, delay)

    def time_qag_vec(self, delay, b):
        quad.qag_vec(
            self.f, self.a, b, limit=self.limit, key=self.key
        )

    def time_qag_par(self, dalay, b):
        quad.qag_par(
            self.f,
            self.a,
            b,
            limit=self.limit,
            key=self.key,
            number_of_thread=self.number_of_thread,
        )
