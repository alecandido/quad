import unittest
import quad
import math
import pytest
import numpy


class MyTestCase(unittest.TestCase):
    def test_invalid(self):
        a = 0.0
        b = 1000.0
        limit = 100000
        key = 6
        epsabs = 0.0
        epsrel = 1.0e-40
        f = lambda x: (math.sin(x), math.cos(x))
        with pytest.raises(Exception) as exc_info:
            quad.qag_vec(f, a, b, epsabs, epsrel, key, limit)
        exception_raised = exc_info.value
        error_message = "The input is invalid, because epsabs <= 0 and epsrel < max(50 * rel.mach.acc.,0.5d-28)"
        assert exception_raised.args[0] == error_message

    def test_max_iteration1(self):
        a = 0.0
        b = 1000.0
        limit = 1
        key = 6
        epsabs = 0.0
        epsrel = 1.0e-3
        f = lambda x: (math.sin(x), math.cos(x))
        with pytest.raises(Exception) as exc_info:
            quad.qag_vec(f, a, b, epsabs, epsrel, key, limit)
        exception_raised = exc_info.value
        error_message = "Maximum number of subdivisions allowed has been achieved. One can allow more subdivisions by increasing the value of limit. However, if this yields no improvement it is rather advised to analyze the integrand in order to determine the integration difficulties. If the position of a local difficulty can be determined(e.g. singularity, discontinuity within the interval) one will probably gain from splitting up the interval at this point and calling the integrator on the subranges. If possible, an appropriate special-purpose integrator should be used which is designed for handling the type of difficulty involved."
        assert exception_raised.args[0] == error_message

    def test_max_iteration2(self):
        a = 0.0
        b = 1000000.0
        limit = 30
        key = 6
        epsabs = 0.0
        epsrel = 1.0e-3
        f = lambda x: (math.sin(x), math.cos(x))
        with pytest.raises(Exception) as exc_info:
            quad.qag_vec(f, a, b, epsabs, epsrel, key, limit)
        exception_raised = exc_info.value
        error_message = "Maximum number of subdivisions allowed has been achieved. One can allow more subdivisions by increasing the value of limit. However, if this yields no improvement it is rather advised to analyze the integrand in order to determine the integration difficulties. If the position of a local difficulty can be determined(e.g. singularity, discontinuity within the interval) one will probably gain from splitting up the interval at this point and calling the integrator on the subranges. If possible, an appropriate special-purpose integrator should be used which is designed for handling the type of difficulty involved."
        assert exception_raised.args[0] == error_message

    def test_key(self):
        a = 0.0
        b = 10000.0
        limit = 10000
        epsabs = 1.0e-3
        epsrel = 0.0
        correct_result = (1.0 - math.cos(10000.0), math.sin(10000.0))
        f = lambda x: (math.sin(x), math.cos(x))

        for key in range(1, 6):
            res = quad.qag_vec(f, a, b, epsabs, epsrel, key, limit)
            assert res.result[0] - correct_result[0] < epsabs and res.result[1] - correct_result[1] < epsabs

    def test_semi_infinite(self):
        a = 0.0
        b = math.inf
        c = - math.inf
        limit = 10000
        epsabs = 1.0e-12
        epsrel = 0.0
        key = 6
        correct_result = (0.4, 0.6)

        f = lambda x: (math.sin(x) * math.sin(x) / math.exp(abs(x)), math.cos(x) * math.cos(x) / math.exp(abs(x))) \
            if math.fabs(x) <= 300.0 else (0.0, 0.0)

        res1 = quad.qag_vec(f, a, b, epsabs, epsrel, key, limit)
        res2 = quad.qag_vec(f, c, a, epsabs, epsrel, key, limit)

        assert res1.result[0] - correct_result[0] < epsabs and res1.result[1] - correct_result[1] < epsabs
        assert res2.result[0] - correct_result[0] < epsabs and res2.result[1] - correct_result[1] < epsabs

    def test_double_infinite(self):
        a = - math.inf
        b = math.inf
        limit = 10000
        epsabs = 1.0e-10
        epsrel = 0.0
        key = 6
        correct_result = (1.2879903316984565533522585284072106913, 1.5974)

        f = lambda x: (math.sin(x) * math.sin(x) / numpy.exp2(abs(x)), math.cos(x) * math.cos(x) / numpy.exp2(abs(x))) \
            if math.fabs(x) <= 300.0 else (0.0, 0.0)

        res = quad.qag_vec(f, a, b, epsabs, epsrel, key, limit)

        assert res.result[0] - correct_result[0] < epsabs and res.result[1] - correct_result[1] < epsabs

    def test_additional_points(self):
        a = 0.0
        b = 1.0
        limit = 10000
        epsabs = 1.0
        epsrel = 0.0
        key = 6
        points = (0.0,0.2,0.4,0.6,0.8,1.0)
        sub_interval_expected = [(0.0,0.2),(0.2,0.4),(0.4,0.6),(0.6,0.8),(0.8,1.0)]
        more_info = True

        f = lambda x: (math.cos(x),math.sin(x))

        res = quad.qag_vec(f, a, b, epsabs, epsrel, key, limit,points,more_info)
        alist = []
        blist = []
        for i in range(5):
            alist.append(res.more_info[2][i][0])
            blist.append(res.more_info[2][i][1])
        alist.sort()
        blist.sort()
        sub_interval = []
        for i in range(5):
            sub_interval.append((alist[i],blist[i]))

        assert sub_interval == sub_interval_expected


if __name__ == '__main__':
    unittest.main()
