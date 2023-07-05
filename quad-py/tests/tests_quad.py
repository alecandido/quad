import unittest
import quad
import math
import pytest


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
            quad.qag_vec(f, a, b, epsabs, epsrel,key,limit)
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
            quad.qag_vec(f, a, b, epsabs, epsrel,key,limit)
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
            quad.qag_vec(f, a, b, epsabs, epsrel,key,limit)
        exception_raised = exc_info.value
        error_message = "Maximum number of subdivisions allowed has been achieved. One can allow more subdivisions by increasing the value of limit. However, if this yields no improvement it is rather advised to analyze the integrand in order to determine the integration difficulties. If the position of a local difficulty can be determined(e.g. singularity, discontinuity within the interval) one will probably gain from splitting up the interval at this point and calling the integrator on the subranges. If possible, an appropriate special-purpose integrator should be used which is designed for handling the type of difficulty involved."
        assert exception_raised.args[0] == error_message



if __name__ == '__main__':
    unittest.main()
