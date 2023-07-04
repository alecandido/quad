import unittest
import quad
import math
import pytest


class MyTestCase(unittest.TestCase):
    def test_max_iteration1(self):
        a = 0.0
        b = 1000.0
        limit = 1
        key = 2
        epsabs = 1.0
        epsrel = 1.0
        f = lambda x: (math.cos(x), math.sin(x),)
        with pytest.raises(TypeError) as exc_info:
            quad.qag_vec(f, a, b, epsabs, epsrel,key,limit)

        exception_raised = exc_info.value
        print(exception_raised)


if __name__ == '__main__':
    unittest.main()
