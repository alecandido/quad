import unittest
import quad
import math


class MyTestCase(unittest.TestCase):
    def test_max_iteration1(self):
        a = 0.0
        b = 1000.0
        limit = 100000000
        key = 2
        epsabs = 1.0
        epsrel = 1.0
        f = lambda x: (math.cos(x), math.sin(x),)
        print(quad.qag_vec(f, a, b, epsabs, epsrel,key,limit))
        self.assertRaises(Exception, quad.qag_vec, (f, a, b, epsabs, epsrel,key,limit))


if __name__ == '__main__':
    unittest.main()
