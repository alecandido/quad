import quad
import math
import time



max = 1
var = range(max)
a = 0.0
b = 10000.0
limit = 1000000

for k in var:

    f = lambda x: (math.cos(x),math.sin(x),math.cos(x),math.cos(x),math.sin(x),math.cos(x),)


    start = time.time()
    res = quad.qag_par(f,a,b,limit = limit)
    end = time.time()
    print("par time is" ,end - start)
    print("par res : ",res.result)
    print("par err : ",res.abserr)

    start = time.time()
    res = quad.qag(f,a,b, limit = limit)
    end = time.time()
    print("serial time is", end - start)
    print("serial res : ",res.result)
    print("serial err : ",res.abserr)
