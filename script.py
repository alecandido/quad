import quad
import math
import time
import cProfile
import numpy as np
from scipy.integrate import quad_vec
from scipy.integrate import quad as q

var = range(2)
for a in var:
    f = lambda x : (math.cos(x),math.sin(x))
    epsabs = 1.0e-3
    epsrel = 0.0
    a = 0.0
    b = 10000.0
    #start = time.time()
    #res = quad.qag(f,a,b,epsabs,epsrel,1,2)
    #end = time.time()
    #print(end - start)
    #print(res)

    start = time.time()
    res = quad.qag(f,a,b,epsabs,epsrel,6,2)
    end = time.time()
    print(end - start)
    print(res)

    start = time.time()
    res = quad.qags(f,a,b,epsabs,epsrel,6,10000,2)
    end = time.time()
    print(end - start)
    print(res)

    #start = time.time()
    #res = quad.qagv2(f,a,b,epsabs,epsrel,6,2)
    #end = time.time()
    #print(end - start)
    #print(res)

    start = time.time()
    res = quad.qagv3(f,a,b,epsabs,epsrel,6,2)
    end = time.time()
    print(end - start)
    print(res)

    #f = lambda x : np.array([math.cos(1/(x*x)),math.sin(1/(x*x))])
#
    #start = time.time()
    #y , err = quad_vec(f,a,b,epsabs=epsabs,epsrel=epsrel,limit=1000000)
    #end = time.time()
    #print(end - start)
    #print(y)
    #print(err)

    f = lambda x : math.cos(x)
    start = time.time()
    y , err = q(f,a,b,epsabs=epsabs,epsrel=epsrel,limit=1000000)
    f = lambda x : math.sin(x)
    y , err = q(f,a,b,epsabs=epsabs,epsrel=epsrel,limit=1000000)
    end = time.time()
    print(end - start)
    print(y)
    print(err)


#cProfile.run('quad.qag(f,a,b,epsabs,epsrel,6,8)','restats')
#cProfile.run("quad_vec(f,a,b,epsabs=epsabs,epsrel=epsrel,limit=1000000)",'restats2')
