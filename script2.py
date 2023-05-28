import quad
import math
import time
from scipy.integrate import quad_vec
import numpy as np
from scipy.integrate import quad as q

var = range(1)
for a in var:
    f = lambda x : (math.cos(x),math.sin(x))
    epsabs = 1.0
    epsrel = 0.0
    a = 0.0
    b = 10000.0

    f = lambda x : np.array([math.cos(1/(x*x)),math.sin(1/(x*x))])

    start = time.time()
    y , err = quad_vec(f,a,b,epsabs=epsabs,epsrel=epsrel,limit=1000000)
    end = time.time()
    print(end - start)
    print(y)
    print(err)

    #f = lambda x : math.cos(x)
    #start = time.time()
    #y , err = q(f,a,b,epsabs=epsabs,epsrel=epsrel,limit=1000000)
    #f = lambda x : math.sin(x)
    #y , err = q(f,a,b,epsabs=epsabs,epsrel=epsrel,limit=1000000)
    #end = time.time()
    #print(end - start)
    #print(y)
    #print(err)






    f = lambda x : [math.cos(x)]
    res = quad.qags2(f,a,b,epsabs,epsrel,6,10000,1, more_info = True)
    print(res.more_info)

    f = lambda x : [math.cos(x),math.cos(x),math.cos(x),math.cos(x)]
    res = quad.qags2(f,a,b,epsabs,epsrel,6,10000,4, more_info = True)
    print(res.more_info)