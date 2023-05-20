import python_bindings
import math
import time
import numpy as np
from scipy.integrate import quad_vec

f = lambda x : (x,x*x)
python_bindings.integrate_vecc(f,0.0,1.0)
f = lambda x : (math.cos(1/(x*x)),math.cos(1/(x*x)),math.sin(1/(x*x)),math.sin(1/(x*x)))
epsabs = 0.01
epsrel = 0.0
a = 0.0
b = 1.0
start = time.time()
y , err = python_bindings.integrate_vecx(f,a,b,epsabs,epsrel)
end = time.time()
print(end - start)
print(y)
print(err)

f = lambda x : np.array([math.cos(1/(x*x)),math.cos(1/(x*x)),math.sin(1/(x*x)),math.sin(1/(x*x))])


start = time.time()
y , err = quad_vec(f,a,b,epsabs=epsabs,epsrel=epsrel,limit=1000000)
end = time.time()
print(end - start)
print(y)
print(err)


