from  scipy.linalg import solve_banded

import numpy as np
import pylab as plt

A = 201e-6;
E = 1.4e11;
slen = 0.1;
Kbond = 8;
Sbond = 2.8e5;
n = 5;
Yield = 180e3;
Rupture = 0.2;
ds = 1e7;
dy = 62.8e3;
dr = 0.41;
theta = 40.0;
disp = 4e-3;

n = 2*n
Kbond = 10**Kbond

c1 = E*A/slen
c2 = Kbond*slen
c3 = Sbond*slen
c2i = np.ones(n)* c2

db0 = 0.0;
db1 = disp;
ab      = np.zeros((3,n))
b       = np.zeros(n)


for i in range(n):
    if i>int(n/2-1):
        b[i] = -db1*c2i[i]
    else:
        b[i] = -db0*c2i[i]
    if i==0:
        ab[0][0] = 0
        ab[1][0] = -c1-c2i[i]
        ab[2][0] = c1
    elif i==(n-1):
        ab[0][i] = c1
        ab[1][i] = -c1-c2i[i]
        ab[2][i] = 0
    else:
        ab[0][i] = c1
        ab[1][i] = -2*c1-c2i[i]
        ab[2][i] = c1
ndis = solve_banded((1,1),ab,b)

x = np.linspace(0, n*slen, n)

print ndis
plt.plot(x,ndis)
plt.show()
