# import uniform distribution
from scipy.stats import uniform
import math
import random

r = 1
unit_delta = 0.0001


tipe = "uniform"
n = 2500000
start = -1
width = 2
x = uniform.rvs(size=n, loc = start, scale=width)
delta = uniform.rvs(size=n, loc = -unit_delta, scale=2*unit_delta)

y = list()
for z in x:
  signum = random.choice((-1, 1))
  Y = math.sqrt(pow(r,2) - pow(z,2)) * int(signum)
  y.append(Y)

x_delta = list()
y_delta = list()
for i in range(len(x)):
    x_delta.append(x[i]+delta[i])
    y_delta.append(y[i]+delta[i])


namefile= "../Data/data.txt"
f = open(namefile, "a")

for i in range(len(x_delta)):
  if(i != len(x_delta)-1):
    f.write("{"+str(x_delta[i])+", "+str(y_delta[i])+"}"+"\n")
  else:
    f.write("{"+str(x_delta[i])+", "+str(y_delta[i])+"}")
f.close()
