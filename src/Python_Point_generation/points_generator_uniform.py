# import uniform distribution
from scipy.stats import uniform
# random numbers from uniform distribution



tipe = "uniform"
n = 400
start = 0
width = 1
uniform_x = uniform.rvs(size=n, loc = start, scale=width)
uniform_y = uniform.rvs(size=n, loc = start, scale=width)

namefile= "../Data/data.txt"
f = open(namefile, "a")

for i in range(len(uniform_x)):
  if(i != len(uniform_x)-1):
    f.write("{"+str(uniform_x[i])+", "+str(uniform_y[i])+"}"+"\n")
  else:
    f.write("{"+str(uniform_x[i])+", "+str(uniform_y[i])+"}")
f.close()
