from scipy.stats import norm
# generate random numbers from N(0,1)
# data_normal = norm.rvs(size=10000,loc=0,scale=1)

n = 10
mu = 0
sigma = 1
norm_x = norm.rvs(size=n, loc = mu, scale=sigma)
norm_y = norm.rvs(size=n, loc = mu, scale=sigma)

namefile= "../Data/data.txt"
f = open(namefile, "a")

for i in range(len(norm_x)):
  if(i != len(norm_x)-1):
    f.write("{"+str(norm_x[i])+", "+str(norm_y[i])+"}"+"\n")
  else:
    f.write("{"+str(norm_x[i])+", "+str(norm_y[i])+"}")
f.close()
