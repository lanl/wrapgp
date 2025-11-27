import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

seed = int(sys.argv[1])
method = int(sys.argv[2])
n = int(sys.argv[3])

# load appropriate python script
if method == 1: # geomstats
  sys.path.append("../../python")
  from wrapgp import fit_wrapgp
  
np.random.seed(10)

# load training and testing data
train_name = "data/train_n" + str(n) + "_seed" + str(seed) + ".csv"
test_name = "data/test_seed" + str(seed) + ".csv"
train = pd.read_csv(train_name)
test = pd.read_csv(test_name)
X = np.array(train["X"]).reshape(-1, 1)
Y = np.column_stack((np.cos(train["Y"]), np.sin(train["Y"])))
XX = np.array(test["X"]).reshape(-1, 1)
YY = test["Y"]

# fit specified model
if method == 1:
  mu, sd, draws = fit_wrapgp(X, Y, XX)
  name = "mallasto"

#import matplotlib.pyplot as plt
#plt.plot(X, train["Y"], "o", XX, mu, "o", XX, mu - 1.96*sd, "o", XX, mu + 1.96*sd, "o")
#plt.show()

# save predictions/draws
df = pd.DataFrame({'X': XX[:,0], 'Y': YY, 'm': mu, 's2': sd**2})
pred_name = "pred/" + str(name) + "_n" + str(n) + "_seed" + str(seed) + ".csv"
df.to_csv(pred_name, index = False)
draws_name = "draws/" + str(name) + "_n" + str(n) + "_seed" + str(seed) + ".csv"
pd.DataFrame(draws).to_csv(draws_name, index = False)
  
