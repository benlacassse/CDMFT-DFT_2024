import numpy as np

layers = 4
liste1 = np.arange(layers)[::-1]
liste2 = np.arange(1, layers)
liste = list(np.hstack((liste1, liste2)))

for l in range(layers):
    for k in range(layers):
        print(liste[layers-1+k-l])