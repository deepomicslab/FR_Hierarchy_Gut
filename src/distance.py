# function to compute distance of two columns
import numpy as np
import copy
import pandas as pd

def Jaccard_distance(x, y):
    union = sum(np.maximum(x, y))
    intersection = sum(np.minimum(x, y))
    # print(intersection, union)
    if union != 0:
        return 1 - intersection/union
    else:
        # print("union size is 0")
        return 1

def corr_distance(x, y):
    numerator = np.dot(x, y)
    denominator = np.linalg.norm(x)*np.linalg.norm(y)
    if denominator == 0:
        return 1
    else:
        return 1 - numerator/denominator

def sorensen_distance(x, y):
    intersection = sum(np.minimum(x, y))
    denominator = sum(x) + sum(y)
    if denominator == 0:
        return 0
    else:
        1 - 2*intersection/denominator

def hamming_distance(x, y):
    # x, y must be binary
    d = sum(abs(x-y))/len(x)
    return d

def mutual_info_distance(x, y):
    hx = entropy(x)
    hy = entropy(y)
    Ixy = hx + hy - cross_entropy(x, y)
    nmi = 2*Ixy/(hx + hy)
    return 1- nmi
    
def entropy(x):
    x_types = list(set(x))
    x_types.sort()
    hx = 0
    unit = 1/len(x)
    entropy_dict = {}
    for xtype in x_types:
        entropy_dict[xtype] = 0
        pxtype = 0
    for xi in x:
        entropy_dict[xi] += unit
    for xtype, pxtype in entropy_dict.items():
        hx -= pxtype*np.log(pxtype)
    return hx

def cross_entropy(x, y):
    unit = 1/len(x)
    xlist = list(set(x))
    xlist.sort()
    ylist = list(set(y))
    ylist.sort()
    hxy = 0
    data = np.zeros(shape=(len(xlist), len(ylist)))
    entropy_df = pd.DataFrame(data, index=xlist, columns=ylist)

    for xi in x:
        for yi in y:
            #print(xi,yi)
            entropy_df.loc[int(xi), int(yi)] += unit

    for xtype in xlist:
        for ytype in ylist:
            cell = entropy_df.loc[int(xtype), int(ytype)]
            hxy -= cell*np.log(cell)
    return hxy
            
def NMI(X):
    n = X.shape[0]
    #print(n)
    result = np.zeros((n,n))
    hx_list = [0]*n
    for i in range(n):
        x = X[i, :]
        hx_list[i] = entropy(x)
        #print(hx_list)
        for j in range(i, 0, -1):
            y = X[j, :]
            hxy = cross_entropy(x, y)
            Ixy = hx_list[i] + hx_list[j] - hxy
            nmi = 2*Ixy/(hx_list[i] + hx_list[j])
            result[j,i] = nmi
    return result
