# author:高金磊
# datetime:2021/10/29 10:51
#Activation function
import numpy as np
def sigmoid(x):
    return 1 / (1 + np.power(np.e,-x))

def powersum(x):
    return x+np.power(x,3)+np.power(x,5)

def linear(x):
    return x

def powersigmoid(X):
    xi=4
    p=3
    output=(1 + np.exp(-xi)) / (1 - np.exp(-xi)) *np.divide(1 - np.power(np.e,-xi * X) , 1 + np.power(np.e,-xi * X))
    for i in range(len(X)):
        if X[i]<=-1 or X[i]>=1:
                output[i]=np.power(X[i],p)
    return output
def Bound(X):
    for i in range(len(X)):
        if abs(X[i])>1:
            X[i]/=abs(X[i])
    return X
# if nargin==1
#     xi=4;
#     p=3;
# elseif nargin==2
#     p=3;
# end
# %论文2.21
# output=(1+exp(-xi))/(1-exp(-xi))*(1-exp(-xi*X))./(1+exp(-xi*X));
# i=find(abs(X)>=1);
# output(i)=X(i).^p;

if __name__ == '__main__':
    print(Bound([1,-1,1.1,-1.1,0.5]))

