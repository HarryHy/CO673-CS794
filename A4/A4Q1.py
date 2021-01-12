#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
import matplotlib.pyplot as plt


# In[2]:


#OR
X = [[1,1],[1,-1],[-1,1],[-1,-1]]
X = np.array(X)

Y = np.array([1,1,1,-1])


# In[3]:


#XOR
X = [[1,1],[1,-1],[-1,1],[-1,-1]]
X = np.array(X)

Y = np.array([-1,1,1,-1])


# In[4]:


a = np.array([[0,0,0],[0,0,0],[0,0,0],[0,0,0]])
w0 = np.zeros(3)


# In[5]:


def dykstra(X,y,w,a,maxiter=1000):
    n,d = X.shape
    #X = np.append(X,np.ones((n,1)),axis = 1)
    for t in range(maxiter):
        i = t % n
        z = w + a[i]
        if y[i]==1:
            w_t = z-(X[i]/ np.linalg.norm(X[i])**2)*(np.inner(X[i],z)) - (X[i]/ np.linalg.norm(X[i])**2)* (max(1,np.inner(X[i],z)) )

        else:
            w_t = z-(X[i]/ np.linalg.norm(X[i])**2 )*(np.inner(X[i],z)) - (X[i]/ np.linalg.norm(X[i])**2)*(min(-1,np.inner(X[i],z)) )
        a[i] = a[i] + w - w_t
        w = w_t
    return w


# In[6]:


def dykstra(X,y,w,a,maxiter = 1000):
    [n,d] = X.shape
    #X = np.append(X,np.ones((n,1)),axis=1)
    w_norm = []
    for t in range(maxiter):
        i = t % n
        z = w + a[i]
        if y[i] == 1:
            w_t = z - (X[i]/np.linalg.norm(X[i])**2)*(np.inner(X[i],z)-max(1,np.inner(X[i],z)))
        else:
            w_t = z - (X[i]/np.linalg.norm(X[i])**2)*(np.inner(X[i],z)-min(-1,np.inner(X[i],z)))
        a[i] = a[i] + w - w_t
        w = w_t
        w_norm.append(np.linalg.norm(w))
    return w, w_norm


# In[7]:


def altproj(X,y,w,eta,maxiter=1000):
    w_norm = []
    n,d = X.shape
    X = np.append(X,np.ones((n,1)),axis = 1)
    for t in range(maxiter):
        i = t % n
        if y[i] == 1:
            w_t = w - (X[i]/np.linalg.norm(X[i])**2)*(np.inner(X[i],w)-max(1,np.inner(X[i],w)))
            w = (1-eta)*w + eta*(w_t)
        else:
            w_t = w - (X[i]/np.linalg.norm(X[i])**2)*(np.inner(X[i],w)-min(-1,np.inner(X[i],w)))
            w = (1-eta)*w + eta*(w_t)
        w_norm.append(np.linalg.norm(w))
    return w, w_norm


if __name__ == "__main__":
	#1.4 OR
	X = [[1,1,1],[1,-1,1],[-1,1,1],[-1,-1,1]]
	X = np.array(X)

	Y = np.array([1,1,1,-1])
	w_norm = []
	w, w_norm = dykstra(X,Y,w0,a,maxiter = 1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[0]
	plt.title("1.4 OR")
	plt.plot(x,y)
	plt.show()
	print(w)


	# In[9]:


	#1.5 XOR
	# This Algorithm doesn't blow up
	X = [[1,1,1],[1,-1,1],[-1,1,1],[-1,-1,1]]
	X = np.array(X)

	Y = np.array([-1,1,1,-1])

	w_norm = []
	w, w_norm = dykstra(X,Y,w0,a,maxiter = 1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[1]
	plt.title("1.5 XOR")
	plt.plot(x,y)
	plt.show()
	print(w)
	plt.title("1.5 XOR norm")
	plt.plot(range(1000), w_norm)
	plt.show()


	# In[10]:


	#1.6 OR
	X = [[1,1],[1,-1],[-1,1],[-1,-1]]
	X = np.array(X)

	Y = np.array([1,1,1,-1])

	wa = np.array([0,0,0])
	wb = np.array([1/3,1/3,1/3])
	wc = np.array([1,0,0])
	w_norm = []

	w, w_norm = altproj(X,Y,wa,eta = 1,maxiter=1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[1]
	plt.title("1.6 OR, eta = 1 plot 1")
	plt.plot(x,y)
	plt.show()
	print(w)

	w, w_norm = altproj(X,Y,wb,eta = 1,maxiter=1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[1]
	plt.title("1.6 OR, eta = 1 plot 2")
	plt.plot(x,y)
	plt.show()
	print(w)

	w, w_norm = altproj(X,Y,wc,eta = 1,maxiter=1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[1]
	plt.title("1.6 OR, eta = 1 plot 3")
	plt.plot(x,y)
	plt.show()
	print(w)


	w, w_norm = altproj(X,Y,wa,eta = 2,maxiter=1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[1]
	plt.title("1.6 OR, eta = 2 plot 1")
	plt.plot(x,y)
	plt.show()


	w, w_norm = altproj(X,Y,wb,eta = 2,maxiter=1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[1]
	plt.title("1.6 OR, eta = 2 plot 2")
	plt.plot(x,y)
	plt.show()

	w, w_norm = altproj(X,Y,wc,eta = 2,maxiter=1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[1]
	plt.title("1.6 OR, eta = 2 plot 3")
	plt.plot(x,y)
	plt.show()


	# In[11]:


	#1.7 XOR
	# In Alternating projection eta = 1, the algorithm doesn't blow up
	X = [[1,1],[1,-1],[-1,1],[-1,-1]]
	X = np.array(X)

	Y = np.array([-1,1,1,-1])

	wa = np.array([0.0001,-0.0001,0.0001])
	wb = np.array([0.001,0.001,-0.001])
	wc = np.array([0.01,0.01,0.01])

	w_norm = []
	w, w_norm = altproj(X,Y,wa,eta = 1,maxiter=1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[1]
	plt.title("1.7 XOR, eta = 1 plot 1")
	plt.plot(x,y)
	plt.show()
	print(w)
	plt.title("1.7 XOR, eta = 1 plot 1 w_norm")
	plt.plot(range(1000), w_norm)
	plt.show()

	X = [[1,1],[1,-1],[-1,1],[-1,-1]]
	X = np.array(X)

	Y = np.array([-1,1,1,-1])
	w_norm = []
	w, w_norm = altproj(X,Y,wb,eta = 1,maxiter=1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[1]
	plt.title("1.7 XOR, eta = 1 plot 2")
	plt.plot(x,y)
	plt.show()
	print(w)
	plt.title("1.7 XOR, eta = 1 plot 2 w_norm")
	plt.plot(range(1000), w_norm)
	plt.show()


	X = [[1,1],[1,-1],[-1,1],[-1,-1]]
	X = np.array(X)

	Y = np.array([-1,1,1,-1])
	w_norm = []
	w, w_norm = altproj(X,Y,wc,eta = 1,maxiter=1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[1]
	plt.title("1.7 XOR, eta = 1 plot 3")
	plt.plot(x,y)
	plt.show()
	print(w)
	plt.title("1.7 XOR, eta = 1 plot 3 w_norm")
	plt.plot(range(1000), w_norm)
	plt.show()


	# In[12]:


	#1.7 XOR
	# In Alternating projection eta = 2, the algorithm blows up
	X = [[1,1],[1,-1],[-1,1],[-1,-1]]
	X = np.array(X)

	Y = np.array([-1,1,1,-1])

	wa = np.array([0.0001,-0.0001,0.0001])
	wb = np.array([0.001,0.001,-0.001])
	wc = np.array([0.01,0.01,0.01])

	w_norm = []
	w, w_norm = altproj(X,Y,wa,eta = 2,maxiter=1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[1]
	plt.title("1.7 XOR, eta = 2 plot 1 ")
	plt.plot(x,y)
	plt.show()
	print(w)
	plt.title("1.7 XOR, eta = 2 plot 1 w_norm")
	plt.plot(range(1000), w_norm)
	plt.show()

	X = [[1,1],[1,-1],[-1,1],[-1,-1]]
	X = np.array(X)

	Y = np.array([-1,1,1,-1])
	w_norm = []
	w, w_norm = altproj(X,Y,wb,eta = 2,maxiter=1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[1]
	plt.title("1.7 XOR, eta = 2 plot 2 ")
	plt.plot(x,y)
	plt.show()
	print(w)
	plt.title("1.7 XOR, eta = 2 plot 2 w_norm")
	plt.plot(range(1000), w_norm)
	plt.show()


	X = [[1,1],[1,-1],[-1,1],[-1,-1]]
	X = np.array(X)

	Y = np.array([-1,1,1,-1])
	w_norm = []
	w, w_norm = altproj(X,Y,wc,eta = 2,maxiter=1000)
	x_point = np.array([1,1,-1,-1])
	y_point = np.array([1,-1,1,-1])
	label = Y
	plt.scatter(x_point, y_point, c = label)
	x = np.linspace(-3,3,120)
	y = -w[0]/w[1]*x -w[2]/w[1]
	plt.title("1.7 XOR, eta = 2 plot 3 ")
	plt.plot(x,y)
	plt.show()
	print(w)
	plt.title("1.7 XOR, eta = 2 plot 3 w_norm")
	plt.plot(range(1000), w_norm)
	plt.show()




