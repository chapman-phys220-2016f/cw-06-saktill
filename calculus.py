#! bin/usr/env python

import numpy as np

def make_diff_matrix(h,N):
    """This function takes in a value h which is the distance between 2 sample points and the amount of points it would sample
    N. It will then use these parameters to create an NxN matrix that is used to creat N central difference
    equations to get the derivative of a function."""
    D = np.zeros((N+1,N+1))
    D[0][0] = -1.0/h
    D[0][1] = 1.0/h
    D[N][N-1] = -1.0/h
    D[N][N] = 1.0/h
    for i in range(1,N):
        for j in range(1,N+1):
			if j == i - 1:
				D[i][j] = -1.0/(2.0*h)
			elif j == i+1:
				D[i][j] = 1/(2.0*h)
	return D

def vec_central_first_diff(f,a,b,N):
	"""Takes a function and an array of equidistant inputs (all discrete)
	and returns the derivative. Achieves this using a derivative matrix.
	Takes args: function f, N sample intervals, [a,b]"""
	x_i = np.linspace(a,b,N+1)
	f_x = f(x_i)
	l = vec_central_first_diff1(x_i, f_x)
	return l

def vec_central_first_diff1(x_i, f_x):
    """This function takes in a list of points x_i and the function that is produced by a function acting on x_i.
    Then it takes the difference matrix defined above and applies it to the function to get the derivative points of the
    function."""
    N = len(x_i) - 1
    h = x_i[1]-x_i[0]
    D = make_diff_matrix(h,N)
    d1 = np.dot(D,f_x)
    return x_i, f_x, d1

def vec_central_second_diff1(x_i,f_x):
    """This function does almost the same thing as vec_central_first_diff1, except this is to get the second derivative.
    So the only difference here is that we multiply the difference matrix to itself and then apply it to the function."""
    N = len(x_i) - 1
    h = x_i[1] - x_i[0]
    D = np.dot(make_diff_matrix(h,N), make_diff_matrix(h,N))
    d1 = np.dot(D,f_x)
    return x_i, f_x, d1

def vec_central_second_diff(f, a, b, N):
    """This is the application function for the second derivative function. Takes in a mathematical function an interval
    [a,b] and number of sample points N, and returns the second derivative."""
    x_i = np.linspace(a,b,N+1)
    f_x = f(x_i)
    l = vec_central_second_diff1(x_i, f_x)
    return l

def vec_central_diff(f,a,b,N,order=1):
    """Application function for vcd. Takes in function f, interval [a,b] and sample points N and the order of the derivatives
    (default to 1)."""
    x_i = np.linspace(a,b,N+1)
    f_x = f(x_i)
    l = vcd(x_i,f_x,order)
    return l


def vcd(x_i,f_x,order):
    """Takes in list of points x_i, function points f_x and a derivative order. Will iterate 'order' amount of times,
    multiplying the difference matrix to itself every time and then applying D^(order) to the function points f_x."""
    N = len(x_i) - 1
    h = x_i[1] - x_i[0]
    D = make_diff_matrix(h,N)
    if order >= 2:
		for i in range(1,order):
			D = np.dot(make_diff_matrix(h,N),D)
    d1 = np.dot(D, f_x)
    return x_i, f_x, d1


