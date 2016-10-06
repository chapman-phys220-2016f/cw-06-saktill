#! bin/usr/env python

import numpy as np

def make_diff_matrix(h,N):
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
	N = len(x_i) - 1
	h = x_i[1]-x_i[0]
	D = make_diff_matrix(h,N)
	d1 = np.dot(D,f_x)
	return x_i, f_x, d1

def vec_central_second_diff1(x_i,f_x):
	N = len(x_i) - 1
	h = x_i[1] - x_i[0]
	D = np.dot(make_diff_matrix(h,N), make_diff_matrix(h,N))
	d1 = np.dot(D,f_x)
	return x_i, f_x, d1

def vec_central_second_diff(f, a, b, N):
	x_i = np.linspace(a,b,N+1)
	f_x = f(x_i)
	l = vec_central_second_diff1(x_i, f_x)
	return l

def vec_central_diff(f,a,b,N,order=1):
	x_i = np.linspace(a,b,N+1)
	f_x = f(x_i)
	l = vcd(x_i,f_x,order)
	return l


def vcd(x_i,f_x,order):
	N = len(x_i) - 1
	h = x_i[1] - x_i[0]
	D = make_diff_matrix(h,N)
	if order >= 2:
		for i in range(1,order):
			D = np.dot(make_diff_matrix(h,N),D)
	d1 = np.dot(D, f_x)
	return x_i, f_x, d1


