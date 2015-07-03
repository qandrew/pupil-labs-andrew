"""
	Andrew Xia playing around with porting c++ code to python
	This file was originally located in geometry/solve.py, now in the parent folder just as solve.py
	HUDING
	Created July 2 2015

	The inputs to the solve equation are tuples a,b,c

"""

import numpy as np
import scipy
import logging
logging.info('Starting logger for...') 
logger = logging.getLogger(__name__)

def solve(a):
	if (a == 0):
		return 0
	else:
		logger.error("No Solution")
		return

def solve_two(a,b):
	if (a == 0):
		return solve(b)
	return -b/a

def solve_three(a,b,c):
	if (a == 0):
		root = solve_two(b,c)
		return [root,root]

	# http://www.it.uom.gr/teaching/linearalgebra/NumericalRecipiesInC/c5-6.pdf
	# Pg 184
	det = np.square(b) - 4*a*c
	if (det < 0):
		logger.error("No Solution")
		return
	sqrtdet = np.sqrt(det)
	if b >= 0:
		q = -0.5*(b + np.sqrt(det))
	else:
		q = -0.5*(b - np.sqrt(det))
	return [q/a, q/c]

def solve_four(a,b,c,d):

	if (a == 0):
		roots = solve_three(b,c,d)
		return [roots[0],roots[1],roots[1]]

	# http://www.it.uom.gr/teaching/linearalgebra/NumericalRecipiesInC/c5-6.pdf
	# http://web.archive.org/web/20120321013251/http://linus.it.uts.edu.au/~don/pubs/solving.html
	p = b/a
	q = c/a
	r = d/a
	u = q - np.square(p)/3
	v = r - p*q/3 + 2*p*p*p/27
	j = 4*u*u*u/27 * v*v

	if (b == 0 and c == 0):
		return [np.cbrt(-d),np.cbrt(-d),np.cbrt(-d)]
	elif (abs(p) > 10e100):
		return [-p,-p,-p]
	elif (abs(q) > 10e100):
		return [-np.cbrt(v),-np.cbrt(v),-np.cbrt(v)]
	elif (abs(u) > 10e100): #some big number
		return [np.cbrt(4)*u/3,np.cbrt(4)*u/3,np.cbrt(4)*u/3]

	if (j > 0):
		#one real root
		w = sqrt(j)
		if (v > 0):
			y = (u / 3)*np.cbrt(2 / (w + v)) - np.cbrt((w + v) / 2) - p / 3;
		else:
			y = np.cbrt((w - v) / 2) - (u / 3)*np.cbrt(2 / (w - v)) - p / 3;
		return 
	else:
		s = np.sqrt(-u/3)
		t = -v/ (2*s*s*s)
		k = np.arccos(t)/3

		y1 = 2 * s*np.cos(k) - p / 3;
		y2 = s*(-np.cos(k) + np.sqrt(3.)*np.sin(k)) - p / 3;
		y3 = s*(-np.cos(k) - np.sqrt(3.)*np.sin(k)) - p / 3;
		return [float(y1[0]),float(y2[0]),float(y3[0])]

if __name__ == '__main__':

	#testing
	solve(5)
