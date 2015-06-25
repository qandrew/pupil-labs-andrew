"""
	Andrew Xia playing around with porting c++ code to python
	I want to port intersect.h into this python script here
	June 25 2015

"""

#THIS FILE IS INCOMPLETE

#ACTUALLY JUST IGNORE THIS FILE. USE SYMPY


import numpy as np
import sympy

def intersect_2D_lines(line1,line2):
	#finds intersection of 2 lines in 2D. the original intersect() function
	"""	PASS
		x1 = line1.origin[0]
		y1 = line1.origin[1]
		x2 = (line1.origin + line1.direction)[0]
		y2 = (line1.origin + line1.direction)[1]
		x3 = line2.origin[0]
		y3 = line2.origin[1]
		x4 = (line2.origin + line2.direction)[0]
		y4 = (line2.origin + line2.direction)[1]

		denom = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)

		print x1,x2,x3,x4,denom
		px = ((x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*x4 - y3*x4))/denom
		py = ((x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*x4 - y3*x4))/denom

		return np.matrix('%s;%s' % (px,py))
	"""

def nearest_intersect(lines):
	#finds the learest intersection of many lines (which may not be a real intersection)
	#the original nearest_intersect(const Range& lines) function
	Vector = np.array()
	Matrix = np.matrix()


def nearest_intersect_3D_lines(line1, line2):
	#if 2D, then just use intersect_2D_line code
	if (dimension == 2):
		return intersect_2D_lines(line1,line2)
	#for lines that are not 2D, was the original nearest() function for general-dimension lines
	else:
		lines = [line1,line2]
		return nearest_intersect(lines)

def sphere_intersect(line,sphere):
	#intersection between a line and a sphere, originally called intersect(line,sphere)
	Vector = []

a1 = sympy.Line((0,0,0),(3,3))
a2 = sympy.Line((0,0),(3,10))
c = a1.intersection(sympy.Point(4,4))
if c:
	print "ay"