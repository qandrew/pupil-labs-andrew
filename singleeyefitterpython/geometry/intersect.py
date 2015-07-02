"""
	Andrew Xia playing around with porting c++ code to python
	I want to port intersect.h into this python script here
	June 25 2015

"""

#THIS FILE IS INCOMPLETE

#ACTUALLY JUST IGNORE THIS FILE. USE SYMPY


import numpy as np
import scipy

class Line2D:
	def __init__(self, origin = [0,0], direction = None):
		self.origin = np.asarray(origin)
		self.direction = [direction[0]/np.linalg.norm(direction),direction[1]/np.linalg.norm(direction)]
		self.direction = np.asarray(self.direction)
		#to be save, convert all to np arrays

	def __str__(self):
		return "Line { origin:" + str(self.origin) + " direction: " + str(self.direction) + " }"

	# other functions from the eigen class that exist, but may not be used
	def distance(self,point):
		# the distance of a point p to its projection onto the line
		pass
	def intersection_hyperplane(self,hyperplane):
		# the parameter value of intersection between this and given hyperplane
		pass
	def intersection_point(self, hyperplane):
		# returns parameter value of intersection between this and given hyperplane
		pass
	def projection(self,point):
		# returns projection of a point onto the line
		pass
	def pointAt(self,x):
		# point at x along this line
		pass

class Line3D:
	def __init__(self, origin = [0,0,0], direction = None):
		self.origin = np.asarray(origin)
		self.direction = [direction[0]/np.linalg.norm(direction),
			direction[1]/np.linalg.norm(direction),
			direction[2]/np.linalg.norm(direction)]
		self.direction = np.asarray(origin)


	def __str__(self):
		return "Line { origin:" + str(self.origin) + " direction: " + str(self.direction) + " }"

def intersect_2D_lines(line1,line2):
	#finds intersection of 2 lines in 2D. the original intersect() function
	x1 = line1.origin[0]
	y1 = line1.origin[1]
	x2 = line1.origin[0] + line1.direction[0]
	y2 = line1.origin[1] + line1.direction[1]
	x3 = line2.origin[0]
	y3 = line2.origin[1]
	x4 = line2.origin[0] + line2.direction[0]
	y4 = line2.origin[1] + line2.direction[1]

	denom = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
	if (abs(denom) <= 10e-15 ):
		# rounding errors by python since it isn't perfect.
		# though this is sketchy math :P
		denom = 0

	if (denom == 0): #edge case
		#they have the same slope
		if (line1.direction[0] == 0):
			#vertical line, give it some big value
			slope = None
		else:
			slope = line1.direction[1]/line1.direction[0]
		if (x3 == x1):
			x1 = x2 #switch vars
			x3 = x4
			if (x3 == x4):
				if (y3 == y1):
					print "Inputs are same lines, here is one of many points of intersection"
					return np.matrix('%s;%s' % (x1,y1))
				else:
					print "Parallel Lines, no intersect"
					return
		if ((y3-y1)/(x3-x1) == slope):
			#is the same line
			print "Inputs are same lines, here is one of many points of intersection"
			return np.matrix('%s;%s' % (x1,y1))
		else:
			#not the same line
			print "Parallel Lines, no intersect"
			return
	else:
		#there exists an intersection
		px = ((x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*x4 - y3*x4))/denom
		py = ((x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*x4 - y3*x4))/denom
		return np.matrix('%s;%s' % (px,py))

def nearest_intersect_3D_lines(line1, line2):
	# just use the generic nearest intersection function
	return nearest_intersect([line1, line2])

def nearest_intersect_3D(lines):
	#finds the learest intersection of many lines (which may not be a real intersection)
	#the original nearest_intersect(const Range& lines) function
	A = np.zeros((3,3))
	b = np.zeros((3,1))
	Ivv = [] #vector of matrices
	for line in lines:
		vi = np.asmatrix(line.direction)
		vi = np.reshape(vi,(3,1))
		pi = np.asmatrix(line.origin)
		pi = np.reshape(pi,(3,1))

		Ivivi = np.identity(3) - vi*vi.T
		Ivv.append(Ivivi)

		A += Ivivi
		b += Ivivi *pi

	# x = A.partialPivLu.solve(b)
	#not sure if partialPivLu actually does anything...

	return np.linalg.solve(A,b)

def nearest_intersect_2D(lines):
	#finds the learest intersection of many lines (which may not be a real intersection)
	#the original nearest_intersect(const Range& lines) function
	A = np.zeros((2,2))
	b = np.zeros((2,1))
	Ivv = [] #vector of matrices
	for line in lines:
		vi = np.asmatrix(line.direction)
		vi = np.reshape(vi,(2,1))
		pi = np.asmatrix(line.origin)
		pi = np.reshape(pi,(2,1))

		Ivivi = np.identity(2) - vi*vi.T
		Ivv.append(Ivivi)

		A += Ivivi
		b += Ivivi *pi

	#correct to here

	# x = A.partialPivLu.solve(b) #WHAT?
	#not sure if partialPivLu actually does anything...

	return np.linalg.solve(A,b)


def sphere_intersect(line,sphere):
	#intersection between a line and a sphere, originally called intersect(line,sphere)
	"""I have not tested if this function works correctly"""

	v = line.direction
	p = line.origin #put p at origin
	c = sphere.centre - p
	r = sphere.radius

	# from wikipedia :)
	vcvc_cc_rr = np.sq(np.dot(v,c) - np.dot(c,c) + np.sq(r))
	if (vcvc_cc_rr < 0):
		print "NO INTERSECTION between line and sphere"
		return
	s1 = np.dot(v,c) - np.sqrt(vcvc_cc_rr)
	s1 = np.dot(v,c) + np.sqrt(vcvc_cc_rr)

	p1 = p + s1*v
	p2 = p + s2*v
	return (p1,p2)




################################################
if __name__ == '__main__':

	#testing stuff
	huding = Line2D([5,7],[10,10])
	huding2 = Line2D([3,5],[-1,-1])
	print huding
	print intersect_2D_lines(huding, huding2)

	#testing nearest_intersect_3D
	# huding = Line3D([0.835233,3.67143,20], [-0.303085,-0.54173,-0.784008])
	# huding2 = Line3D([0, 0, 0], [-0.0843631, 0.00459802, 0.996424])
	# print nearest_intersect_3D([huding, huding2])

	# huding = Line3D([ 0.,  0, 0], [ -0.427425,-0.293268,-0.855162])
	# huding2 = Line3D([0, 0, 0], [ -0.150507,0.109365,0.982541])
	# print nearest_intersect_3D([huding, huding2])

	# huding = Line3D([ 0.837533,  3.67339, 20], [ -0.427425,-0.293268,-0.855162])
	# huding2 = Line3D([0, 0, 0], [ -0.150507,0.109365,0.982541])
	# print nearest_intersect_3D([huding, huding2])

	# huding = Line3D([0.835451, 3.67313, 20],[0.0687859, -0.42695, -0.901655])
	# huding2 = Line3D([0, 0, 0],[0.0852856, 0.0771611, 0.993364])
	# print nearest_intersect_3D([huding, huding2])

	print "done"