"""
	Andrew Xia playing around with porting c++ code to python
	I want to port Conic.h into this python script here
	June 25 2015

"""

import numpy as np

class Conic:

	def __init__(self, ellipse=None):
		#a Conic is defined by 6 scalar parameters A,B,C,D,E,F

		if ellipse != None:
			#extracting information from ellipse
			ax = np.cos(ellipse.angle)
			ay = np.sin(ellipse.angle)

			a2 = np.square(ellipse.major_radius)
			b2 = np.square(ellipse.minor_radius)

			#scalars
			self.A = (ax*ax)/a2 + (ay*ay)/b2
			self.B = 2*(ax*ay)/a2 - 2*(ax*ay)/b2
			self.C = (ay*ay)/a2 +(ax*ax)/b2
			self.D = (-2*ax*ay*ellipse.centre[1] - 2*ax*ax*ellipse.centre[0])/a2 + (2*ax*ay*ellipse.centre[1] - 2*ay*ay*ellipse.centre[0])/b2
			self.E = (-2*ax*ay*ellipse.centre[0] - 2*ay*ay*ellipse.centre[1])/a2 + (2*ax*ay*ellipse.centre[0] - 2*ax*ax*ellipse.centre[1])/b2
			self.F = (2*ax*ay*ellipse.centre[0]*ellipse.centre[1]+ax*ax*ellipse.centre[0]*ellipse.centre[0]+ay*ay*ellipse.centre[1]*ellipse.centre[1])/a2+ (-2*ax*ay*ellipse.centre[0]*ellipse.centre[1]+ ay*ay*ellipse.centre[0]*ellipse.centre[0]+ax*ax*ellipse.centre[1]*ellipse.centre[1])/b2-1

		else:
			self.A = float()
			self.B = float()
			self.C = float()
			self.D = float()
			self.E = float()
			self.F = float()

	def operator(self,x,y):
		#this function returns the conic based on coordinates x and y
		return A*x*x + B*x*y + C*y*y + D*x + E*y + F

	def __str__(self):
		return "Conic { " + str(self.A) + "x^2 + " + str(self.B) + 'xy + ' + str(self.C) + 'y^2 + ' + str(self.D) + 'x + ' + str(self.E) + 'y + ' + str(self.F) + " = 0 }"

	def transformed(self,a,t):
		#this function returns the transformed conic

		#a is a Eigen::MatrixBase, the affine transform matrix
		#t is a Eigen::MatrixBase, the translation matrix
		""" We map x,y to a new space using [x y] -> affine*[x y] + translation
		    
		     Using a for affine and t for translation:
		         x -> a_00*x + a01*y + t0
		         y -> a_10*x + a11*y + t1
		    
		     So
		         Ax^2 + Bxy + Cy^2 + Dx + Ey + F
		     becomes
		           A(a_00*x + a01*y + t0)(a_00*x + a01*y + t0)
		         + B(a_00*x + a01*y + t0)(a_10*x + a11*y + t1)
		         + C(a_10*x + a11*y + t1)(a_10*x + a11*y + t1)
		         + D(a_00*x + a01*y + t0)
		         + E(a_10*x + a11*y + t1)
		         + F
		    
		     Collecting terms gives:
		"""

		return (
                A*np.square(a(0, 0)) + B*a(0, 0)*a(1, 0) + C*np.square(a(1, 0)),
                2 * A*a(0, 0)*a(0, 1) + B*a(0, 0)*a(1, 1) + B*a(0, 1)*a(1, 0) + 2 * C*a(1, 0)*a(1, 1),
                A*np.square(a(0, 1)) + B*a(0, 1)*a(1, 1) + C*np.square(a(1, 1)),
                2 * A*a(0, 0)*t(0) + B*a(0, 0)*t(1) + B*a(1, 0)*t(0) + 2 * C*a(1, 0)*t(1) + D*a(0, 0) + E*a(1, 0),
                2 * A*a(0, 1)*t(0) + B*a(0, 1)*t(1) + B*a(1, 1)*t(0) + 2 * C*a(1, 1)*t(1) + D*a(0, 1) + E*a(1, 1),
                A*np.square(t(0)) + B*t(0)*t(1) + C*np.square(t(1)) + D*t(0) + E*t(1) + F
                )
