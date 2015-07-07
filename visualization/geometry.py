"""
	Andrew Xia playing around with porting c++ code to python
	This file contains the Circle.py, Conic,py, Ellipse.py, Conicoic.py, and Sphere.py
	Having all 5 geometry classes in one file allows for better maintainability.
	Created July 2 2015

"""

import numpy as np
import scipy
import logging
logging.info('Starting logger for...') 
logger = logging.getLogger(__name__)

class Circle3D:
	#originally from Circle.py, see older version on github for further details of the original file
	def __init__(self,center=[0,0,0],normal=[0,0,0],radius=0):
		self.center = np.array(center)
		self.normal = np.array(normal)
		self.radius = radius

	def __str__(self):
		return "Circle { center: " + str(self.center) + ", normal: " + str(self.normal) + ", radius: " + str(self.radius) + " } "

	def is_same(self,circle):
		return  (self.center == circle.center and self.normal == circle.normal and self.radius == circle.radius)

class Conic:
	#originally from Circle.py, see older version on github for further details of the original file
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
			self.D = (-2*ax*ay*ellipse.center[1] - 2*ax*ax*ellipse.center[0])/a2 + (2*ax*ay*ellipse.center[1] - 2*ay*ay*ellipse.center[0])/b2
			self.E = (-2*ax*ay*ellipse.center[0] - 2*ay*ay*ellipse.center[1])/a2 + (2*ax*ay*ellipse.center[0] - 2*ax*ax*ellipse.center[1])/b2
			self.F = (2*ax*ay*ellipse.center[0]*ellipse.center[1]+ax*ax*ellipse.center[0]*ellipse.center[0]+ay*ay*ellipse.center[1]*ellipse.center[1])/a2+ (-2*ax*ay*ellipse.center[0]*ellipse.center[1]+ ay*ay*ellipse.center[0]*ellipse.center[0]+ax*ax*ellipse.center[1]*ellipse.center[1])/b2-1

		else:
			self.A = float()
			self.B = float()
			self.C = float()
			self.D = float()
			self.E = float()
			self.F = float()

	def operator(self,x,y):
		#this function returns the conic based on coordinates x and y
		if (x == None or y == None):
			logger.error("Inputs x or y is none")
			return
		return self.A*x*x + self.B*x*y + self.C*y*y + self.D*x + self.E*y + self.F

	def __str__(self):
		return "Conic { " + str(self.A) + "x^2 + " + str(self.B) + 'xy + ' + str(self.C) + 'y^2 + ' + str(self.D) + 'x + ' + str(self.E) + 'y + ' + str(self.F) + " = 0 }"

	def transformed(self,a,t):
		#this function returns the transformed conic

		# HAVENT TESTED THIS FUNCTION TO CONFIRM I WROTE IT WITHOUT SYNTAX ERROR

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

		toreturn = Conic(A*np.square(a(0, 0)) + B*a(0, 0)*a(1, 0) + C*np.square(a(1, 0)),
                2 * A*a(0, 0)*a(0, 1) + B*a(0, 0)*a(1, 1) + B*a(0, 1)*a(1, 0) + 2 * C*a(1, 0)*a(1, 1),
                A*np.square(a(0, 1)) + B*a(0, 1)*a(1, 1) + C*np.square(a(1, 1)),
                2 * A*a(0, 0)*t(0) + B*a(0, 0)*t(1) + B*a(1, 0)*t(0) + 2 * C*a(1, 0)*t(1) + D*a(0, 0) + E*a(1, 0),
                2 * A*a(0, 1)*t(0) + B*a(0, 1)*t(1) + B*a(1, 1)*t(0) + 2 * C*a(1, 1)*t(1) + D*a(0, 1) + E*a(1, 1),
                A*np.square(t(0)) + B*t(0)*t(1) + C*np.square(t(1)) + D*t(0) + E*t(1) + F
                )
		return toreturn

class Conicoid:

	def __init__(self,A=0.0,B=0.0,C=0.0,F=0.0,G=0.0,H=0.0,U=0.0,V=0.0,W=0.0,D=0.0,conic = None, vertex = None):
		self.A = A
		self.B = B
		self.C = C
		self.D = D
		self.F = F
		self.G = G
		self.H = H
		self.U = U
		self.V = V
		self.W = W

		if (conic != None and vertex != None):
			self._initialize_by_conic(conic,vertex)

	def __str__(self):
		#return "Conicoid, A: " + str(self.A) + ", B: " + str(self.B) + ", C: "  + str(self.C) + ", F: " + str(self.F) + ", G: " + str(self.G) + ", H: " + str(self.H) + ", U: " + str(self.U) + ", V: " + str(self.V) + ", W: " + str(self.W) + ", D: " + str(self.D)
		return "Conicoid { " + str(self.A) + "x^2 + " + str(self.B) + "y^2 + "  + str(self.C) + "z^2 + " + str(self.F) + "yz + 2*" + str(self.G) + "zx +2*" + str(self.H) + "xy + " + str(self.U) + "x + 2*" + str(self.V) + "y + 2*" + str(self.W) + "z + " + str(self.D) + " = 0 }"

	def _initialize_by_conic(self,conic,vertex):
		#private method
		alpha = vertex[0]
		beta = vertex[1]
		gamma = vertex[2]

		self.A = np.square(gamma)*conic.A
		self.B = np.square(gamma)*conic.C
		self.C = np.square(alpha)*conic.A + alpha*beta*conic.B + np.square(beta)*conic.C + conic.D*alpha + conic.E*beta + conic.F
		self.F = -gamma * (conic.C * beta + conic.B / 2 * alpha + conic.E / 2)
		self.G = -gamma * (conic.B / 2 * beta + conic.A * alpha + conic.D / 2)
		self.H = np.square(gamma)*conic.B/2
		self.U = np.square(gamma)*conic.D/2
		self.V = np.square(gamma)*conic.E/2
		self.W = -gamma * (conic.E / 2 * beta + conic.D / 2 * alpha + conic.F)
		self.D = np.square(gamma)*conic.F

	def operator(self,x,y,z):
		return self.A*np.square(x) + self.B*np.square(y) + self.C*np.square(z) + 2 * self.F*y*z + 2 * self.G*x*z + 2 * self.H*x*y + 2 * self.U*x + 2 * self.V*y + 2 * self.W*z + self.D

	def intersectZ(self,z=0):
		# finds a conic at a given z intersection
		# Ax^2 + By^2 + Cz^2 + 2Fyz + 2Gzx + 2Hxy + 2Ux + 2Vy + 2Wz + D = 0
		# becomes
		# Ax^2 + Bxy + Cy^2 + Fx + Ey + D = 0
		toreturn = Conic()
		toreturn.A = self.A
		toreturn.B = 2*self.H
		toreturn.C = self.C
		toreturn.D = 2*self.G*z + 2*self.U
		toreturn.E = 2*self.F*z + 2*self.V
		toreturn.F = self.C*np.square(z) + 2*self.W*z + self.D
		return toreturn

	def transformed(a,t):
		"""We map x,y,z to a new space using [x y z] -> affine*[x y z] + translation
			
			Using a for affine and t for translation:
				x -> a_00*x + a01*y + a02*z + t0
				y -> a_10*x + a11*y + a12*z + t1
				z -> a_20*x + a21*y + a22*z + t2
			
			So
				Ax^2 + By^2 + Cz^2 + 2Fyz + 2Gzx + 2Hxy + 2Ux + 2Vy + 2Wz + D
			becomes
				  A(a_00*x + a01*y + a02*z + t0)(a_00*x + a01*y + a02*z + t0)
				+ B(a_10*x + a11*y + a12*z + t1)(a_10*x + a11*y + a12*z + t1)
				+ C(a_20*x + a21*y + a22*z + t2)(a_20*x + a21*y + a22*z + t2)
				+ 2F(a_10*x + a11*y + a12*z + t1)(a_20*x + a21*y + a22*z + t2)
				+ 2G(a_20*x + a21*y + a22*z + t2)(a_00*x + a01*y + a02*z + t0)
				+ 2H(a_00*x + a01*y + a02*z + t0)(a_10*x + a11*y + a12*z + t1)
				+ 2U(a_00*x + a01*y + a02*z + t0)
				+ 2V(a_10*x + a11*y + a12*z + t1)
				+ 2W(a_20*x + a21*y + a22*z + t2)
				+ D
			
			Collecting terms gives:
		"""
		a = self.A*np.square(a(0, 0)) + self.B*np.square(a(1, 0)) + self.C*np.square(a(2, 0)) + 2*self.F*a(1, 0)*a(2, 0) + 2*self.G*a(0, 0)*a(2, 0) + 2*self.H*a(0, 0)*a(1, 0)
		b = self.A*np.square(a(0, 1)) + self.B*np.square(a(1, 1)) + self.C*np.square(a(2, 1)) + 2*self.F*a(1, 1)*a(2, 1) + 2*self.G*a(0, 1)*a(2, 1) + 2*self.H*a(0, 1)*a(1, 1),
		c = self.A*np.square(a(0, 2)) + self.B*np.square(a(1, 2)) + self.C*np.square(a(2, 2)) + 2*self.F*a(1, 2)*a(2, 2) + 2*self.G*a(0, 2)*a(2, 2) + 2*self.H*a(0, 2)*a(1, 2),
		f = self.A*a(0, 1)*a(0, 2) + self.B*a(1, 1)*a(1, 2) + self.C*a(2, 1)*a(2, 2) + self.F*a(1, 1)*a(2, 2) + self.F*a(1, 2)*a(2, 1) + self.G*a(0, 1)*a(2, 2) + self.G*a(0, 2)*a(2, 1) + self.H*a(0, 1)*a(1, 2) + self.H*a(0, 2)*a(1, 1),
		g = self.A*a(0, 2)*a(0, 0) + self.B*a(1, 2)*a(1, 0) + self.C*a(2, 2)*a(2, 0) + self.F*a(1, 2)*a(2, 0) + self.F*a(1, 0)*a(2, 2) + self.G*a(0, 2)*a(2, 0) + self.G*a(0, 0)*a(2, 2) + self.H*a(0, 2)*a(1, 0) + self.H*a(0, 0)*a(1, 2),
		h = self.A*a(0, 0)*a(0, 1) + self.B*a(1, 0)*a(1, 1) + self.C*a(2, 0)*a(2, 1) + self.F*a(1, 0)*a(2, 1) + self.F*a(1, 1)*a(2, 0) + self.G*a(0, 0)*a(2, 1) + self.G*a(0, 1)*a(2, 0) + self.H*a(0, 0)*a(1, 1) + self.H*a(0, 1)*a(1, 0),
		u = self.A*a(0, 0)*t(0) + self.B*a(1, 0)*t(1) + self.C*a(2, 0)*t(2) + self.F*a(1, 0)*t(2) + self.F*a(2, 0)*t(1) + self.G*a(0, 0)*t(2) + self.G*a(2, 0)*t(0) + self.H*a(0, 0)*t(1) + self.H*a(1, 0)*t(0) + self.U*a(0, 0) + self.V*a(1, 0) + self.W*a(2, 0),
		v = self.A*a(0, 1)*t(0) + self.B*a(1, 1)*t(1) + self.C*a(2, 1)*t(2) + self.F*a(1, 1)*t(2) + self.F*a(2, 1)*t(1) + self.G*a(0, 1)*t(2) + self.G*a(2, 1)*t(0) + self.H*a(0, 1)*t(1) + self.H*a(1, 1)*t(0) + self.U*a(0, 1) + self.V*a(1, 1) + self.W*a(2, 1),
		w = self.A*a(0, 2)*t(0) + self.B*a(1, 2)*t(1) + self.C*a(2, 2)*t(2) + self.F*a(1, 2)*t(2) + self.F*a(2, 2)*t(1) + self.G*a(0, 2)*t(2) + self.G*a(2, 2)*t(0) + self.H*a(0, 2)*t(1) + self.H*a(1, 2)*t(0) + self.U*a(0, 2) + self.V*a(1, 2) + self.W*a(2, 2),
		d = self.A*np.square(t(0)) + B*np.square(t(1)) + self.C*np.square(t(2)) + 2*self.F*t(1)*t(2) + 2*self.G*t(0)*t(2) + 2*self.H*t(0)*t(1) + 2*self.U*t(0) + 2*V*t(1) + 2*self.W*t(2) + self.D
		return Conicoid(a,b,c,d,f,g,h,u,v,w,d)

class Ellipse:

	def __init__(self,center=[0,0],major_radius=0.0,minor_radius=0.0,angle=0.0, conic = None):
		self.center = np.array(center)
		self.major_radius = major_radius
		self.minor_radius = minor_radius
		self.angle = angle
		if conic:
			self._initialize_by_conic(conic)

	def _initialize_by_conic(self, conic):
		#to initialize with a conic, create the ellipse first (temp = Ellipse())
		#then call with the conic: (temp.initialize_conic(some_conic)). 
		#this is because python can only have 1 __init__() function.

		self.angle = 0.5*np.arctan2(conic.B,conic,A - conic.C)
		cost = np.cos(angle)
		sint = np.sin(angle)
		cos_squared = np.square(cost)
		sin_squared = np.square(sint)

		Ao = conic.F
		Au = conic.D*cost + conic.E*sint
		Av = -conic.D*sint + conic.E*cost
		Auu= conic.A*cos_squared +conic.C*sin_squared +conic.B*sint*cost
		Avv= conic.A*sin_squared +conic.C*sin_squared -conic.B*sint*cost
		#ROTATED = [Ao Au Av Auu Avv]
		tuCenter = -Au / (2.0*Auu)
		tvCenter = -Av / (2.0*Auu)
		wCenter = Ao - Auu*np.square(tuCenter) - Avv*np.square(tvCenter)

		self.center = [0,0]
		self.center[0] = tuCenter*cost - tvCenter*sint
		self.center[1] = tuCenter*sint - tvCenter*cost

		self.major_radius = np.sqrt(abs(-wCenter/Auu))
		self.minor_radius = np.sqrt(abs(-wCenter/Avv))

		if (self.major_radius < self.minor_radius):
			major_radius,minor_radius = minor_radius,major_radius
			self.angle = self.angle + scipy.pi/2

		if (self.angle > scipy.pi):
			self.angle = self.angle - scipy.pi

	def __str__(self):
		return "Ellipse center: " + str(self.center) + ", major_radius: " + str(self.major_radius) + ", minor_radius: " + str(self.minor_radius) + ", angle: " + str(self.angle)

	def major_axis(self):
		return [self.major_radius*np.sin(self.angle),self.major_radius*np.cos(self.angle)]

	def minor_axis(self):
		return [self.minor_radius*np.cos(self.angle),self.minor_radius*np.sin(self.angle)]

	def is_same(self,ellipse):
		return (self.center == ellipse.center and self.major_radius == ellipse.major_radius and self.minor_radius == ellipse.minor_radius and self.angle == ellipse.angle)

	def scale(self,scale):
		self.center = [self.center[0]*scale,self.center[1]*scale]
		self.major_radius = self.major_radius*scale
		self.minor_radius = self.minor_radius*scale
		self.angle = self.angle*scale

	def pointAlongEllipse(self, theta):
		#theta is the angle
		xt = self.center[0] + self.major_radius*np.cos(self.angle)*np.cos(theta) - self.minor_radius*np.sin(self.angle)*np.sin(theta)
		yt = self.center[1] + self.major_radius*np.sin(self.angle)*np.cos(theta) + self.major_radius*np.cos(self.angle)*np.sin(theta)
		return np.array([xt,yt])

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

class Sphere:
	def __init__(self,center=[0,0,0],radius=0):
		self.center = np.array(center)
		self.radius = radius

	def __str__(self):
		return "Sphere center: " + str(self.center) + " ,radius: " + str(self.radius)

	def is_same(self,sphere):
		return (self.center == sphere.center and self.radius == sphere.radius)

if __name__ == '__main__':
	#testing if modules here work correctly
	print "yay testing"
	huding = Ellipse((-141.07,72.6412),46.0443, 34.5685, 0.658744*scipy.pi)
	hucon = Conic(huding)
	hucon.operator(None, None)