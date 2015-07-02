"""
	Andrew Xia playing around with porting c++ code to python
	I want to port Ellipse.h into this python script here
	June 25 2015

"""

import numpy as np
import scipy

class Ellipse:

	def __init__(self,centre=[0,0],major_radius=0.0,minor_radius=0.0,angle=0.0, conic = None):
		self.centre = np.array(centre)
		self.major_radius = major_radius
		self.minor_radius = minor_radius
		self.angle = angle
		if conic:
			self._initialize_conic(conic)

	def _initialize_conic(self, conic):
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
		tuCentre = -Au / (2.0*Auu)
		tvCentre = -Av / (2.0*Auu)
		wCentre = Ao - Auu*np.square(tuCentre) - Avv*np.square(tvCentre)

		self.centre = [0,0]
		self.centre[0] = tuCentre*cost - tvCentre*sint
		self.centre[1] = tuCentre*sint - tvCentre*cost

		self.major_radius = np.sqrt(abs(-wCentre/Auu))
		self.minor_radius = np.sqrt(abs(-wCentre/Avv))

		if (self.major_radius < self.minor_radius):
			major_radius,minor_radius = minor_radius,major_radius
			self.angle = self.angle + scipy.pi/2

		if (self.angle > scipy.pi):
			self.angle = self.angle - scipy.pi

	def __str__(self):
		return "Ellipse center: " + str(self.centre) + ", major_radius: " + str(self.major_radius) + ", minor_radius: " + str(self.minor_radius) + ", angle: " + str(self.angle)

	def major_axis(self):
		return [self.major_radius*np.sin(self.angle),self.major_radius*np.cos(self.angle)]

	def minor_axis(self):
		return [self.minor_radius*np.cos(self.angle),self.minor_radius*np.sin(self.angle)]

	def is_same(self,ellipse):
		return (self.centre == ellipse.centre and self.major_radius == ellipse.major_radius and self.minor_radius == ellipse.minor_radius and self.angle == ellipse.angle)

	def scale(self,scale):
		self.centre = [self.centre[0]*scale,self.centre[1]*scale]
		self.major_radius = self.major_radius*scale
		self.minor_radius = self.minor_radius*scale
		self.angle = self.angle*scale

	def pointAlongEllipse(self, theta):
		#theta is the angle
		xt = self.centre[0] + self.major_radius*np.cos(self.angle)*np.cos(theta) - self.minor_radius*np.sin(self.angle)*np.sin(theta)
		yt = self.centre[1] + self.major_radius*np.sin(self.angle)*np.cos(theta) + self.major_radius*np.cos(self.angle)*np.sin(theta)
		return [xt,yt]