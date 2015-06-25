"""
	Andrew Xia playing around with porting c++ code to python
	I want to port Circle.h into this python script here
	June 25 2015

"""

import numpy as np #to store lists as mathematical arrays

class Circle3D:
	def __init__(self,centre=[0,0,0],normal=[0,0,0],radius=0):
		self.centre = np.array(centre)
		self.normal = np.array(normal)
		self.radius = radius

	def __str__(self):
		return "Circle center: " + str(self.centre) + ", normal: " + str(self.normal) + ", radius: " + str(self.radius)

	def is_same(self,circle):
		return  (self.centre == circle.centre and self.normal == circle.normal and self.radius == circle.radius)
