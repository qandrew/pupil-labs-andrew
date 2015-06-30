"""
	Andrew Xia playing around with porting c++ code to python
	I want to port Sphere.h into this python script here
	June 25 2015

"""

import numpy as np

class Sphere:
	def __init__(self,centre=[0,0,0],radius=0):
		self.centre = np.array(centre)
		self.radius = radius

	def __str__(self):
		return "Sphere centre: " + str(self.centre) + " ,radius: " + str(self.radius)

	def is_same(self,sphere):
		return (self.centre == sphere.centre and self.radius == sphere.radius)
