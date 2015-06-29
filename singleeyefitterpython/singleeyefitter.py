"""
	Andrew Xia playing around with porting c++ code to python
	I want to port singgleeyefitter.cpp into this python script here
	Probably just singgleeyefitter.cpp
	June 29 2015

"""

#importing stuff
import numpy as np
import math
import cv2
import scipy
from singleeyefitter import Ellipse

def bounding_box(ellipse):
	ux = ellipse.major_radius * np.cos(ellipse.angle)
	uy = ellipse.major_radius * np.sin(ellipse.angle)
	vx = ellipse.minor_radius * np.cos(ellipse.angle + scipy.pi/2)
	vy = ellipse.minor_radius * np.sin(ellipse.angle + scipy.pi/2)
	bbox_halfwidth = np.sqrt(ux*ux + vx*vx)
	bbox_halfheight = np.sqrt(uy*uy + vy*vy)

	toreturn = cv2.rect(math.floor(ellipse.centre[0] - bbox_halfwidth),
		math.floor(ellipse.centre[1] - bbox_halfheight),
		2*math.ceil(bbox_halfwidth) + 1,
		2*math.ceil(bbox_halfheight) + 1)
	return toreturn

def norm(x,y):
	return np.sqrt(x*x + y*y)

def getXCrossing(conic,y):
	# Calculates the x crossings of a conic at a given y value. Returns the number of crossings (0, 1 or 2)
	a = Conic.A
	b = Conic.B*y + conic.D
	c = conic.C*y*y + conic.E*y + conic.F
	det = b*b - 4*a*c
	if (det == 0):
		x1 = -b/(2*a)
		return [1, [x1]] #return[0] says how many, return[1] are the values
	elif (det < 0):
		return [0]
	else: 
		rtdet = np.sqrt(det)
		return [2, [(-b-rtdet)/(2*a),(-b+rtdet)/(2*a)]]

def scaledMajorRadius(ellipse, target_radius):
	toreturn = Ellipse.Ellipse(ellipse.centre,target_radius,
		target_radius*ellipse.minor_radius/ellipse.major_radius,ellipse.angle)
	return toreturn

def ellipseGoodness (ellipse,eye,band_width, step_epsilon):	
	"""	Calculates the "goodness" of an ellipse.
		See the singleeyefitter.cpp file for more documentation

		INCOMPLETE!!!

	"""
	# band_width	the width of each band, inner and outer
	# step_epsilon	the epsilon of the soft step function

	outerEllipse = scaledMajorRadius(ellipse, ellipse.major_radius + ((band_width + step_epsilon) + 0.5))
	innerEllipse = scaledMajorRadius(ellipse, ellipse.major_radius - ((band_width + step_epsilon) + 0.5))

	sum_inner = 0
	count_inner = 0
	sum_outer = 0
	count_outer = 0

	bb = bounding_box(outerEllipse)
	bb &= cv2.Rect(-eye.cols/2,-eye.rows/2,eye.cols,eye.rows)

def angleDiffGoodness(theta1, psi1, theta2, psi2, sigma):
	if (theta1 == theta2 and psi1 == psi2):
		return 1
	temp = np.sq(np.sin((theta1-theta2)/2))
	dist = 2*np.arcsin(np.sqrt(  temp + np.cos(theta1)*np.cos(theta2)*np.sq(np.sin(psi1-psi2/2))))
	return np.exp(-np.sq(dist)/np.sq(sigma))

def circleOnSphere(sphere, theta, psi, circle_radius):
	radial = [1, theta, psi]
	return Sphere.Sphere(sphere.centre + sphere.radius * radial, radial, circle_radius)

#structure EllipseGoodnessFunction

#structure EllipsePointDistanceFunction

#not sure what these functions do/how much used.
# def Heaviside(val, epsilon):
# 	return smootherstep(-epsilon,epslion, val)

# def smootherstep(edge0, edge1, val):
# 	return 

class EllipseDistCalculator:

	def __init__(self, ellipse = None): #THIS NEEDS TO BE FIXED
		if ellipse != None:
			r = ellipse.major_radius #?????????/
			self.rA = [r*np.cos(ellipse.angle)/ellipse.major_radius, r*np.sin(ellipse.angle)/ellipse.major_radius,
				-r*np.sin(ellipse.angle)/ellipse.minor_radius, r*np.cos(ellipse.angle)/ellipse.minor_radius]
			self.rAt = self.rA*ellipse.centre
		else:
			print "ERROR: need ellipse input"
			self.rAt = None

	def calculate_ceres(self,x,y,):
		rAxt = [self.rA[0][0]*x.a + self.rA[0][1]*y.a - self.rAt[0], self.rA[0][0] * x.v + self.rA[0][1] * y.v ]
		rAyt = [self.rA[1][0]*x.a + self.rA[1][1]*y.a - self.rAt[1], self.rA[1][0] * x.v + self.rA[1][1] * y.v ]
		xy_dist = np.linalg.norm(rAxt, rAyt)
		return r - xy_dist

	def calculate(self,x,y):
		rAxt = self.rA[0][0]*x + self.rA[0][1]*y - self.rAt[0]
		rAxt = self.rA[1][0]*x + self.rA[1][1]*y - self.rAt[1]
		xy_dist = norm(rAxt, rAyt)
		return r - xy_dist

class EllipseDistanceResidualFunction:

	def __init__(self, eye_image = None, pupil_inliers = None, eye_radius = None, focal_length = None):
		self.eye_image = eye_image
		self.pupil_inliers = pupil_inliers
		self.eye_radius = eye_radius
		self.focal_length = focal_length

	def operator():
		print "huding, must be done"





if __name__ == '__main__':

	#testing uproject
	ellipse = Ellipse.Ellipse((150.254,122.157),68.3431,33.9279,0.958216*scipy.pi)
	ellipsedist = EllipseDistCalculator(ellipse)