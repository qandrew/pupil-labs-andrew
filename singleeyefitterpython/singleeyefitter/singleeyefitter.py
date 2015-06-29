"""
	Andrew Xia playing around with porting c++ code to python
	I want to port singgleeyefitter.cpp and singgleeyefitter.h into this python script here
	Probably just singgleeyefitter.h
	June 25 2015

"""

#importing stuff
import numpy as np
import cv2
import scipy
import Ellipse


print "imported single eye fitter"

#global functions
def toEigen(point): pass
def toPoint2f(point): pass
def toPoint(point): pass

def toRotatedRect(ellipse):
	toreturn = cv2.RotatedRect(ellipse.centre,
		cv2.Size2f(2*ellipse.major_radius),
		cv2.Size2f(2*ellipse.minor_radius),
		ellipse.angle * 180 / scipy.pi)
	return toreturn

def toEllipse(rotatedrect):
	toreturn = Ellipse.Ellipse(rect.center,
		rect.size.width/2,
		rect.size.height/2,
		rect.angle*scipy.pi/180)
	return toreturn

#the class
class EyeModelFitter():

	def __init__(self, focal_length = 0, region_band_width = 0, region_step_epsilon = 0):

		self.focal_length = focal_length
		self.region_band_width = region_band_width
		self.region_step_epsilon = region_step_epsilon
		self.region_scale = 0

		self.camera_centre = np.array([0,0,0])
		self.index = []
		self.eye = Sphere.Sphere()
		self.pupils = []

		self.model_version = 0

	def add_ovservation_int(self, image, pupil, n_pseudo_inliers = 0):
		#add an observation by number of pseudo inliers
		pass

	def add_observation_vector(self,image,pupil,pupil_inliers):
		pass

	def reset():
		pass

	def unproject_observations(self,pupil_radius = 1, eye_z = 20, use_ransac = True):
		pass

	def initialize_model(self):
		pass

	def refine_with_region_contrast(self,CallbackFunction):
		pass

	def refine_with_inliers(self,CallbackFunction):
		pass

	def unproject_single_observation(self,id, pupil_radius = 1):
		pass

	def initialise_single_observation(self,id):
		pass

	def refine_single_with_contrast(self,id):
		pass

	def print_single_constant_metric(self,id):
		pass

	def circleFromParams(params):
		pass

	def circleFromParams_eye(eye,params):
		pass





