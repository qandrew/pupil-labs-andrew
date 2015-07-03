"""
	Andrew Xia playing around with porting c++ code to python
	I want to port singgleeyefitter.cpp and singgleeyefitter.h into this python script here
	June 25 2015

	EDIT June 30 2015: This file will contain only the class EyeModelFitter renamed spherefitter, that was 
	originally in both singleeyefitter.h and singleeyefitter.cpp.

"""

#importing stuff
import numpy as np
import cv2
import scipy

import auxiliary_functions
import geometry
import projection
import intersect
import logging
logging.info('Starting logger for...') 
logger = logging.getLogger(__name__)

#global functions
def toRotatedRect(ellipse):
	toreturn = cv2.RotatedRect(ellipse.centre,
		cv2.Size2f(2*ellipse.major_radius),
		cv2.Size2f(2*ellipse.minor_radius),
		ellipse.angle * 180 / scipy.pi)
	return toreturn

def toEllipse(rotatedrect):
	#input is a cv2 rotated rect
	toreturn = geometry.Ellipse(rect.center,
		rect.size.width/2,
		rect.size.height/2,
		rect.angle*scipy.pi/180)
	return toreturn

def circleFromParams_eye(eye,params):
	if (params.radius == None):
		logger.warning("dafaq, gimme params pls")
		return geometry.Circle3D()

	radial = auxiliary_functions.sph2cart(1, params.theta, params.psi)
	return geometry.Circle3D(eye.centre + eye.radius*radial,radial, params.radius)

class PupilParams: #was a structure in C
	def __init__(self, theta = 0, psi = 0, radius = 0):
		self.theta = theta
		self.psi = psi
		self.radius = radius

	def __str__(self):
		return "PupilParams Class: Theta" + str(self.theta) + " psi " + str(self.psi) + " r " + str(self.radius)

class Pupil: #data structure for a pupil
	def __init__(self, ellipse = geometry.Ellipse()):
		self.ellipse = ellipse
		self.circle = geometry.Circle3D()
		self.params = PupilParams()
		self.init_valid = False

	def __str__(self):
		return "Pupil Class: " + str(self.ellipse) + str(self.circle) + " " + str(self.params) + " init_valid: " + str(self.init_valid)


#the class
class Sphere_Fitter():

	def __init__(self, focal_length = 879.193, region_band_width = 5, region_step_epsilon = 0.5):
		#initialize based on what was done in singleeyefitter.cpp
		self.focal_length = focal_length
		self.region_band_width = region_band_width
		self.region_step_epsilon = region_step_epsilon
		self.region_scale = 1
		# self.index = []
		self.camera_centre = np.array([0,0,0])
		self.eye = geometry.Sphere() #model of our eye
		self.pupil_ellipse_array = [] #array containing elements in pupil class
		self.model_version = 0

	def add_observation(self,ellipse):
		#ellipse is the ellipse of pupil in camera image
		self.pupil_ellipse_array.append(Pupil(ellipse))

	def reset(self):
		self.pupil_ellipse_array = []
		eye = geometry.Sphere()
		self.model_version += 1

	def circleFromParams(self, params):
		# currently badly written
		# params = angles + diameter (theta, psi, radius)
		return circleFromParams_eye(self.eye,params)

	def unproject_observations(self,pupil_radius = 1, eye_z = 20): 
		# ransac default to false so I skip for loop (haven't implemented it yet)
		# this function for every ellipse from the image creates corresponding circles 
		# unprojected to the pupil sphere model
		if (len(self.pupil_ellipse_array) < 2):
			logger.error("Need at least two observations")
			return
		pupil_unprojection_pairs = [] #each element should be [Circle.Cirle3D, Circle.Circle3D]
		pupil_gazelines_proj = [] #it is a vector<line> !!

		for pupil in self.pupil_ellipse_array:
			""" get pupil circles
				Do a per-image unprojection of the pupil ellipse into the two fixed
				size circles that would project onto it. The size of the circles
				doesn't matter here, only their center and normal does.
			"""
			# print pupil.ellipse_Observation.ellipse
			# print pupil_radius
			# print self.focal_length
			unprojection_pair = projection.unproject(pupil.ellipse, pupil_radius, self.focal_length)

			""" get projected circles and gaze vectors
				Project the circle centers and gaze vectors down back onto the image plane.
				We're only using them as line parameterizations, so it doesn't matter which of the two centers/gaze
				vectors we use, as the two gazes are parallel and the centers are co-linear
			"""

			#why do I default use the 0th one, not the 1st one???

			#here maybe write some function that determines which line is better

			c = np.reshape(unprojection_pair[0].centre, (3,1)) #it is a 3D circle
			v = np.reshape(unprojection_pair[0].normal, (3,1))

			c_proj = projection.project_point(np.reshape(c,(3,1)),self.focal_length)
			c_proj = np.array(c_proj)
			v_proj = projection.project_point(v + c, self.focal_length) - c_proj
			v_proj = v_proj/np.linalg.norm(v_proj) #normalizing

			pupil_unprojection_pairs.append(unprojection_pair)
			line = geometry.Line2D(c_proj, v_proj)
			# print "hudong"
			# print "c " + str(c.T)
			# print "v " + str(v.T)
			# print "c_proj " + str(c_proj.T)
			# print "v_proj " + str(v_proj.T)
			pupil_gazelines_proj.append(line)

		""" Get eyeball center
			Find a least-squares 'intersection' (point nearest to all lines) of
			the projected 2D gaze vectors. Then, unproject that circle onto a
			point a fixed distance away.
			For robustness, use RANSAC to eliminate stray gaze lines
			(This has to be done here because it's used by the pupil circle disambiguation)
		"""
		eye_centre_proj = []
		valid_eye = bool

		# if (use_ransac):
		# 	""" TO BE IMPLEMENTED (or maybe I won't bother since ransac isn't most important part"""
		# 	pass
		# else:
		for pupil in self.pupil_ellipse_array:
			pupil.init_valid = True
		eye_centre_proj = intersect.nearest_intersect_2D(pupil_gazelines_proj)
		eye_centre_proj = np.reshape(eye_centre_proj,(2,))
		valid_eye = True

		if (valid_eye):
			self.eye.centre = [eye_centre_proj[0] * eye_z / self.focal_length,
				eye_centre_proj[1] * eye_z / self.focal_length, eye_z] #force it to be a 3x1 array
			self.eye.radius = 1

			for i in xrange(len(self.pupil_ellipse_array)):
				#disambiguate pupil circles using projected eyeball center
				pupil_pair = pupil_unprojection_pairs[i]
				line = pupil_gazelines_proj[i]
				c_proj = np.reshape(line.origin, (2,))
				v_proj = np.reshape(line.direction, (2,))

				if (np.dot(c_proj - eye_centre_proj, v_proj) >= 0):
					#check if v_proj going away from estimated eye center, take the one going away.
					self.pupil_ellipse_array[i].circle = pupil_pair[0]
				else: 
					self.pupil_ellipse_array[i].circle = pupil_pair[1]
				print self.pupil_ellipse_array[i].circle
		else:
			#no inliers, so no eye
			self.eye = Sphere.Sphere()

			# arbitrarily pick first circle
			for i in xrange(len(self.pupil_ellipse_array)):
				pupil_pair = pupil_unprojection_pairs[i]
				self.pupil_ellipse_yearray[i].circle = pupil_pair[0]

		self.model_version += 1
		print self.eye

""" USLESS FUNCTIONS
	def initialize_model(self): pass
	def refine_with_region_contrast(self,CallbackFunction): pass
	def refine_with_inliers(self,CallbackFunction): pass
	def initialise_single_observation_id(self,id): self.initialise_single_observation(self.pupil_ellipse_array[id])
	def unproject_single_observation(self,id, pupil_radius = 1):
		if (eye.centre[0] == 0 and eye.centre[1] == 0 and eye.centre[2]==0 and eye.radius == 0):
			print "RUNTIME ERROR: need to get eye centre estimate first (by unprojecting multiple observations)"
			return

		unprojection_pair = projection.unproject(self.pupil.observation.ellipse,pupil_radius,self.focal_length)
		c = unprojection_pair[0].centre #it is a 3D circle
		v = unprojection_pair[0].normal

		c_proj = projection.project_point(c,self.focal_length)
		v_proj = projection.project_point(v + c, self.focal_length) - c_proj
		v_proj = v_proj/np.linalg.norm(v_proj)

		eye_centre_proj = projection.project_point(eye.centre, self.focal_length)
		if (np.dot(c_proj - eye_centre_proj, v_proj) >= 0):
			self.pupil.circle = unprojection_pair[0]
		else:
			self.pupil.circle = unprojection_pair[1]
		return self.pupil.circle
	def initialise_single_observation(self,pupil): #to be implemented
		# Ignore pupil circle norma, intersect pupil circle
		# centre projection line with eyeball sphere
		try:
			#meed to implement Intersect.py!!!!!!
			pupil_centre_sphere_intersect = intersect.intersect(Line)
		except:
			print "huding"
	# Local (single pupil) calculations
	def single_contrast_metric(self,id): pass
	def print_single_contrast_metric(self,id): pass
	def refine_single_with_contrast(self,id): pass
	def print_single_contrast_metric_id(self,id): #don't need this functions
		self.print_single_contrast_metric(self.pupil_ellipse_array[id])

	def print_single_contrast_metric(self, pupil): #don't need this function
		if (pupil.circle == False):
			print "No Pupil"
			return 0

		params = [pupil.params.theta, pupil.params.psi, pupil.params.radius]
		varz = [params, 0]

		contrast_term = auxiliary_functions.PupilContrastTerm(eye, focal_length*region_scale,
			cv2.resize(pupil.observation.image,region_scale), region_band_width, region_step_epsilon)
		#contrast_val = contrast_term. 
"""

if __name__ == '__main__':

	#testing stuff
	huding = Sphere_Fitter()

	#testing unproject_observation
	ellipse1 = geometry.Ellipse((-141.07,72.6412),46.0443, 34.5685, 0.658744*scipy.pi)
	ellipse2 = geometry.Ellipse((-134.405,98.3423),45.7818, 36.7225, 0.623024*scipy.pi)
	ellipse3 = geometry.Ellipse((75.7523,68.8315),60.8489, 55.8412, 0.132388*scipy.pi)
	ellipse4 = geometry.Ellipse((-76.9547,52.0801),51.8554, 44.3508, 0.753157*scipy.pi)
	ellipse5 = geometry.Ellipse((-73.8259,5.54398),64.1682, 48.5875, 0.810757*scipy.pi)
	ellipse6 = geometry.Ellipse((-62.2873,-60.9237),41.1463, 23.5819, 0.864127*scipy.pi)

	huding.add_ellipse(ellipse1)
	huding.add_ellipse(ellipse2)
	huding.add_ellipse(ellipse3)
	huding.add_ellipse(ellipse4)
	huding.add_ellipse(ellipse5)
	huding.add_ellipse(ellipse6)

	huding.unproject_observations()
	print huding.eye

	print huding.circleFromParams(PupilParams(1,1,50))

	lol = huding.eye
	print circleFromParams_eye(lol, PupilParams(1,1,12))