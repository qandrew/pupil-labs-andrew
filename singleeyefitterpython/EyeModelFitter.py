"""
	Andrew Xia playing around with porting c++ code to python
	I want to port singgleeyefitter.cpp and singgleeyefitter.h into this python script here
	June 25 2015

	EDIT June 30 2015: This file will contain only the class EyeModelFitter, that was 
	originally in both singleeyefitter.h and singleeyefitter.cpp.

"""

#importing stuff
import numpy as np
import cv2
import scipy
from geometry import Ellipse
from geometry import Sphere
from geometry import Circle
from geometry import projection
from geometry import intersect
import singleeyefitter

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

class Observation: #was a structure in C
	def __init__(self, image = None, ellipse = Ellipse.Ellipse(), inliers = []):
		self.image = image #cv.Mat
		self.ellipse = ellipse
		self.inliers = inliers #containing cv2.Point2f

	def __str__(self):
		return "Observation: Image (file too big), {" +  str(self.ellipse) + "}, inliers (file too big)."

class PupilParams: #was a structure in C
	def __init__(self, theta = 0, psi = 0, radius = 0):
		self.theta = theta
		self.psi = psi
		self.radius = radius

class Pupil:
	def __init__(self, observation = Observation()):
		self.observation = observation
		self.circle = Circle.Circle3D()
		self.params = PupilParams()


#the class
class EyeModelFitter():

	def __init__(self, focal_length = 1030.3, region_band_width = 5, region_step_epsilon = 0.5):
		#initialize based on what was done in singleeyefitter.cpp
		self.focal_length = focal_length
		self.region_band_width = region_band_width
		self.region_step_epsilon = region_step_epsilon
		self.region_scale = 1

		self.camera_centre = np.array([0,0,0])
		self.index = []
		self.eye = Sphere.Sphere()
		self.pupils = []

		self.observation = Observation()
		self.pupilParams = PupilParams()
		self.pupil = Pupil()
		self.model_version = 0

	def add_observation_ellipse(self,pupil_ellipse):
		self.pupils.append(Observation(None,pupil_ellipse,None))

	def add_observation_int(self, image, pupil_ellipse, n_pseudo_inliers = 0):
		#add an observation by number of pseudo inliers
		pupil_inliers = []
		for i in xrange(n_pseudo_inliers):
			p = singleeyefitter.pointAlongEllipse(pupil, i*2*scipy.pi/n_pseudo_inliers)
			pupil_inliers.append([p[0],p[1]])
		return self.add_observation_vector (image, pupil_ellipse, pupil_inliers)


	def add_observation_vector(self,image,pupil_ellipse,pupil_inliers):
		if (image.channels() == 1 and image.depth == CV_8U):
			self.pupils.append(Observation(image, pupil_ellipse, pupil_inliers))
			return len(self.pupils) -1
		else:
			print "ERROR: image channels != 1 or image depth != CV_8U, terminated"
			return

	def reset():
		self.pupils = []
		eye = Sphere.Sphere()
		model_version += 1

	def circleFromParams(self, params):
		return self.circleFromParams_eye(self.eye,params)
		
	def circleFromParams_eye(self,eye,params):
		if (params.radius == 0):
			return Circle.Circle3D()

		radial = singleeyefitter.sph2cart(1, params.theta, params.psi)
		return Circle.Circle3D(eye.centre + eye.radius*radial,radial, params.radius)

	def print_single_contrast_metric_id(self,id): #don't need this functions
		self.print_single_contrast_metric(self.pupils[id])

	def print_single_contrast_metric(self, pupil): #don't need this function
		if (pupil.circle == False):
			print "No Pupil"
			return 0

		params = [pupil.params.theta, pupil.params.psi, pupil.params.radius]
		varz = [params, 0]

		contrast_term = singleeyefitter.PupilContrastTerm(eye, focal_length*region_scale,
			cv2.resize(pupil.observation.image,region_scale), region_band_width, region_step_epsilon)
		#contrast_val = contrast_term. """HUDING"""

	def unproject_observations(self,pupil_radius = 1, eye_z = 20, use_ransac = False): 
		#default to false so I skip for loop
		if (len(self.pupils) < 2):
			print "RUNTIME ERROR: Need at least two observations"
			return
		pupil_unprojection_pairs = [Circle.Circle3D(), Circle.Circle3D()]
		pupil_gazelines_proj = [] #it is a vector<line> !!

		for pupil in self.pupils:
			# get pupil circles
			# Do a per-image unprojection of the pupil ellipse into the two fixed
			# size circles that would project onto it. The size of the circles
			# doesn't matter here, only their center and normal does.
			print "huding"
			print pupil.ellipse
			print pupil_radius
			print self.focal_length
			unprojection_pair = projection.unproject(pupil.ellipse, pupil_radius, self.focal_length)

			# get projected circles and gaze vectors
			# Project the circle centers and gaze vectors down back onto the image plane.
			# We're only using them as line parameterizations, so it doesn't matter which of the two centers/gaze
			# vectors we use, as the two gazes are parallel and the centers are co-linear
			c = unprojection_pair[0].centre #it is a 3D circle
			v = unprojection_pair[0].normal
			c_proj = projection.project_point(np.reshape(c,(3,1)),self.focal_length)
			temp = v + c
			temp = np.reshape(temp,(3,1))
			v_proj = projection.project_point(temp, self.focal_length)# - c_proj
			v_proj = np.asvector() #HERE!!!
			v_proj = v_proj/np.linalg.norm(v_proj) #normalizing

			pupil_unprojection_pairs.append(unprojection_pair)
			line = intersect.Line3D(c_proj, v_proj)
			pupil_gazelines_proj.append(line)

		# Get eyeball center
		# Find a least-squares 'intersection' (point nearest to all lines) of
		# the projected 2D gaze vectors. Then, unproject that circle onto a
		# point a fixed distance away.
		# For robustness, use RANSAC to eliminate stray gaze lines
		# (This has to be done here because it's used by the pupil circle disambiguation)
		eye_centre_proj = []
		valid_eye = bool

		if (use_ransac):
			""" TO BE IMPLEMENTED (or maybe I won't bother since ransac isn't most important part"""
			pass
		else:
			for pupil in self.pupils:
				self.pupil.init_valid = True
			eye_centre_proj = intersect.nearest_intersect(pupil_gazelines_proj) #need to implement intersect!!!
			valid_eye = True

		if (valid_eye):
			eye_centre = eye_centre_proj * eye_z / focal_length, eye_z
			self.eye.radius = 1

			for i in xrange(len(self.pupils)):
				#disambiguate pupil circles using projected eyeball center
				pupil_pair = pupil_unprojection_pairs[i]
				line = pupil_gazelines_proj[i]
				c_proj = line.origin() #what?
				v_proj = line.direction() #definitely doesn't work in python...

				if (np.dot(c_proj - eye_centre_proj, v_proj) >= 0):
					#check if v_proj going away from estimated eye center, take the one going away.
					pupils[i].circle = pupil_pair[0]
				else: 
					pupils[i].circle = pupil_pair[1]
		else:
			#no inliers, so no eye
			eye = Sphere.Sphere()

			# arbitrarily pick first circle
			for i in xrange(len(self.pupils)):
				pupil_pair = pupil_unprojection_pairs[i]
				self.pupils[i].circle = pupil_pair[0]

	def initialize_model(self):
		pass

	def refine_with_region_contrast(self,CallbackFunction):
		pass

	def refine_with_inliers(self,CallbackFunction):
		pass

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


	def initialise_single_observation_id(self,id):
		self.initialise_single_observation(self.pupils[id])

	def initialise_single_observation(self,pupil): #to be implemented
		# Ignore pupil circle norma, intersect pupil circle
		# centre projection line with eyeball sphere
		try:
			"""Need to implement Intersect.py!!!!!!"""
			pupil_centre_sphere_intersect = intersect.intersect(Line)
		except:
			print "huding"

	# Local (single pupil) calculations
	def single_contrast_metric(self,id): pass
	def print_single_contrast_metric(self,id): pass
	def refine_single_with_contrast(self,id): pass



if __name__ == '__main__':

	#testing stuff
	huding = EyeModelFitter()

	#test data
	ellipse1 = Ellipse.Ellipse((-141.07,72.6412),46.0443, 34.5685, 0.658744*scipy.pi)
	ellipse2 = Ellipse.Ellipse((-134.405,98.3423),45.7818, 36.7225, 0.623024*scipy.pi)
	ellipse3 = Ellipse.Ellipse( (75.7523,68.8315),60.8489, 55.8412, 0.132388*scipy.pi)
	ellipse4 = Ellipse.Ellipse((-76.9547,52.0801),51.8554, 44.3508, 0.753157*scipy.pi)
	ellipse5 = Ellipse.Ellipse((-73.8259,5.54398),64.1682, 48.5875, 0.810757*scipy.pi)
	ellipse6 = Ellipse.Ellipse((-62.2873,-60.9237),41.1463, 23.5819, 0.864127*scipy.pi)

	huding.add_observation_ellipse(ellipse1)
	huding.add_observation_ellipse(ellipse2)
	huding.add_observation_ellipse(ellipse3)
	huding.add_observation_ellipse(ellipse4)
	huding.add_observation_ellipse(ellipse5)
	huding.add_observation_ellipse(ellipse6)

	huding.unproject_observations()

	#print huding.circleFromParams(PupilParams(1,1,10))