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
import geometry
import logging
logging.info('Starting logger for...') 
logger = logging.getLogger(__name__)

def smootherstep(edge0,edge1,x):
	if (x >= edge1):
		return 1
	elif (x <= edge0):
		return 0
	else:
		x = (x- edge0)/(edge1 - edge0)
		return x*x*x*(x*(x*6.0 - 15) +10 ) #WHAT IS T in C++? I think it's just a auto type thing

def Heaviside(val, epsilon):
	return smootherstep(-epsilon,epsilon,val)

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
	toreturn = geometry.Ellipse(ellipse.centre,target_radius,
		target_radius*ellipse.minor_radius/ellipse.major_radius,ellipse.angle)
	return toreturn

def ellipseGoodness (ellipse,eye,band_width, step_epsilon, scalar_tag):	 #or what exactly is scalar_tag??
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


def sph2cart(r,theta,psi):
	toreturn = np.matrix([np.sin(theta)*np.cos(psi), np.cos(theta), np.sin(theta)*np.sin(psi)]) 
	#np is column major, so initializing this should give a 3x1 matrix (or vector)
	return r*toreturn

def angleDiffGoodness(theta1, psi1, theta2, psi2, sigma):
	if (theta1 == theta2 and psi1 == psi2):
		return 1
	temp = np.sq(np.sin((theta1-theta2)/2))
	dist = 2*np.arcsin(np.sqrt(  temp + np.cos(theta1)*np.cos(theta2)*np.sq(np.sin(psi1-psi2/2))))
	return np.exp(-np.sq(dist)/np.sq(sigma))

def circleOnSphere(sphere, theta, psi, circle_radius):
	radial = [1, theta, psi]
	return geometry.Circle3D(sphere.centre + sphere.radius*np.asmatrix(radial), radial, circle_radius)

#structure EllipsePointDistanceFunction

class EllipseDistCalculator:

	def __init__(self, ellipse = None): #THIS NEEDS TO BE FIXED
		if ellipse != None:
			r = ellipse.major_radius #?????????/
			self.rA = np.matrix([[r*np.cos(ellipse.angle)/ellipse.major_radius, 
				r*np.sin(ellipse.angle)/ellipse.major_radius],
				[-r*np.sin(ellipse.angle)/ellipse.minor_radius, 
				r*np.cos(ellipse.angle)/ellipse.minor_radius]]).T #np is row major, not column major :P
			self.rAt = self.rA*np.reshape(ellipse.centre,(2,1))
		else:
			logger.error("need ellipse input")
			self.rA = None
			self.rAt = None

	def operator(x,y):
		return calculate(x,y)

	def calculate_ceres(self,x,y,):
		rAxt = [self.rA[0][0]*x.a + self.rA[0][1]*y.a - self.rAt[0], self.rA[0][0] * x.v + self.rA[0][1] * y.v ]
		rAyt = [self.rA[1][0]*x.a + self.rA[1][1]*y.a - self.rAt[1], self.rA[1][0] * x.v + self.rA[1][1] * y.v ]
		xy_dist = [rAxt, rAyt]
		#xy_dist = np.linalg.norm(rAxt, rAyt) figure out normalizing later, but not sure if this function is used
		return r - xy_dist

	def calculate(self,x,y):
		rAxt = self.rA[0][0]*x + self.rA[0][1]*y - self.rAt[0]
		rAxt = self.rA[1][0]*x + self.rA[1][1]*y - self.rAt[1]
		xy_dist = norm(rAxt, rAyt)
		return r - xy_dist

def EllipseGoodnessFunctionOperator(self,eye,theta,psi,pupil_radius,focal_length, band_width, step_epsilon, mEye):
	#originally a structure in the cpp code but I might as well just use a function call

	self.camera_centre = np.matrix([0,0,0])		

	if (pupil_radius <= 0):
		# pupil radius must be positive
		# Return -255 for radius == 0, and even lower values for
		# radius < 0
		# This should push the gradient towards positive radius,
		# rather than just returning flat -255
		return -255 + pupil_radius

	pupil_circle = circleOnSphere(eye,theta,psi,pupil_radius)
	normalDotPos = np.dot(pupil_circle.normal,camera_centre - pupil_radius.centre)
	
	if (normalDotPos <= 0):
		# Ellipse normal must point towards camera
		# Return -255 for normalDotPos == 0, and even lower values for
		# normalDotPos < 0
		# This should push the gradient towards positive normalDotPos,
		# rather than just returning flat -255
		return -255 + normalDotPos

	# Angles should be in the range
	# theta: 0 -> pi
	# psi: -pi -> 0
	# If we're outside of this range AND radialDotEye > 0, then we must
	# have gone all the way around, so just return worst case (i.e as bad
	# as radialDotEye == -1) with additional penalty for how far out we
	# are, again to push the gradient back inwards.
	if (theta < 0 or theta > scipy.pi or psi < -scipy.pi or psi > 0):
		normalized = (camera_centre - pupil_circle.centre)/np.linalg.norm(camera_centre - pupil_circle.centre)
		ret = -255 - normalized
		if (theta < 0):
			ret += theta
		elif (theta > scipy.pi):
			ret -= theta - scipy.pi
		if (psi < -scipy.pi):
			ret -= -scipy.pi - psi
		elif (psi > 0):
			ret -= psi

	# now that everything should look good, calculate the actual goodness.
	pupil_ellipse = geometry.Ellipse()
	pupil_ellipse.initialize_conic(projection.project_sphere(pupil_circle,focal_length))
	return ellipseGoodness(pupil_ellipse, mEye, band_width, step_epsilon)

class EllipseDistanceResidualFunction:

	def __init__(self, eye_image = None, pupil_inliers = None, eye_radius = None, focal_length = None):
		self.eye_image = eye_image
		self.pupil_inliers = pupil_inliers
		self.eye_radius = eye_radius
		self.focal_length = focal_length

	def operator(self,eye_param, pupil_param,e_array):
		Const = None
		eye_pos = np.matrix([eye_param[0],eye_param[1],eye_param[2]])
		eye = geometry.Sphere(eye_pos, eye_radius)

		pupil_ellipse = geometry.Ellipse()
		temp = projection.project_sphere(circleOnSphere(eye, pupil_param[0],pupil_param[1],pupil_param[2]), focal_length)
		pupil_ellipse.initialize_conic(temp)

		for i in xrange(pupil_inliers.size()):
			inlier = cv.Point2f(pupil_inliers[i])
			e_array[i] = ellipDist(inlier.x, inlier.y) #e_array will be modified by the function

		return True

class EllipsePointDistanceFunction:
	def __init__(self,ellipse = 0, x=0, y=0):
		self.ellipse = ellipse 
		self.x = x
		self.y = y

	def operator(t, e):
		pt = self.ellipse.pointAlongEllipse(t[0])
		temp = [x - pt[0],y - pt[1]]
		e[0] = temp.np.linalg.norm(temp) #pt[0] is x, pt[1] is y. See Ellipse.py for further info
		return True

#class PupilContrastTerm(spii.Term): #what do I want to do about this??
class PupilContrastTerm(): #DONT DO
	def __init__(self, init_eye = None, focal_length = None, eye_image = None, band_width = None, has_eye_var = True):
		self.has_eye_var = has_eye_var

		self.init_eye = init_eye
		self.focal_length = focal_length
		self.eye_image = eye_image
		self.band_width = band_width
		self.step_epsilon = step_epsilon

	def change_has_eye_var(self):
		self.has_eye_var = not self.has_eye_var

	def eye_var_idx(self): return 0 if self.has_eye_var else -1
	def pupil_var_index(self): return 1 if self.has_eye_var else 0	
	def number_of_variables(self):	return 2 if self.has_eye_var else 1

	def evaluate(self,varz):
		pupil_vars = varz[self.pupil_var_idx()]
		eye = self.init_eye
		if (has_eye_var):
			eye_vars = varz[self.eye_var_idx()]
			eye.centre = [eye_vars[0],eye_vars[1],eye_vars[2]]
		
		#goodnessFunction = ellipseGoodness



#class PupilAnthroTerm(spii.Term):
class PupilAnthroTerm:
	def __init__(self, mean = 0., sigma = 0., scale = 0.):
		self.mean = mean
		self.sigma = sigma
		self.scale = scale

	def number_of_variables(self): return 1
	def variable_dimension(self,var): return 3 if var == 0 else -1

	def evaluate(vars):
		r = vars[0][2]
		goodness = np.exp(-np.sq(r - self.mean)/np.sq(self.sigma))
		return -goodness*self.scale #the cost

	def evaluate_with_gradient(self, vars, gradient):
		pass

	def evaluate_with_hessian(self, variables, gradient, hessian):
		#NOT IMPLEMENTED
		pass

if __name__ == '__main__':

	#testing stuff
	smootherstep(1,2,1.5)

	ellipse = geometry.Ellipse((150.254,122.157),68.3431,33.9279,0.958216*scipy.pi)
	ellipsedist = EllipseDistCalculator(ellipse)

	sphere = geometry.Sphere([0,0,0],10)
	print circleOnSphere(sphere, 1,2,3)