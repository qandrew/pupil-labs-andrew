"""
	Andrew Xia playing around with porting c++ code to python
	I want to port PupilTracker.cpp into this python script here
	June 24 2015

"""

import cv2

class TrackerParams():

	def __init__(self):

		#initializing the variables
		self.Radius_min = None #int
		self.Radius_max = None #int

		self.CannyBlur = None #double
		self.CannyThreshold1 = None #double
		self.CannyThreshold2 = None #double
		self.StarburstPoints = None #int

		self.PercentageInliers = None #int
		self.InlierIterations = None #int
		self.ImageAwareSupport = True #boolean
		self.EarlyTerminationPercentage = None #int
		self.EarlyRejection = None
		self.seed = None

	def fildPupilEllipse(self):
		params = None #const pupiltracker::TrackerParams& params,
		m = cv.Mat()


class findPupilEllipse_out():

	def __init__(self):
		self.roiHaarPupil = None #cv2.Rect()
		self.mHaarPupil = [] #can I store it as an array?
		self.histPupil = [] 
		self.threshold = None #double
		self.mPupilThresh = [] #cv::Mat_<uchar> mHaarPupil;
		self.bbPupilThresh = None #cv2.Rect()
		self.roiPupil = None #cv2.Rect()
		self.mPupil = []
		self.mPupilOpened = []
		self.mPupilBlurred = []
		self.mPupilEdges = []
		self.mPupilSobelX = []
		self.mPupilSobelY = []

		self.edgePounts = [] #it should be a vector
		self.inliers = None #cv2.Point2f()
		self.ransacIterations = None #int
		self.earlyRejections = None #int
		self.earlyTermination = None #bool

		pPupil = None #cv2.Point2f()
		elPupil = None #cv2.RotatedRect()







####################################################################
print "imported PupilTracker"