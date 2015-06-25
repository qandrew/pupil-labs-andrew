"""
	Andrew Xia playing around with porting c++ code to python
	I want to port singleeyefitter.cpp into this python script here
	June 24 2015

"""


###########################
#The actual Code here

#Importing Stuff
import numpy as np #needed for math calculations
import cv2 #open source computer vision
import glob, os #for scanning through directories
import boost #boost provides free peer reviewed portable C++ source library
import minieigen

#Pupil Tracker Stuff
from pupiltrackings import PupilTracker
from singleeyefitter import singleeyefitter


#Initialize Variables (around like 300 in singleeyefitter.cpp)
ids = []
obs_eye_images = []
obs_pupil_ellipses = [] #need ellipse.py
obs_pupil_inliers = [] #need 


###########Step 1: Load all of the images, generate cache files

img_dir = '/home/andrew/singleeyefitter/img_andrew_eye/' #for now, I'll hard code the image directory
os.chdir(img_dir)

stop = 6

for image_file in glob.glob('*.png'): #hardcoded to read png files right now
	#temporary, don't want to deal with too many files for now
	print image_file
	if (stop == 0):
		break
	stop -= 1

	cache_valid = False #we assume there isnt a cache file or the current cache file is wrong
	cache_version = 8
	cache_path = image_file[:-4] + ".cache"

	if (os.path.isfile(cache_path)):
		read_cache = open(cache_path,'r')
		for line in read_cache:
			ia_cache_version = int(line[29]) #i feel like this is unsafe

		if (cache_version == ia_cache_version):
			print "Loading: " + image_file
			cache_valid = True

			ellipse = line
			inlier_pts = line
		

	if (cache_valid == False):
		print "Pupil Tracking " + image_file

		pupil_tracker_params = PupilTracker.TrackerParams()
        pupil_tracker_params.Radius_Min = 20
        pupil_tracker_params.Radius_Max = 70
        pupil_tracker_params.CannyThreshold1 = 20
        pupil_tracker_params.CannyThreshold2 = 40
        pupil_tracker_params.CannyBlur = 1.6
        pupil_tracker_params.EarlyRejection = True
        pupil_tracker_params.EarlyTerminationPercentage = 95
        pupil_tracker_params.PercentageInliers = 20
        pupil_tracker_params.InlierIterations = 2
        pupil_tracker_params.ImageAwareSupport = True
        pupil_tracker_params.StarburstPoints = 0

        pupil_tracker_out = PupilTracker.findPupilEllipse_out()

        """
        #code below has not been confirmed to work
        found = PupilTracker.findPupilEllipse(pupil_tracker_params, eye, pupil_tracker_out, log)

        if (found):
        	el = pupil_tracker_out.elPupil
        	el.centre -= Eigen

        	for inlier in pupil_tracker_out.inliers:
        		inlier_pts.push_back(cv2.Point2f(pupil_tracker_out.roiPupil.x + inlier.x - eye.cols/2, pupil_tracker_out.roiPupil.y + inlier.y - eye.rows/2))
        else:
        	el = None
        read_cache.write(el)
        read_cache.close()
        """
        el = True
        if (el):
        	ids.append(int(image_file[11:-4]))
        	obs_eye_images.append(eye)
        	obs_pupil_ellipses.append(el)
        	obs_pupil_inliers.append(move(inlier_pts))
        	print "Pupil detected"

	sorted(ids)


""" One way of printing all of the images

	for image_file in glob.glob('*.png'): #hardcoded to read png files right now
		pngdata = cv2.imread(image_file,cv2.CV_LOAD_IMAGE_GRAYSCALE)
		cv2.imshow('img',pngdata)
		cv2.waitKey(1)
		#print image_file

	cv2.destroyAllWindows()
"""

""" Another way of printing all the images
	i = 1
	while (True):
		pngdata = cv2.imread('/home/andrew/singleeyefitter/img_andrew_eye/andrew_eye_%000i.png'%i,0)
		cv2.imshow('img',pngdata)
		cv2.waitKey(1)
		#i += 1		
	cv2.destroyAllWindows()
"""

############step 2: THE COOL STUFF