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

#Pupil Tracker Stuff
from pupiltrackings import PupilTracker


#Initialize Variables (around like 300 in singleeyefitter.cpp)
ids = []
obs_eye_images = None
obs_pupil_ellipses = None #need ellipse.py
obs_pupil_inliers = None #need 


###########Step 1: Load all of the images
img_dir = '/home/andrew/singleeyefitter/img_andrew_eye/' #for now, I'll hard code the image directory
os.chdir(img_dir)

stop = 6

for image_file in glob.glob('*.png'): #hardcoded to read png files right now
	#temporary, don't want to deal with too many files for now
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

############step 2: Generate all of the cache files