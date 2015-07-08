"""
	Andrew Xia working on visualizing data.
	I want to use opengl to display the 3d sphere and lines that connect to it.
	This file is in pupil-labs-andrew/sphere_fitter, so it is the prototype version
	July 6 2015

"""
import logging
from glfw import *
from OpenGL.GL import *
from OpenGL.GLUT import *

# create logger for the context of this function
logger = logging.getLogger(__name__)
from pyglui import ui

from pyglui.cygl.utils import init
from pyglui.cygl.utils import RGBA
from pyglui.cygl.utils import *
from pyglui.cygl import utils as glutils
from pyglui.pyfontstash import fontstash as fs
from trackball import Trackball

import numpy as np
import scipy
import geometry #how do I find this folder?
import cv2

def convert_fov(fov,width):
	fov = fov*scipy.pi/180
	focal_length = (width/2)/np.tan(fov/2)
	return focal_length

class Visualizer():
	def __init__(self,name = "unnamed", run_independently = False, width = 1280, height = 720, focal_length = 554.25625):
		self.name = name
		self.sphere = geometry.Sphere()
		self.ellipses = [] #collection of ellipses to display
		self.projected_lines = [] #collection of projected lines to display

		self.frame = None #video frame from eye
		self._window = None
		self.width = width
		self.height = height
		self.focal_length = focal_length
		self.input = None
		self.trackball = None
		self.run_independently = run_independently

		self.window_should_close = False

		self.video_frame = (np.linspace(0,1,num=(400*400*4))*255).astype(np.uint8).reshape((400,400,4)) #the randomized image
		self.test_sphere = geometry.Sphere([0,5,0],1)
		self.test_ellipse = geometry.Ellipse((0,3),5,3,0)

	############## DRAWING FUNCTIONS ##############################

	def draw_frustum(self,f, scale=1):
	    # average focal length
	    #f = (K[0, 0] + K[1, 1]) / 2
	    # compute distances for setting up the camera pyramid
	    W = 0.5*self.width
	    H = 0.5*self.height
	    Z = f
	    # scale the pyramid
	    W *= scale
	    H *= scale
	    Z *= scale
	    # draw it
	    glColor4f( 1, 0.5, 0, 0.5 )
	    glBegin( GL_LINE_LOOP )
	    glVertex3f( 0, 0, 0 )
	    glVertex3f( -W, H, Z )
	    glVertex3f( W, H, Z )
	    glVertex3f( 0, 0, 0 )
	    glVertex3f( W, H, Z )
	    glVertex3f( W, -H, Z )
	    glVertex3f( 0, 0, 0 )
	    glVertex3f( W, -H, Z )
	    glVertex3f( -W, -H, Z )
	    glVertex3f( 0, 0, 0 )
	    glVertex3f( -W, -H, Z )
	    glVertex3f( -W, H, Z )
	    glEnd( )

	def draw_coordinate_system(self,l=1):
		# Draw x-axis line. RED
		glLineWidth(2)
		glColor3f( 1, 0, 0 )
		glBegin( GL_LINES )
		glVertex3f( 0, 0, 0 )
		glVertex3f( l, 0, 0 )
		glEnd( )

		# Draw z-axis line. BLUE
		glColor3f( 0, 0,1 )
		glBegin( GL_LINES )
		glVertex3f( 0, 0, 0 )
		glVertex3f( 0, 0, l )
		glEnd( )

		# Draw y-axis line. GREEN. #not working... why? 
		glColor3f( 0, 1, 0 )
		glBegin( GL_LINES )
		glVertex3f( 0, 0, 0 )
		glVertex3f( 0, l, 0 )
		glEnd( )

	def draw_sphere(self,sphere):
		# this function draws the location of the eye sphere
		glPushMatrix()
		glColor3f(0.0, 0.0, 1.0)  
		glTranslate(sphere.center[0], sphere.center[1], sphere.center[2])
		glutWireSphere(sphere.radius,20,20)
		glPopMatrix()

	def draw_ellipse(self,ellipse):
		glPushMatrix()  
		glTranslate(ellipse.center[0], ellipse.center[1], 0)
		glBegin(GL_LINE_LOOP)
		for i in xrange(360):
			rad = i*2*scipy.pi/360.
			glVertex2f(np.cos(rad)*ellipse.major_radius,np.sin(rad)*ellipse.minor_radius)
		glEnd()
		glPopMatrix()

	def draw_projected_line(self,line):
		#draw a line from projected sphere center to the ellipse on frame.
		""" TO BE IMPLEMENTED """
		pass

	def draw_video_screen(self,frame):
		glPushMatrix()
		tex_id = create_named_texture(frame.shape)
		update_named_texture(tex_id,frame) #since image doesn't change, do not need to put in while loop
		draw_named_texture(tex_id)
		glPopMatrix()

	########## Setup functions I don't really understand ############

	def basic_gl_setup(self):
		glEnable(GL_POINT_SPRITE )
		glEnable(GL_VERTEX_PROGRAM_POINT_SIZE) # overwrite pointsize
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
		glEnable(GL_BLEND)
		glClearColor(.8,.8,.8,1.)
		glEnable(GL_LINE_SMOOTH)
		# glEnable(GL_POINT_SMOOTH)
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
		glEnable(GL_LINE_SMOOTH)
		glEnable(GL_POLYGON_SMOOTH)
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)

	def adjust_gl_view(self,w,h):
		"""
		adjust view onto our scene.
		"""

		glViewport(0, 0, w, h)
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		glOrtho(0, w, h, 0, -1, 1)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()

	def clear_gl_screen(self):
		glClearColor(.9,.9,0.9,1.)
		glClear(GL_COLOR_BUFFER_BIT)

	########### Open, update, close #####################

	def open_window(self):
		if not self._window:
			self.input = {'down':False, 'mouse':(0,0)}
			self.trackball = Trackball()

			# get glfw started
			if self.run_independently:
				glfwInit()
			window = glfwGetCurrentContext()					
			self._window = glfwCreateWindow(self.width, self.height, self.name, None, window)
			glfwMakeContextCurrent(self._window)

			if not self._window:
				exit()

			glfwSetWindowPos(self._window,2000,0)
			# Register callbacks window
			glfwSetFramebufferSizeCallback(self._window,self.on_resize)
			glfwSetWindowIconifyCallback(self._window,self.on_iconify)
			glfwSetKeyCallback(self._window,self.on_key)
			glfwSetCharCallback(self._window,self.on_char)
			glfwSetMouseButtonCallback(self._window,self.on_button)
			glfwSetCursorPosCallback(self._window,self.on_pos)
			glfwSetScrollCallback(self._window,self.on_scroll)
			glfwSetWindowCloseCallback(self._window,self.on_close)

			# get glfw started
			if self.run_independently:
				init()

			glutInit()
			self.basic_gl_setup()

			# self.gui = ui.UI()
			self.on_resize(self._window,*glfwGetFramebufferSize(self._window))

	def update_window(self):
		if self.window_should_close:
			self.close_window()
		if self._window != None:
			glfwMakeContextCurrent(self._window)
			self.clear_gl_screen()

			self.trackball.push()

			#THINGS I NEED TO DRAW
			self.draw_sphere(self.test_sphere) #draw the 
			self.draw_ellipse(self.test_ellipse)
			self.draw_video_screen(self.video_frame)
			self.draw_frustum(self.focal_length, scale = .01)
			self.draw_coordinate_system(20)

			self.trackball.pop()
			glfwSwapBuffers(self._window)
			glfwPollEvents()
			return True

	def close_window(self):
		if self.window_should_close == True:
			glfwDestroyWindow(self._window)
			if self.run_independently:
				glfwTerminate()
			self._window = None
			self.window_should_close = False
			logger.debug("Process done")

	############ window callbacks #################
	def on_resize(self,window,w, h):
		h = max(h,1)
		w = max(w,1)
		self.trackball.set_window_size(w,h)

		active_window = glfwGetCurrentContext()
		glfwMakeContextCurrent(window)
		self.adjust_gl_view(w,h)
		glfwMakeContextCurrent(active_window)

	def on_iconify(self,window,x,y): pass

	def on_key(self,window, key, scancode, action, mods): pass
		#self.gui.update_key(key,scancode,action,mods)

	def on_char(window,char): pass
		# self.gui.update_char(char)

	def on_button(self,window,button, action, mods):
		# self.gui.update_button(button,action,mods)
		if action == GLFW_PRESS:
			self.input['down'] = True
			self.input['mouse'] = glfwGetCursorPos(window)
		if action == GLFW_RELEASE:
			self.input['down'] = False

		# pos = normalize(pos,glfwGetWindowSize(window))
		# pos = denormalize(pos,(frame.img.shape[1],frame.img.shape[0]) ) # Position in img pixels

	def on_pos(self,window,x, y):
		hdpi_factor = float(glfwGetFramebufferSize(window)[0]/glfwGetWindowSize(window)[0])
		x,y = x*hdpi_factor,y*hdpi_factor
		# self.gui.update_mouse(x,y)
		if self.input['down']:
			old_x,old_y = self.input['mouse']
			self.trackball.drag_to(x-old_x,y-old_y)
			self.input['mouse'] = x,y

	def on_scroll(self,window,x,y):
		# self.gui.update_scroll(x,y)
		self.trackball.zoom_to(y)

	def on_close(self,window=None):
		self.window_should_close = True

if __name__ == '__main__':
 	# huding = Visualizer("huding", run_independently = True)
 	# huding.open_window()
 	# a = 0
 	# while huding.update_window():
 	# 	a += 1
 	# huding.close_window()
 	# print a

 	print convert_fov(60,640)
