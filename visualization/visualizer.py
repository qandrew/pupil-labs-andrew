"""
	Andrew Xia working on visualizing data.
	I want to use opengl to display the 3d sphere and lines that connect to it.
	This file is in pupil-labs-andrew/sphere_fitter, so it is the prototype version
	July 6 2015

"""

from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

window = 0                                             # glut window number
width, height = 500, 400                               # window size

############################### THE DUMB WAY OF NOT USING A CLASS ################################
def init():
	window = 0                                             # glut window number
	width, height = 500, 400                               # window size   
	# initialization
	glutInit()                                             # initialize glut
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)
	glutInitWindowSize(width, height)                      # set window size
	glutInitWindowPosition(2000, 0)                        # set window position
	window = glutCreateWindow("noobtuts.com")              # create window with title
	glutDisplayFunc(draw)                                  # set draw function callback
	glutIdleFunc(draw)                                     # draw all the time
	glutMainLoop()                                         # start everything

def draw():                                            # ondraw is called all the time
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) # clear the screen
    glLoadIdentity()                                   # reset position
    refresh2d(width, height)
        
    glColor3f(0.0, 0.0, 1.0)                           # set color to blue
    draw_rect(10, 10, 200, 100)                        # rect at (10, 10) with width 200, height 100
    
    glutSwapBuffers()   

def refresh2d(width, height):
    glViewport(0, 0, width, height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glOrtho(0.0, width, 0.0, height, 0.0, 1.0)
    glMatrixMode (GL_MODELVIEW)
    glLoadIdentity()

def draw_rect(x, y, width, height):
    glBegin(GL_QUADS)                                  # start drawing a rectangle
    glVertex2f(x, y)                                   # bottom left point
    glVertex2f(x + width, y)                           # bottom right point
    glVertex2f(x + width, y + height)                  # top right point
    glVertex2f(x, y + height)                          # top left point
    glEnd()    

################### ABOVE IS CODE THAT WORKS BUT IS NOT WRITTEN WELL #############################

class Visualizer: #trying to make this the main class

	def __init__(self, name = "unnamed", width = 400, height = 500):
		self.name = name
		self.width = width
		self.height = height

	def run(self):
		glutInit()                                             # initialize glut
		glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)
		glutInitWindowSize(self.width, self.height)                      # set window size
		glutInitWindowPosition(2000, 0)                           # set window position
		window = glutCreateWindow(self.name)              # create window with title
		glutDisplayFunc(self.draw())                                  # set draw function callback
		glutIdleFunc(self.draw())                                     # draw all the time
		glutMainLoop()                                         # start everything

	def draw(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) # clear the screen
		glLoadIdentity()                                   # reset position
		self.refresh2d(width,height)

		glColor3f(0.0, 0.0, 1.0)                           # set color to blue
		self.draw_rect(10, 10, 200, 100)                        # rect

		glutSwapBuffers()                                  # important for double buffering

	def draw_rect(self,x,y,width,height):
		#draws a rectangle.
		glBegin(GL_QUADS)                                  # start drawing a rectangle
		glVertex2f(x, y)                                   # bottom left point
		glVertex2f(x + width, y)                           # bottom right point
		glVertex2f(x + width, y + height)                  # top right point
		glVertex2f(x, y + height)                          # top left point
		glEnd()  

	def refresh2d(self,width, height):
		glViewport(0, 0, width, height)
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		glOrtho(0.0, width, 0.0, height, 0.0, 1.0)
		glMatrixMode (GL_MODELVIEW)
		glLoadIdentity()

if __name__ == '__main__':
	print "yay local file"

	#the dumb way works.
	init()

	#Testing Visualizer class
	huding = Visualizer("huding",400,500)
	#huding.run()