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


width, height = (1280,720)

import numpy as np
import scipy
import geometry #how do I find this folder?
import cv2

def draw_coordinate_system(l=1):
    # Draw x-axis line. RED
    glColor3f( 1, 0, 0 )
    glBegin( GL_LINES )
    glVertex3f( 0, 0, 0 )
    glVertex3f( l, 0, 0 )
    glEnd( )

    # Draw y-axis line. GREEN. #not working... why? 
    glColor3f( 0, 1, 0 )
    glBegin( GL_LINES )
    glVertex3f( 0, 0, 0 )
    glVertex3f( 0, l, 0 )
    glEnd( )

    # Draw z-axis line. BLUE
    glColor3f( 0, 0,1 )
    glBegin( GL_LINES )
    glVertex3f( 0, 0, 0 )
    glVertex3f( 0, 0, l )
    glEnd( )

def draw_sphere(sphere):
    # this function draws the location of the eye sphere
    glPushMatrix()
    glColor3f(0.0, 0.0, 1.0)  
    glTranslate(sphere.center[0], sphere.center[1], sphere.center[2])
    glutWireSphere(sphere.radius,20,20)
    glPopMatrix()

def draw_ellipse(ellipse):
    glPushMatrix()
    # glColor3f(0.0, 1., 0)  
    # pts = cv2.ellipse2Poly( (int(ellipse.center[0]),
    #             int(ellipse.center[1])),
    #             (int(ellipse.major_radius/2),
    #             int(ellipse.minor_radius/2)),
    #             int(ellipse.angle),0,360,15)
    # for i in xrange(len(pts)-1):
    #     glutils.draw_polyline((pts[i],pts[i+1]))    
    glTranslate(ellipse.center[0], ellipse.center[1], 0)
    glBegin(GL_LINE_LOOP)
    for i in xrange(360):
        rad = i*2*scipy.pi/360.
        glVertex2f(np.cos(rad)*ellipse.major_radius,np.sin(rad)*ellipse.minor_radius)
    glEnd()
    glPopMatrix()

def draw_projected_line(line):
    #draw a line from projected sphere center to the ellipse on frame.
    pass

def draw_video_screen(frame):
    glPushMatrix()
    tex_id = create_named_texture(frame.shape)
    update_named_texture(tex_id,frame) #since image doesn't change, do not need to put in while loop
    draw_named_texture(tex_id)
    glPopMatrix()

def basic_gl_setup():
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

def adjust_gl_view(w,h):
    """
    adjust view onto our scene.
    """

    glViewport(0, 0, w, h)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glOrtho(0, w, h, 0, -1, 1)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()

def clear_gl_screen():
    glClearColor(.9,.9,0.9,1.)
    glClear(GL_COLOR_BUFFER_BIT)

def demo(video_frame = None):
    user_input = {'down':False, 'mouse':(0,0)}

    # Callback functions
    def on_resize(window,w, h):
        h = max(h,1)
        w = max(w,1)
        gui.update_window(w,h)
        active_window = glfwGetCurrentContext()
        glfwMakeContextCurrent(window)
        # norm_size = normalize((w,h),glfwGetWindowSize(window))
        # fb_size = denormalize(norm_size,glfwGetFramebufferSize(window))
        adjust_gl_view(w,h)
        glfwMakeContextCurrent(active_window)
        track.set_window_size(w,h)

    def on_iconify(window,iconfied):
        pass

    def on_key(window, key, scancode, action, mods):
        gui.update_key(key,scancode,action,mods)

    def on_char(window,char):
        gui.update_char(char)

    def on_button(window,button, action, mods):
        gui.update_button(button,action,mods)
        if action == GLFW_PRESS:
            user_input['down'] = True
            user_input['mouse'] = glfwGetCursorPos(window)
        if action == GLFW_RELEASE:
            user_input['down'] = False

        # pos = normalize(pos,glfwGetWindowSize(window))
        # pos = denormalize(pos,(frame.img.shape[1],frame.img.shape[0]) ) # Position in img pixels

    def on_pos(window,x, y):
        hdpi_factor = float(glfwGetFramebufferSize(window)[0]/glfwGetWindowSize(window)[0])
        x,y = x*hdpi_factor,y*hdpi_factor
        gui.update_mouse(x,y)
        if user_input['down']:
            old_x,old_y = user_input['mouse']
            track.drag_to(x-old_x,y-old_y)
            user_input['mouse'] = x,y

    def on_scroll(window,x,y):
        gui.update_scroll(x,y)
        track.zoom_to(y)

    # get glfw started
    glfwInit()
    window = glfwCreateWindow(width, height, "pyglui demo", None, None)
    glfwMakeContextCurrent(window)

    if not window:
        exit()

    glfwSetWindowPos(window,0,0)
    # Register callbacks window
    glfwSetFramebufferSizeCallback(window,on_resize)
    glfwSetWindowIconifyCallback(window,on_iconify)
    glfwSetKeyCallback(window,on_key)
    glfwSetCharCallback(window,on_char)
    glfwSetMouseButtonCallback(window,on_button)
    glfwSetCursorPosCallback(window,on_pos)
    glfwSetScrollCallback(window,on_scroll)

    init()
    glutInit()
    basic_gl_setup()

    gui = ui.UI()
    track = Trackball()
    on_resize(window,*glfwGetFramebufferSize(window))

    if video_frame == None:
        video_frame = (np.linspace(0,1,num=(400*400*4))*255).astype(np.uint8).reshape((400,400,4)) #the randomized image

    test_sphere = geometry.Sphere([0,5,0],1)
    test_ellipse = geometry.Ellipse((0,3),5,3,0)
    #return window

#def demo_update(window):
    while not glfwWindowShouldClose(window):
        clear_gl_screen()

        track.push()

        #glutils.draw_polyline3d([(0,0,2),(1,1,2)],color=RGBA(0.4,0.5,0.3,0.5))
        draw_sphere(test_sphere)
        draw_ellipse(test_ellipse)
        draw_video_screen(video_frame)
        draw_coordinate_system(20)

        track.pop()
        glfwSwapBuffers(window)
        glfwPollEvents()

#def demo_close(window):
    glfwDestroyWindow(window)
    glfwTerminate()
    logger.debug("Process done")

if __name__ == '__main__':
 #   w = demo()
  #  w = demo_update(w)
   # demo_close(w)
    demo()