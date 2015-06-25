"""
	Andrew Xia playing around with porting c++ code to python
	I want to port Projection.h into this python script here
	June 25 2015

"""

#THIS FILE IS INCOMPLETE!

import numpy as np
import scipy

import Circle
import Conic
import Sphere
import Ellipse
import Conicoid

def project_circle(circle,focal_length):
	c = circle.centre
	n = circle.normal
	r = circle.radius
	#f is the focal length

	"""INSTRUCTIONS
		Construct cone with circle as base and vertex v = (0,0,0).
		
		For the circle,
			|p - c|^2 = r^2 where (p-c).n = 0 (i.e. on the circle plane)
		
		A cone is basically concentric circles, with center on the line c->v.
		For any point p, the corresponding circle centre c' is the intersection
		of the line c->v and the plane through p normal to n. So,
		
			d = ((p - v).n)/(c.n)
			c' = d c + v
		
		The radius of these circles decreases linearly as you approach 0, so
		
			|p - c'|^2 = (r*|c' - v|/|c - v|)^2
		
		Since v = (0,0,0), this simplifies to
		
			|p - (p.n/c.n)c|^2 = (r*|(p.n/c.n)c|/|c|)^2
		
			|(c.n)p - (p.n)c|^2         / p.n \^2
			------------------- = r^2 * | --- |
				  (c.n)^2               \ c.n /
		
			|(c.n)p - (p.n)c|^2 - r^2 * (p.n)^2 = 0
		
		Expanding out p, c and n gives
		
			|(c.n)x - (x*n_x + y*n_y + z*n_z)c_x|^2
			|(c.n)y - (x*n_x + y*n_y + z*n_z)c_y|   - r^2 * (x*n_x + y*n_y + z*n_z)^2 = 0
			|(c.n)z - (x*n_x + y*n_y + z*n_z)c_z|
		
			  ((c.n)x - (x*n_x + y*n_y + z*n_z)c_x)^2
			+ ((c.n)y - (x*n_x + y*n_y + z*n_z)c_y)^2
			+ ((c.n)z - (x*n_x + y*n_y + z*n_z)c_z)^2
			- r^2 * (x*n_x + y*n_y + z*n_z)^2 = 0
		
			  (c.n)^2 x^2 - 2*(c.n)*(x*n_x + y*n_y + z*n_z)*x*c_x + (x*n_x + y*n_y + z*n_z)^2 c_x^2
			+ (c.n)^2 y^2 - 2*(c.n)*(x*n_x + y*n_y + z*n_z)*y*c_y + (x*n_x + y*n_y + z*n_z)^2 c_y^2
			+ (c.n)^2 z^2 - 2*(c.n)*(x*n_x + y*n_y + z*n_z)*z*c_z + (x*n_x + y*n_y + z*n_z)^2 c_z^2
			- r^2 * (x*n_x + y*n_y + z*n_z)^2 = 0
		
			  (c.n)^2 x^2 - 2*(c.n)*c_x*(x*n_x + y*n_y + z*n_z)*x
			+ (c.n)^2 y^2 - 2*(c.n)*c_y*(x*n_x + y*n_y + z*n_z)*y
			+ (c.n)^2 z^2 - 2*(c.n)*c_z*(x*n_x + y*n_y + z*n_z)*z
			+ (x*n_x + y*n_y + z*n_z)^2 * (c_x^2 + c_y^2 + c_z^2 - r^2)
		
			  (c.n)^2 x^2 - 2*(c.n)*c_x*(x*n_x + y*n_y + z*n_z)*x
			+ (c.n)^2 y^2 - 2*(c.n)*c_y*(x*n_x + y*n_y + z*n_z)*y
			+ (c.n)^2 z^2 - 2*(c.n)*c_z*(x*n_x + y*n_y + z*n_z)*z
			+ (|c|^2 - r^2) * (n_x^2*x^2 + n_y^2*y^2 + n_z^2*z^2 + 2*n_x*n_y*x*y + 2*n_x*n_z*x*z + 2*n_y*n_z*y*z)
		
		Collecting conicoid terms gives
		
			  [xyz]^2 : (c.n)^2 - 2*(c.n)*c_[xyz]*n_[xyz] + (|c|^2 - r^2)*n_[xyz]^2
		   [yzx][zxy] : - 2*(c.n)*c_[yzx]*n_[zxy] - 2*(c.n)*c_[zxy]*n_[yzx] + (|c|^2 - r^2)*2*n_[yzx]*n_[zxy]
					  : 2*((|c|^2 - r^2)*n_[yzx]*n_[zxy] - (c,n)*(c_[yzx]*n_[zxy] + c_[zxy]*n_[yzx]))
				[xyz] : 0
					1 : 0
	"""
	
	cn = np.dot(n,c)
	c2r2 = np.dot(c,c) - np.square(r)
	ABC = (np.square(cn) - 2.0*cn*c*n + c2r2*np.square(n))
	F = 2*(c2r2*n[1]*n[2] - cn*(n[1]*c[2] + n[2]*c[1]))
	G = 2*(c2r2*n[2]*n[0] - cn*(n[2]*c[0] + n[0]*c[2]))
	H = 2*(c2r2*n[0]*n[1] - cn*(n[0]*c[1] + n[1]*c[0]))

	conic = Conic.Conic()
	conic.A = ABC[0]
	conic.B = H
	conic.C = ABC[1]
	conic.D = G*focal_length
	conic.E = F*focal_length
	conic.F = ABC[2]*np.square(focal_length)
	return conic

def project_sphere(sphere,focal_length):
	return Ellipse.Ellipse(
			focal_length * sphere.centre / sphere.centre[2], #ellipse centre
			focal_length * sphere.radius / sphere.centre[2], #major
			focal_length * sphere.radius / sphere.centre[2], #minor
			0 #angle
		)

def project_point(point,focal_length):
	return focal_length*np.array(point)/point[2]

def unproject(ellipse,circle_radius,focal_length):
	circle = Circle.Circle3D()
	Matrix3 = np.zeros((3,3))
	RowArray3 = np.zeros((1,3))
	Translation3 = np.zeros((1,3))

	conic = Conic.Conic(ellipse)
	cam_centre_in_ellipse = np.array([[0],[0],[-focal_length]])
	pupil_cone = Conicoid.Conicoid(conic,cam_centre_in_ellipse)

	a = pupil_cone.A
	b = pupil_cone.B
	c = pupil_cone.C
	f = pupil_cone.F
	g = pupil_cone.G
	h = pupil_cone.H
	u = pupil_cone.U
	v = pupil_cone.V
	w = pupil_cone.W
	d = pupil_cone.D

	""" Get canonical conic form:

		lambda(1) X^2 + lambda(2) Y^2 + lambda(3) Z^2 = mu
		Safaee-Rad 1992 eq (6)
		Done by solving the discriminating cubic (10)
		Lambdas are sorted descending because order of roots doesn't
		matter, and it later eliminates the case of eq (30), where
		lambda(2) > lambda(1)
	"""

	#where does this solve thing in the projection.h file come from???????
	lamb = solve(1., -(a + b + c), (b*c + c*a + a*b - f*f - g*g - h*h), -(a*b*c + 2 * f*g*h - a*f*f - b*g*g - c*h*h))
	if (lamb[1] >= lamb[1]):
		print "Lambda 0 > Lambda 1, die"
		break
	if (lamb[1] > 0):
		print "Lambda 1 > 0, die"
		break
	if (lamb[2] < 0):
		print "Lambda 2 < 0, die"
		break

	#Calculate l,m,n of plane
	n = np.sqrt(lamb[1] - lamb[2])/(lamb[0]-lamb[2])
	m = 0.0
	l = np.sqrt(lamb[0] - lamb[1])/(lamb[0]-lamb[2])

	#Safaee-Rad 1992 Eq 8
	T1 = np.zeros((3,3))
	li = T1[0]
	mi = T1[1]
	ni = T1[2]
	#Safaee-Rad 1992 Eq 12
	t1 = (b - lamb)*g - f*h
	t2 = (a - lamb)*f - g*h
	t3 = -(a - lamb)*(t1/t2)/g - h/g

	#If li,mi,ni follow the left hand rule, flip their signs
	if (np.dot(np.cross(li,mi),ni) < 0):
		li = -li
		mi = -mi
		ni = -ni

	#Calculate t2


#testing
circle = Circle.Circle3D([1,1,1],[2,2,2],8)
print project_circle(circle,10)

sphere = Sphere.Sphere(centre=[0,1,2],radius=5)
ellipse = project_sphere(sphere,10)

print unproject(ellipse,5,5)