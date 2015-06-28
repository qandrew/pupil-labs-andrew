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
import solve

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
			[focal_length * sphere.centre[0] / sphere.centre[2],focal_length * sphere.centre[1] / sphere.centre[2]], #ellipse centre
			focal_length * sphere.radius / sphere.centre[2], #major
			focal_length * sphere.radius / sphere.centre[2], #minor
			0 #angle
		)

def project_point(point,focal_length):
	return [focal_length*np.array(point)[0]/point[2],focal_length*np.array(point)[1]/point[2]]

def unproject(ellipse,circle_radius,focal_length):
	circle = Circle.Circle3D()
	Matrix3 = np.zeros((3,3))
	RowArray3 = np.zeros((1,3))
	Translation3 = np.zeros((1,3)) #see T2 for actual implementation

	conic = Conic.Conic(ellipse)
	cam_centre_in_ellipse = np.array([[0],[0],[-focal_length]])
	pupil_cone = Conicoid.Conicoid()
	pupil_cone.initialize_conic(conic,cam_centre_in_ellipse) #this step is fine

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

	"""FIX HERE IT DOESNT WORK"""
	lamb = solve.solve_four(1., 
		-(a + b + c), 
		(b*c + c*a + a*b - f*f - g*g - h*h), 
		-(a*b*c + 2 * f*g*h - a*f*f - b*g*g - c*h*h) )
	print lamb
	if (lamb[0] < lamb[1]):
		print "Lambda 0 > Lambda 1, die"
		return
	if (lamb[1] <= 0):
		print "Lambda 1 > 0, die"
		return
	if (lamb[2] >= 0):
		print "Lambda 2 < 0, die"
		return

	#Calculate l,m,n of plane
	n = np.sqrt(lamb[1] - lamb[2])/(lamb[0]-lamb[2])
	m = 0.0
	l = np.sqrt(lamb[0] - lamb[1])/(lamb[0]-lamb[2])

	# print "n: " + str(n) + ", m: " + str(m) + ", l:" + str(l)

	#Safaee-Rad 1992 Eq 8

	#Safaee-Rad 1992 Eq 12
	t1 = (b - lamb)*g - f*h
	t2 = (a - lamb)*f - g*h
	t3 = -(a - lamb)*(t1/t2)/g - h/g

	mi = 1 / np.sqrt(1 + np.square(t1 / t2) + np.square(t3))
	li = (t1 / t2) * mi
	ni = t3 * mi

	print "t1: " +str(t1)
	print "t2: " +str(t2)
	print "t3: " +str(t3)
	print "li: " +str(li) 
	print "mi: " +str(mi)
	print "ni: " +str(ni)

	#If li,mi,ni follow the left hand rule, flip their signs
	li = np.reshape(li,(3,))
	mi = np.reshape(mi,(3,))
	ni = np.reshape(ni,(3,))

	if (np.dot(np.cross(li,mi),ni) < 0):
		li = -li
		mi = -mi
		ni = -ni

	T1 = np.zeros((3,3))
	T1[:,0] = li
	T1[:,1] = mi
	T1[:,2] = ni

	print T1

	#Calculate t2 a translation transformation from the canonical
	#conic frame to the image space in the canonical conic frame
	#Safaee-Rad 1992 eq (14)

	"""3 or 4???"""
	T2 = np.identity(3) #the transformation matrix
	temp = -(u*li + v*mi + w*ni) / lamb
	print temp
	T2[0,2] = temp[0]
	T2[1,2] = temp[1]
	T2[2,2] = temp[2]

	# print "T2: " + str(T2)

	solutions = Circle.Circle3D()
	ls = [l, -l]

	solutions = [0,0] #two solutions for the circles that we will return

	for i in (0,1):
		l = ls[i]
		# T1 = [[-0.0343973,0.988934,-0.144312],
  # 			[0.987195,0.0561263,0.149318],
  # 			[-0.155766,0.137328,0.978201]]

		gaze = T1 * np.matrix([[l],[m],[n]])
		# print gaze

		#calculate t3, rotation from frame where Z is circle normal

		T3 = np.zeros((3,3))
		if (l == 0):
			if (n == 1):
				print "Warning: l == 0"
				break
			T3 = np.matrix([[0,-1,0],
				[1,0,0],
				[0,0,1]])
		else:
			T3 = np.matrix([[0,-n*np.sign(l),0], 
				[np.sign(l),0,0],
				[0,abs(l),n]])

		#calculate circle center 
		#Safaee-Rad 1992 eq (38), using T3 as defined in (36)
		lamb =  np.reshape(lamb,(3,))
		T30 = abs(np.array([T3[0,0],T3[0,1],T3[0,2]]))
		T31 = np.cross(np.array([T3[0,0],T3[0,1],T3[0,2]]), np.array([T3[2,0],T3[2,1],T3[2,2]]))
		T32 = np.cross(np.array([T3[1,0],T3[1,1],T3[1,2]]), np.array([T3[2,0],T3[2,1],T3[2,2]]))
		T33 = abs(np.array([T3[2,0],T3[2,1],T3[2,2]]))

		A = np.dot(lamb,T30)
		B = np.dot(lamb,T31)
		C = np.dot(lamb,T32)
		D = np.dot(lamb,T33)

		# Safaee-Rad 1992 eq 41
		centre_in_Xprime = np.zeros((3,1))
		centre_in_Xprime[2] = A*circle_radius/ np.sqrt(np.square(B) + np.square(C) - A*D)
		centre_in_Xprime[0] = -B / A * centre_in_Xprime[2]
		centre_in_Xprime[1] = -C / A * centre_in_Xprime[2]

		# Safaee-Rad 1992 eq 34
		"""3 or 4???"""
		T0 = np.matrix(np.identity(3)) #the transformation matrix
		T0[2,2] = focal_length

		# Safaee-Rad 1992 eq 42 using eq 35
		centre = T0*T1*T2*T3*centre_in_Xprime

		if (centre[2] < 0):
			centre_in_Xprime = -centre_in_Xprime
			centre = T0*T1*T2*T3*centre_in_Xprime #make sure z is positive

		gaze = np.reshape(gaze,(3,))

		if (np.dot(gaze,centre) > 0):
			gaze = -gaze
		np.linalg.norm(gaze)

		centre = np.reshape(centre,3)

		solutions[i] = Circle.Circle3D(centre,gaze,circle_radius)

	#print T3
	return solutions

if __name__ == '__main__':

	#testing uproject
	ellipse = Ellipse.Ellipse((-152.295,157.418),46.7015,32.4274,0.00458883*scipy.pi)
	huding = unproject(ellipse,1,1030.3) 
	print huding[0] #should get { c: (1.06199,-0.715622,7.39044), n: (-0.697086,0.603144,-0.38767), r: 1 }

	# testing project_circle
	# circle = Circle.Circle3D([1.35664,-0.965954,9.33736],[0.530169,-0.460575,-0.711893],1)
	# print project_circle(circle,1030.3) #correct
	# circle = Circle.Circle3D([-1.36026,0.415986,4.9436],[-0.701632,0.102242,-0.705166],1)
	# print project_circle(circle,1030.3) #correct
	# circle = Circle.Circle3D((0.68941,0.58537,4.71737),(0.0464665,0.794044,-0.606082),1)
	# print project_circle(circle,1030.3) #correct
	# circle = Circle.Circle3D([-11.3327,11.8036,76.7857],[0.0843856,0.548524,-0.831865],5.07561)
	# print project_circle(circle,1030.3) #correct
	# circle = Circle.Circle3D([-21.9309,6.70676,79.7034],[-0.798796,0.123788,-0.588729],16.1225)
	# print project_circle(circle,1030.3) #correct

	#testing project_point
	# point = [0.493976,-0.376274,4.35446]
	# print project_point(point,1030.3) #good

	#testing project_sphere
	# sphere = Sphere.Sphere(centre=[-12.3454,5.22129,86.7681],radius=12)
	# print project_sphere(sphere,1030.3) #GOOD
	# sphere = Sphere.Sphere(centre=(10.06,-6.20644,86.8967),radius=12)
	# print project_sphere(sphere,1030.3) #GOOD