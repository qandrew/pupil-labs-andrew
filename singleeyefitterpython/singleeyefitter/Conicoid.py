"""
	Andrew Xia playing around with porting c++ code to python
	I want to port Concoid.h into this python script here
	June 25 2015

"""

import numpy as np
import Conic

class Conicoid:

	def __init__(self,A=0.0,B=0.0,C=0.0,F=0.0,G=0.0,H=0.0,U=0.0,V=0.0,W=0.0,D=0.0):
		#
		self.A = A
		self.B = B
		self.C = C
		self.D = D
		self.F = F
		self.G = G
		self.H = H
		self.U = U
		self.V = V
		self.W = W

	def __str__(self):
		#return "Conicoid, A: " + str(self.A) + ", B: " + str(self.B) + ", C: "  + str(self.C) + ", F: " + str(self.F) + ", G: " + str(self.G) + ", H: " + str(self.H) + ", U: " + str(self.U) + ", V: " + str(self.V) + ", W: " + str(self.W) + ", D: " + str(self.D)
		return "Conicoid { " + str(self.A) + "x^2 + " + str(self.B) + "y^2 + "  + str(self.C) + "z^2 + " + str(self.F) + "yz + 2*" + str(self.G) + "zx +2*" + str(self.H) + "xy + " + str(self.U) + "x + 2*" + str(self.V) + "y + 2*" + str(self.W) + "z + " + str(self.D) + " = 0 }"

	def initialize_conic(self,conic,vertex):
		alpha = vertex[0]
		beta = vertex[1]
		gamma = vertex[2]

		self.A = np.square(gamma)*conic.A
		self.B = np.square(gamma)*conic.C
		self.C = np.square(alpha)*conic.A + alpha*beta*conic.B + np.square(beta)*conic.C + conic.D*alpha + conic.E*beta + conic.F
		self.F = -gamma * (conic.C * beta + conic.B / 2 * alpha + conic.E / 2)
		self.G = -gamma * (conic.B / 2 * beta + conic.A * alpha + conic.D / 2)
		self.H = np.square(gamma)*conic.B/2
		self.U = np.square(gamma)*conic.D/2
		self.V = np.square(gamma)*conic.E/2
		self.W = -gamma * (conic.E / 2 * beta + conic.D / 2 * alpha + conic.F)
		self.D = np.square(gamma)*conic.F

	def operator(self,x,y,z):
		return self.A*np.square(x) + self.B*np.square(y) + self.C*np.square(z) + 2 * self.F*y*z + 2 * self.G*x*z + 2 * self.H*x*y + 2 * self.U*x + 2 * self.V*y + 2 * self.W*z + self.D

	def intersectZ(self,z=0):
		toreturn = Conic()
		toreturn.A = self.A
		toreturn.B = 2*self.H
		toreturn.C = self.C
		toreturn.D = 2*self.G*z + 2*self.U
		toreturn.E = 2*self.F*z + 2*self.V
		toreturn.F = self.C*np.square(z) + 2*self.W*z + self.D
		return toreturn

	def transformed(a,t):
		"""We map x,y,z to a new space using [x y z] -> affine*[x y z] + translation
			
			Using a for affine and t for translation:
				x -> a_00*x + a01*y + a02*z + t0
				y -> a_10*x + a11*y + a12*z + t1
				z -> a_20*x + a21*y + a22*z + t2
			
			So
				Ax^2 + By^2 + Cz^2 + 2Fyz + 2Gzx + 2Hxy + 2Ux + 2Vy + 2Wz + D
			becomes
				  A(a_00*x + a01*y + a02*z + t0)(a_00*x + a01*y + a02*z + t0)
				+ B(a_10*x + a11*y + a12*z + t1)(a_10*x + a11*y + a12*z + t1)
				+ C(a_20*x + a21*y + a22*z + t2)(a_20*x + a21*y + a22*z + t2)
				+ 2F(a_10*x + a11*y + a12*z + t1)(a_20*x + a21*y + a22*z + t2)
				+ 2G(a_20*x + a21*y + a22*z + t2)(a_00*x + a01*y + a02*z + t0)
				+ 2H(a_00*x + a01*y + a02*z + t0)(a_10*x + a11*y + a12*z + t1)
				+ 2U(a_00*x + a01*y + a02*z + t0)
				+ 2V(a_10*x + a11*y + a12*z + t1)
				+ 2W(a_20*x + a21*y + a22*z + t2)
				+ D
			
			Collecting terms gives:
		"""

		a = self.A*np.square(a(0, 0)) + self.B*np.square(a(1, 0)) + self.C*np.square(a(2, 0)) + 2*self.F*a(1, 0)*a(2, 0) + 2*self.G*a(0, 0)*a(2, 0) + 2*self.H*a(0, 0)*a(1, 0)
		b = self.A*np.square(a(0, 1)) + self.B*np.square(a(1, 1)) + self.C*np.square(a(2, 1)) + 2*self.F*a(1, 1)*a(2, 1) + 2*self.G*a(0, 1)*a(2, 1) + 2*self.H*a(0, 1)*a(1, 1),
		c = self.A*np.square(a(0, 2)) + self.B*np.square(a(1, 2)) + self.C*np.square(a(2, 2)) + 2*self.F*a(1, 2)*a(2, 2) + 2*self.G*a(0, 2)*a(2, 2) + 2*self.H*a(0, 2)*a(1, 2),
		f = self.A*a(0, 1)*a(0, 2) + self.B*a(1, 1)*a(1, 2) + self.C*a(2, 1)*a(2, 2) + self.F*a(1, 1)*a(2, 2) + self.F*a(1, 2)*a(2, 1) + self.G*a(0, 1)*a(2, 2) + self.G*a(0, 2)*a(2, 1) + self.H*a(0, 1)*a(1, 2) + self.H*a(0, 2)*a(1, 1),
		g = self.A*a(0, 2)*a(0, 0) + self.B*a(1, 2)*a(1, 0) + self.C*a(2, 2)*a(2, 0) + self.F*a(1, 2)*a(2, 0) + self.F*a(1, 0)*a(2, 2) + self.G*a(0, 2)*a(2, 0) + self.G*a(0, 0)*a(2, 2) + self.H*a(0, 2)*a(1, 0) + self.H*a(0, 0)*a(1, 2),
		h = self.A*a(0, 0)*a(0, 1) + self.B*a(1, 0)*a(1, 1) + self.C*a(2, 0)*a(2, 1) + self.F*a(1, 0)*a(2, 1) + self.F*a(1, 1)*a(2, 0) + self.G*a(0, 0)*a(2, 1) + self.G*a(0, 1)*a(2, 0) + self.H*a(0, 0)*a(1, 1) + self.H*a(0, 1)*a(1, 0),
		u = self.A*a(0, 0)*t(0) + self.B*a(1, 0)*t(1) + self.C*a(2, 0)*t(2) + self.F*a(1, 0)*t(2) + self.F*a(2, 0)*t(1) + self.G*a(0, 0)*t(2) + self.G*a(2, 0)*t(0) + self.H*a(0, 0)*t(1) + self.H*a(1, 0)*t(0) + self.U*a(0, 0) + self.V*a(1, 0) + self.W*a(2, 0),
		v = self.A*a(0, 1)*t(0) + self.B*a(1, 1)*t(1) + self.C*a(2, 1)*t(2) + self.F*a(1, 1)*t(2) + self.F*a(2, 1)*t(1) + self.G*a(0, 1)*t(2) + self.G*a(2, 1)*t(0) + self.H*a(0, 1)*t(1) + self.H*a(1, 1)*t(0) + self.U*a(0, 1) + self.V*a(1, 1) + self.W*a(2, 1),
		w = self.A*a(0, 2)*t(0) + self.B*a(1, 2)*t(1) + self.C*a(2, 2)*t(2) + self.F*a(1, 2)*t(2) + self.F*a(2, 2)*t(1) + self.G*a(0, 2)*t(2) + self.G*a(2, 2)*t(0) + self.H*a(0, 2)*t(1) + self.H*a(1, 2)*t(0) + self.U*a(0, 2) + self.V*a(1, 2) + self.W*a(2, 2),
		d = self.A*np.square(t(0)) + B*np.square(t(1)) + self.C*np.square(t(2)) + 2*self.F*t(1)*t(2) + 2*self.G*t(0)*t(2) + 2*self.H*t(0)*t(1) + 2*self.U*t(0) + 2*V*t(1) + 2*self.W*t(2) + self.D
		return Conicoid(a,b,c,d,f,g,h,u,v,w,d)
