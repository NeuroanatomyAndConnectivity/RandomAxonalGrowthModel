# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 17:54:28 2015

@author: Jan Zelmer
@mail: jan.zelmer@gmx.net
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 17:33:22 2015

@author: Jan Zelmer
@mail: jzelmer@gmx.net

"""
from math import sqrt
from math import radians
from math import cos
from math import sin
from math import pow
from math import e
from math import pi
from copy import deepcopy
import random
import numpy as np

random.seed()


class n3Point(object):
	
	def __init__ (self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
		self.length = sqrt (x * x + y * y + z * z)

	def __add__ (self, p):
		return self.__class__(self.x + p.x, self.y + p.y, self.z + p.z)

	def __sub__ (self, p):
		return self.__class__(self.x - p.x, self.y - p.y, self.z - p.z)
		
	def __mul__ (self, k):
		return self.__class__(self.x * k, self.y * k, self.z * k)
		
	def __div__ (self, k):
         i = float(k)
         return self.__class__(self.x / i, self.y / i, self.z / i)
		
	def __truediv__ (self, k):
         i = float(k)
         return self.__class__(self.x / i, self.y / i, self.z / i)
		
	def scalar_product (self, v):
		return self.x * v.x + self.y * v.y + self.z * v.z	
	
	def vector_product (self, v):
		return self.__class__(self.y * v.z - self.z * v.y, self.z * v.x - self.x * v.z, self.x * v.y - self.y * v.x)

	def distance_to (self, p):
		return (self - p).length
		
	def __str__(self):
		""" Returns a readable string.
		"""
		return str("(" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")")

	def __repr__(self):
		""" Returns a readable string.
		"""
		return str("(" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")")
		
	def rotate_by_degree (self, alpha, beta, gamma):
		""" Drehung im Uhrzeigersinn
		"""
		alpha = radians (alpha)
		beta = radians (beta)
		gamma = radians (gamma)		
		y1 = self.y * cos (alpha) - self.z * sin (alpha)
		z1 = self.y * sin (alpha)	+ self.z * cos (alpha)		
		x1 = self.x * cos (beta) + z1 * sin (beta)		
		z2 = round(z1 * cos (beta) - self.x * sin (beta), 4)
		x2 = round(x1 * cos(gamma) - y1 * sin (gamma), 4)
		y2 = round(x1 * sin(gamma) + y1 * cos (gamma), 4)
		return self.__class__(x2, y2, z2)
		
		

class n3Line(object):
	
	def __init__(self, *points):
		""" Simple constructor. Pushes the line by iterating over all points
		provided. So the last point is the head and the one before last is 
		the tail.
		"""		
		self.head = points[0]
		for p in points[1:]:
			self.push(p)
		
	def push (self, p): 
		""" Pushes the line by adding the point p. The new head is point p 
			and the new tail is the old head.
		"""
		self.tail = self.head
		self.head = p
		self.direction_vector = self.head - self.tail
		self.direction_vector_norm = self.direction_vector / self.direction_vector.length
		self.middle = self.tail + (self.head - self.tail) / 2 

	def push_delta (self, p):
		""" Pushes the line by adding the coordinates form the point p to the
			actual head as the new head. Finally setting the tail to the old
			head.
		"""
		self.tail = self.head
		self.head = self.head + p
		self.direction_vector = self.head - self.tail
		self.direction_vector_norm = self.direction_vector / self.direction_vector.length
		self.middle = self.tail + (self.head - self.tail) / 2 

	def get_middle (self):
		""" Returns the point representing the middle.
		"""
		return self.middle

	def get_head (self):
		""" Returns the point representing the head.
		"""
		return self.head
		
	def get_tail (self):
		""" Returns the point representing the tail.
		"""
		return self.tail

	def compute_distance (self, p):		
		""" Computes the distance between a point p and this line. It is 
			basically the minimum distance over all the points from the 
			line and the point p. This is done by computing the perpendicular
			and then using euclid.		
		"""		
		v = (self.tail - p) * -1 
		r = self.direction_vector.scalar_product(v) / self.direction_vector.length
		if r < 0:
			return self.tail.distance_to(p)
		elif r > self.direction_vector.length:
			return self.head.distance_to(p)
		else:
			return p.distance_to(self.direction_vector_norm * r + self.tail)

	def get_length (self):
		""" Returns the length of this line
		"""
		return (self.head - self.tail).length
 
 
class n3Sphere (object):
	
	def __init__ (self, point, radius):
		self.middle = point
		self.radius = radius
		self.longest_in_area_line_length = radius * 2
		self.volume = 4 / 3 * pi * radius * radius * radius
		
	def lies_inside (self, p):
		return self.middle.distance_to(p) <= self.radius
		
	def get_bounding_box (self):
		return n3AxisParallelRectangle(n3Point(self.middle.x - self.radius, self.middle.y - self.radius, self.middle.z - self.radius), \
		n3Point(self.middle.x + self.radius, self.middle.y + self.radius, self.middle.z + self.radius))
		
	def translate (self, p):
		self.middle += p
		
	def rotate_by_root (self, alpha, beta, gamma):
		self.middle = self.middle.rotate_by_degree(alpha, beta, gamma)
    
	def rotate_by_self (self, alpha, beta, gamma):
		print "Rotating a sphere with the axis aligned on it's center won't do anything."
				
class n3Ellipsoid (object):
	
	def __init__ (self, middle, rx, ry, rz):
		self.middle = middle
		self.rx = rx
		self.ry = ry
		self.rz = rz
		self.alpha = 0
		self.beta = 0
		self.gamma = 0
		self.volume = 4 / 3 * pi * rx * ry * rz
		self.longest_in_area_line_length = 2 * max (rx, ry, rz)
		
	def lies_inside (self, p):
		rotated_p = p.rotate(-1 * self.alpha, -1 * self.beta, -1 * self.gamma)
		x = rotated_p.x - self.middle.x
		y = rotated_p.y - self.middle.y
		z = rotated_p.z - self.middle.z
		return (x * x / self.rx / self.rx + y * y / self.ry / self.ry + z * z / self.rz / self.rz ) <= 1
		
	def get_bounding_box (self):
		r = n3Rectangle(n3Point(self.middle.x - self.rx, self.middle.y - self.ry, self.middle.z - self.rz), \
		n3Point(self.middle.x + self.rx, self.middle.y + self.ry, self.middle.z + self.rz))
		r.rotate_by_self(self.alpha, self.beta, self.gamma)
		return r.get_bounding_box()
		
	def rotate_by_root (self, alpha, beta, gamma):
		self.middle = self.middle.rotate(alpha, beta, gamma)
		self.alpha += alpha
		self.beta += beta
		self.gamma += gamma
	
	def rotate_by_self (self, alpha, beta, gamma):
		self.alpha += alpha
		self.beta += beta
		self.gamma += gamma
	
	def translate (self, p):
		self.middle += p
	

class n3Pyramid(object):
	#TODO: implement
	pass


class n3Rectangle(object):
    
	def __init__ (self, points): 
		""" Simple constructor with following ordering of the points.
		"""		
		"""Orders points according to: Most left Bottom Point is first, then
			left upper, then right upper and finally right bottom. 	The 
			points are determined via the distance to coordinate origin
			and a point on the x axis lying to the right of the rectangle.
		"""
		self.points = points
		x = 0
		y = 0
		z = 0
		for p in points:
			x = min(x, p.x) 
			y = min(y, p.y)
			z = min(z, p.z)
		translate_vector = n3Point(x, y, z)
		ps = map(lambda n: n - translate_vector, points)
		max_y = 0
		for i in ps:
			if i.y > max_y:
				max_y = i.y
		max_x = 0
		for i in ps:
			if i.x > max_x:
				max_x = i.x
		d = []
		for i in ps:
			d.append(i.distance_to (n3Point(0,0,0)))
		d.sort(lambda x,y: -1 if x < y else 1)			
		if d[0]==d[1]:
			#idiotic special case, when rectangle was axis parallel and got rotated by 45 degrees
			self.left_front_lower, self.right_back_top = self.init_get_extreme_points(ps, n3Point(0, 0.1, 0))
			self.left_back_lower, self.right_front_top = self.init_get_extreme_points(ps, n3Point(0, max_y + 0.1, 0))
			self.right_back_lower, self.left_front_top = self.init_get_extreme_points(ps, n3Point(max_x, max_y + 0.1, 0))
			self.right_front_lower, self.left_back_top = self.init_get_extreme_points(ps, n3Point(max_x, 0.1, 0))
		else:
			self.left_front_lower, self.right_back_top = self.init_get_extreme_points(ps, n3Point(0, 0, 0))
			self.left_back_lower, self.right_front_top = self.init_get_extreme_points(ps, n3Point(0, max_y, 0))
			self.right_back_lower, self.left_front_top = self.init_get_extreme_points(ps, n3Point(max_x, max_y, 0))
			self.right_front_lower, self.left_back_top = self.init_get_extreme_points(ps, n3Point(max_x, 0, 0))
		self.left_front_lower = self.left_front_lower + translate_vector
		self.left_back_lower = self.left_back_lower + translate_vector
		self.right_back_lower = self.right_back_lower + translate_vector
		self.right_front_lower = self.right_front_lower + translate_vector
		self.left_front_top = self.left_front_top + translate_vector
		self.left_back_top = self.left_back_top + translate_vector
		self.right_back_top = self.right_back_top + translate_vector
		self.right_front_top = self.right_front_top + translate_vector
		self.compute_hesse_params()			
		self.inner_radius = ((self.right_front_lower - self.left_front_lower)/2).length
		self.outer_radius = (self.middle - self.left_front_lower).length
		self.diagonal_length = (self.left_front_lower - self.right_back_top).length
		self.longest_in_area_line_length = self.diagonal_length
		self.volume = (self.right_front_lower - self.left_front_lower).length * (self.left_back_lower - self.left_front_lower).length * (self.left_front_top - self.left_front_lower).length
			
	def compute_hesse_params (self):
		self.middle = self.left_front_lower + (self.right_back_top - self.left_front_lower) / 2
		self.norm_vector_side_bottom = (self.right_front_lower - self.left_front_lower).vector_product(self.left_back_lower - self.left_front_lower)
		self.norm_vector_side_bottom = self.norm_vector_side_bottom / self.norm_vector_side_bottom.length
		self.norm_vector_side_left = (self.left_front_lower - self.left_front_top).vector_product(self.left_back_lower - self.left_front_lower)		
		self.norm_vector_side_left = self.norm_vector_side_left / self.norm_vector_side_left.length
		self.norm_vector_side_front = (self.left_front_lower - self.right_front_lower).vector_product(self.left_front_top - self.left_front_lower)				
		self.norm_vector_side_front = self.norm_vector_side_front / self.norm_vector_side_front.length
		self.norm_vector_side_top = self.norm_vector_side_bottom * -1
		self.norm_vector_side_right = self.norm_vector_side_left * -1
		self.norm_vector_side_back = self.norm_vector_side_front * -1
		self.hesse_bottom_d = self.norm_vector_side_bottom.scalar_product(self.left_front_lower)
		self.hesse_left_d = self.norm_vector_side_left.scalar_product(self.left_front_lower)
		self.hesse_front_d = self.norm_vector_side_front.scalar_product(self.left_front_lower)
		self.hesse_back_d = self.norm_vector_side_back.scalar_product(self.right_back_top)
		self.hesse_right_d = self.norm_vector_side_right.scalar_product(self.right_back_top)
		self.hesse_top_d = self.norm_vector_side_top.scalar_product(self.right_back_top)		
		
	def init_get_extreme_points(self, liste, p):
         l = []
         for i in liste:
             l.append((i.distance_to (p), i))
         l.sort(lambda x,y: -1 if x < y else 1)
         return (l[0][1], l[7][1])
  
	def get_ordered_points(self):
		return [self.left_front_lower, self.left_back_lower, self.right_back_lower, self.right_front_lower, self.left_front_top, self.left_back_top, self.right_back_top, self.right_front_top]
		
	def lies_inside (self, p):
		d = self.middle.distance_to(p)
		if d > self.outer_radius:
			return False
		elif d <= self.inner_radius:
			return True
		else :
			return (p.scalar_product(self.norm_vector_side_back) - self.hesse_back_d) >= 0 and (p.scalar_product(self.norm_vector_side_front) - self.hesse_front_d) >= 0 and \
				(p.scalar_product(self.norm_vector_side_left) - self.hesse_left_d) >= 0 and (p.scalar_product(self.norm_vector_side_right) - self.hesse_right_d) >= 0 and \
				(p.scalar_product(self.norm_vector_side_bottom) - self.hesse_bottom_d) >= 0 and (p.scalar_product(self.norm_vector_side_top) - self.hesse_top_d) >= 0
				
	def rotate_by_root (self, alpha, beta, gamma):
		self.left_front_lower = self.left_front_lower.rotate_by_degree(alpha, beta, gamma)
		self.left_back_lower = self.left_back_lower.rotate_by_degree(alpha, beta, gamma)
		self.right_back_lower = self.right_back_lower.rotate_by_degree(alpha, beta, gamma)
		self.right_front_lower = self.right_front_lower.rotate_by_degree(alpha, beta, gamma)
		self.left_front_top = self.left_front_top.rotate_by_degree(alpha, beta, gamma)
		self.left_back_top = self.left_back_top.rotate_by_degree(alpha, beta, gamma)
		self.right_back_top = self.right_back_top.rotate_by_degree(alpha, beta, gamma)
		self.right_front_top = self.right_front_top.rotate_by_degree(alpha, beta, gamma)
		for p in xrange(len(self.points)):
			self.points[p] = self.points[p].rotate_by_degree(alpha, beta, gamma)
		self.compute_hesse_params()
		
	def translate (self, p):
		self.left_front_lower +=  p
		self.left_back_lower += p
		self.right_back_lower += p
		self.right_front_lower += p
		self.left_front_top += p
		self.left_back_top += p
		self.right_back_top += p
		self.right_front_top += p
		for i in xrange(len(self.points)):
			self.points[i] += p
		self.compute_hesse_params()
			
	def rotate_by_self (self, alpha, beta, gamma):
		middle = self.middle
		self.translate(middle * -1)
		self.rotate_by_root(alpha, beta, gamma)
		self.translate(middle)
		
	def get_bounding_box (self):
		min_x, min_y, min_z, max_x, max_y, max_z = self.points[0].x, self.points[0].y, self.points[0].z, self.points[0].x, self.points[0].y, self.points[0].z
		for p in self.points[1:]:
			min_x = min( min_x, p.x)
			min_y = min( min_y, p.y)
			min_z = min( min_z, p.z)
			max_x = max( max_x, p.x)				
			max_y = max( max_y, p.y)
			max_z = max( max_z, p.z)
		return n3AxisParallelRectangle(n3Point(min_x, min_y, min_z), n3Point(max_x, max_y, max_z))
  
		
class n3AxisParallelRectangle (n3Rectangle):
	""" The lines of the rectangle are parallel to the axis.
	"""	
	
	def __init__(self, p1, p7):
		p2 = n3Point(p1.x, p7.y, p1.z)
		p3 = n3Point(p7.x, p7.y, p1.z)
		p4 = n3Point(p7.x, p1.y, p1.z)
		p5 = n3Point(p1.x, p1.y, p7.z)
		p6 = n3Point(p1.x, p7.y, p7.z)
		p8 = n3Point(p7.x, p1.y, p7.z)
		
		self.points = [p1, p2, p3, p4, p5, p6, p7, p8]				
		self.left_front_lower = p1
		self.left_back_lower = p2
		self.right_back_lower = p3
		self.right_front_lower = p4
		self.left_front_top = p5
		self.left_back_top = p6
		self.right_back_top = p7
		self.right_front_top = p8
		self.middle = self.left_front_lower + (self.right_back_top - self.left_front_lower) / 2
		self.diagonal_length = (self.left_front_lower - self.right_back_top).length
		self.longest_in_area_line_length = self.diagonal_length
	
	def lies_inside (self, p1):
		return self.left_front_lower.x <= p1.x <= self.right_front_lower.x and  \
			self.left_front_lower.y <= p1.y <= self.left_back_lower.y and \
			self.left_front_lower.z <= p1.z <= self.left_front_top.z
			
	def get_bounding_box (self):
		return deepcopy(self)


class n2Point(object):
	"""Actually this is a location vector. And therefore some vector operations are supported.
	"""		

	def __init__ (self, x, y):
		""" Simple constructor.
		"""
		self.x = x
		self.y = y
		self.length = sqrt (x * x + y * y)
		
		
	def __add__ (self, p):
		""" Simple vector addition. Returns a new vector lengthend by the vector.
		"""
		return self.addition_with_vector(p)

	def addition_with_vector (self, p):
		""" Simple vector addition. Returns a new vector lengthend by the vector.
		"""
		return self.__class__(self.x + p.x, self.y + p.y)
	
	def __sub__(self,p):
		""" Simple vector subraction. Returns a new vector shortened by the vector.
		"""
		return self.subtraction_from_vector(p)
	
	def subtraction_from_vector (self, p):
		""" Simple vector subraction. Returns a new vector shortened by the vector.
		"""
		return self.__class__(self.x - p.x, self.y - p.y)
		
	def scalar_product (self, p):
		""" Vector with vector multiplication. Returns the scalar product.
		"""
		return self.x * p.x + self.y * p.y

	def __mul__ (self, scalar):
		""" Simple vector multiplication. Returns a new vector lengthend by the scalar.
		"""
		return self.__class__(self.x * scalar, self.y * scalar)

	def multiplikation_with_scalar (self, scalar):
		""" Simple vector multiplication. Returns a new vector lengthend by the scalar.
		"""
		return self.__class__(self.x * scalar, self.y * scalar)

	def __div__ (self, scalar):
		""" Simple vector division. Returns a new vector shortened by the scalar.
		"""	
		return  self.division_from_scalar(scalar)
		
	def __truediv__ (self, scalar):
		""" Simple vector division. Returns a new vector shortened by the scalar.
		"""
		return  self.division_from_scalar(scalar)

	def division_from_scalar (self, scalar):
		""" Simple vector division. Returns a new vector shortened by the scalar.
		"""
		return self.__class__(self.x / float(scalar), self.y / float(scalar))
		
	def __str__(self):
		""" Returns a readable string.
		"""
		return str("(" + str(self.x) + "," + str(self.y) + ")")

	def __repr__(self):
		""" Returns a readable string.
		"""
		return str("(" + str(self.x) + "," + str(self.y) + ")")

	def distance_to (self, p):
		""" Returns the distance between this and the given point by computing
			the length of the resulting vector.			
		"""
		return (self - p).length

	def get_position(self):
		""" Returns its position which is basically itself.
		"""
		return self
		
	def get_right_normal (self):
		return self.__class__(self.y, self.x * -1)
		
	def rotate_by_degree (self, alpha):
		""" Im Uhrzeigersinn
		"""
		alpha *= -1
		alpha = radians (alpha)
		x = round( self.x * cos (alpha) - self.y * sin (alpha), 4)
		y = round( self.x * sin (alpha) + self.y * cos (alpha), 4)
		return self.__class__(x, y)		


class n2Line (object):
	""" This class defines a line. It is defined by 2 points: head and tail.
		For every push head and tail are newly set and the old ones are 
		forgotten. Finally it implements also a function computing the distance
		to a given point.
	"""		
	
	def __init__(self, *points):
		""" Simple constructor. Pushes the line by iterating over all points
		provided. So the last point is the head and the one before last is 
		the tail.
		"""		
		self.head = points[0]
		for p in points[1:]:
			self.push(p)
		
	def push (self, p): 
		""" Pushes the line by adding the point p. The new head is point p 
			and the new tail is the old head.
		"""
		self.tail = self.head
		self.head = p
		self.direction_vector = self.head - self.tail
		self.direction_vector_norm = self.direction_vector / self.direction_vector.length
		self.middle = self.tail + (self.head - self.tail) / 2

	def push_delta (self, p):
		""" Pushes the line by adding the coordinates form the point p to the
			actual head as the new head. Finally setting the tail to the old
			head.
		"""
		self.tail = self.head
		self.head = self.head + p
		self.direction_vector = self.head - self.tail
		self.direction_vector_norm = self.direction_vector / self.direction_vector.length
		self.middle = self.tail + (self.head - self.tail) / 2

	def get_middle (self):
		""" Returns the point representing the middle.
		"""
		return self.middle 

	def get_head (self):
		""" Returns the point representing the head.
		"""
		return self.head
		
	def get_tail (self):
		""" Returns the point representing the tail.
		"""
		return self.tail

	def compute_distance (self, p):		
		""" Computes the distance between a point p and this line. It is 
			basically the minimum distance over all the points from the 
			line and the point p. This is done by computing the perpendicular
			and then using euclid.		
		"""		
		v = (self.tail - p) * -1 
		r = self.direction_vector.scalar_product(v) / self.direction_vector.length
		if r < 0:
			return self.tail.distance_to(p)
		elif r > self.direction_vector.length:
			return self.head.distance_to(p)
		else:
			return p.distance_to(self.direction_vector_norm * r + self.tail)

	def get_length (self):
		""" Returns the length of this line
		"""
		return (self.head - self.tail).length


class n2Sphere (object):
	
	def __init__ (self, point, radius):
		self.middle = point
		self.radius = radius
		self.volume = pi * radius * radius
		self.longest_in_area_line_length = 2 * radius
		
	def lies_in (self, p):
		return self.middle.distance_to(p) <= self.radius
		
	def get_bounding_box (self):
		return n2AxisParallelRectangle(n2Point(self.middle.x - self.radius, self.middle.y - self.radius), \
		n2Point(self.middle.x + self.radius, self.middle.y + self.radius))
		
	def translate (self, p):
		self.middle += p
		
	def rotate_by_root (self, alpha):
		self.middle = self.middle.rotate_by_degree(alpha)
    
	def rotate_by_self (self, alpha):
		print "Rotating a sphere with the axis aligned on it's center won't do anything."
    
				
class n2Ellipsoid (object):
	
	def __init__ (self, middle, rx, ry):
		self.middle = middle
		self.rx = rx
		self.ry = ry
		self.alpha = 0
		self.volume = pi * rx * ry
		self.longest_in_area_line_length = 2 * max (rx, ry)		
		
	def lies_inside (self, p):
		rotated_p = p.rotate(-1 * self.alpha, -1 * self.beta)
		x = rotated_p.x - self.middle.x
		y = rotated_p.y - self.middle.y
		return (x * x / self.rx / self.rx + y * y / self.ry / self.ry) <= 1
		
	def get_bounding_box (self):
		r =  n2Rectangle(n2Point(self.middle.x - self.rx, self.middle.y - self.ry), \
		n2Point(self.middle.x + self.rx, self.middle.y + self.ry))
		return r.get_
		
	def rotate_by_root (self, alpha):
		self.middle = self.middle.rotate(alpha)
		self.alpha += alpha
	
	def rotate_by_self (self, alpha):
		self.alpha += alpha
	
	def translate (self, p):
		self.middle += p
	

class n2Rectangle (object):
	""" This class implements a rectangle and is implemented based on n2Point"""	
	
	def __init__ (self, p1, p2, p3, p4): 
		""" Simple constructor with following ordering of the points.
		"""		
			
		points = [p1, p2, p3, p4]
		self.points = points
		
		max_x = 0
		for i in points:
			if i.y > max_x:
				max_x = i.y

		d = []
		for i in [p1, p2, p3, p4]:
			d.append(i.distance_to (n2Point(0,0)))
		d.sort(lambda x,y: -1 if x < y else 1)			
		if d[0]==d[1]:
			#idiotic special case, when rectangle was axis parallel and got rotated by 45 degrees
			self.left_front, self.right_back = self.init_get_extreme_points(points, n3Point(0, 0.1, 0))
			self.right_front, self.left_back = self.init_get_extreme_points(points, n3Point(max_x, 0.1, 0))
		else:
			self.left_front, self.right_back = self.init_get_extreme_points(points, n3Point(0, 0, 0))
			self.right_front, self.left_back = self.init_get_extreme_points(points, n3Point(max_x, 0, 0))
		self.inner_radius = ((self.right_front - self.left_front)/2).length
		self.middle = self.left_front + (self.right_back - self.left_front) / 2
		self.outer_radius = (self.middle - self.left_front).length
		self.compute_hesse_params()
		self.volume = (self.left_back - self.left_front).length * (self.right_front - self.left_front).length
		
	def compute_hesse_params (self):
		self.middle = self.left_front + (self.right_back - self.left_front) / 2
		self.norm_vector_side_left = (self.left_back - self.left_front).get_right_normal()
		self.norm_vector_side_left = self.norm_vector_side_left / self.norm_vector_side_left.length
		self.norm_vector_side_right = self.norm_vector_side_left * -1
		self.norm_vector_side_front = (self.left_front - self.right_front).get_right_normal()
		self.norm_vector_side_front = self.norm_vector_side_front / self.norm_vector_side_front.length
		self.norm_vector_side_back = self.norm_vector_side_front * -1		
		self.hesse_left_d = self.norm_vector_side_left.scalar_product(self.left_front)
		self.hesse_front_d = self.norm_vector_side_front.scalar_product(self.left_front)
		self.hesse_back_d = self.norm_vector_side_back.scalar_product(self.right_back)
		self.hesse_right_d = self.norm_vector_side_right.scalar_product(self.right_back)
		
	def init_get_extreme_points(self, liste, p):
         l = []
         for i in liste:
             l.append((i.distance_to (p), i))
         l.sort(lambda x,y: -1 if x < y else 1)
         return (l[0][1], l[3][1])
				
	def get_ordered_points (self): 
		""" Simple getter method. Returns the points of the rectangle in a list.
			The first element is the point in the left bottom, the second
			is in the left upper corner, the next in the right upper corner
			and finally the last is in the right bottom.		
		"""
		return [self.left_front, self.left_back, self.right_back, self.right_front,]
	
	def get_middle (self):
		""" Returns the point in the middle of the rectangle
		"""
		return self.left_front + (self.right_back - self.left_front) / 2

	def lies_inside (self, p):
		d = self.middle.distance_to(p)
		if d > self.outer_radius:
			return False
		elif d <= self.inner_radius:
			return True
		else :
			return (p.scalar_product(self.norm_vector_side_back) - self.hesse_back_d) >= 0 and (p.scalar_product(self.norm_vector_side_front) - self.hesse_front_d) >= 0 and \
				(p.scalar_product(self.norm_vector_side_left) - self.hesse_left_d) >= 0 and (p.scalar_product(self.norm_vector_side_right) - self.hesse_right_d) >= 0
				
	def rotate_by_root (self, alpha):
		self.left_front = self.left_front.translate_by_degree(alpha)
		self.left_back = self.left_back.translate_by_degree(alpha)
		self.right_back = self.right_back.translate_by_degree(alpha)
		self.right_front = self.right_front.translate_by_degree(alpha)		
		self.compute_hesse_params()
		
	def translate (self, p):
		self.left_front = self.left_front + p
		self.left_back = self.left_back + p
		self.right_back = self.right_back + p
		self.right_front = self.right_front + p		
		self.compute_hesse_params()
			
	def rotate_by_self (self, alpha):
		middle = self.middle
		self.translate(middle * -1)
		self.rotate_by_root(alpha)
		self.translate(middle)
  
	def get_bounding_box (self):
		min_x, min_y, max_x, max_y = self.points[0].x, self.points[0].y, self.points[0].x, self.points[0].y 
		for p in self.points[1:]:
			min_x = min( min_x, p.x)
			min_y = min( min_y, p.y)			
			max_x = max( max_x, p.x)				
			max_y = max( max_y, p.y)			
		return n2AxisParallelRectangle(n3Point(min_x, min_y), n3Point(max_x, max_y))

class n2AxisParallelRectangle (n2Rectangle):
	""" The lines of the rectangle are parallel to the axis.
	"""	
	
	def __init__(self, p1, p2):
		p3 = n2Point(p1.x, p2.y)
		p4 = n2Point(p2.x, p1.y)
		points = [p1, p2, p3, p4]
		self.points = points
		
		max_x = 0
		for i in points:
			if i.y > max_x:
				max_x = i.y
			
		self.left_front, self.right_back = self.init_get_extreme_points(points, n3Point(0, 0, 0))
		self.right_front, self.left_back = self.init_get_extreme_points(points, n3Point(max_x, 0, 0))
		self.middle = self.left_front + (self.right_back - self.left_front) / 2
	
	def lies_inside (self, p1):
		return self.lbPoint.x <= p1.x <= self.ruPoint.x and self.lbPoint.y <= p1.y <= self.ruPoint.y
		
	def get_bounding_box (self):
		return deepcopy(self)
		

class n3GeoTree(object):
	
	def __init__ (self, limit, cube):
		self.limit = limit
		self.container = []
		self.cube = cube
		self.is_leaf = True
		self.middle_x = cube.get_middle().x
		self.middle_y = cube.get_middle().y
		self.middle_z = cube.get_middle().z

	def add_element (self, element):				
		trie = self
		while not trie.is_leaf:
			if element.x <= trie.middle_x:
				#left
				if element.y <= trie.middle_y:
					#front
					if element.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if element.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if element.y <= trie.middle_y:
					#front
					if element.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if element.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		trie.container.append(element)
		if len(trie.container) > trie.limit:
				#split				
				p = trie.cube.get_points()				
				#bottom
				p.insert(1, p[0] + (p[1] - p[0])/2)
				p.insert(3, p[2] + (p[3] - p[2])/2)
				p.insert(5, p[4] + (p[5] - p[4])/2)
				p.insert(7, p[6] + (p[0] - p[6])/2)
				p.insert(8, p[0] + (p[4] - p[0])/2)
				#top				
				p.insert(10, p[9] + (p[10] - p[9])/2)
				p.insert(12, p[11] + (p[12] - p[11])/2)
				p.insert(14, p[13] + (p[14] - p[13])/2)
				p.append(p[15] + (p[9] - p[15])/2)
				p.append(p[9] + (p[13] - p[9])/2)				
				#middle				
				p.insert(9, p[9] + (p[0] - p[9])/2)
				p.insert(10, p[11] + (p[1] - p[11])/2)
				p.insert(11, p[13] + (p[2] - p[13])/2)
				p.insert(12, p[15] + (p[3] - p[15])/2)
				p.insert(13, p[17] + (p[4] - p[17])/2)
				p.insert(14, p[19] + (p[5] - p[19])/2)
				p.insert(15, p[21] + (p[6] - p[21])/2)
				p.insert(16, p[23] + (p[7] - p[23])/2)
				p.insert(17, p[25] + (p[8] - p[25])/2)
				
				trie.cube_left_front_lower = n3GeoTree(self.limit, n3AxisParallelRectangle(p[0], p[17]))
				trie.cube_left_back_lower  = n3GeoTree(self.limit, n3AxisParallelRectangle(p[1], p[12]))
				trie.cube_right_back_lower = n3GeoTree(self.limit, n3AxisParallelRectangle(p[8], p[13]))
				trie.cube_right_front_lower= n3GeoTree(self.limit, n3AxisParallelRectangle(p[7], p[14]))
				trie.cube_left_front_top    = n3GeoTree(self.limit, n3AxisParallelRectangle(p[9], p[26]))
				trie.cube_left_back_top     = n3GeoTree(self.limit, n3AxisParallelRectangle(p[10], p[21]))
				trie.cube_right_back_top    = n3GeoTree(self.limit, n3AxisParallelRectangle(p[17], p[22]))
				trie.cube_right_front_top   = n3GeoTree(self.limit, n3AxisParallelRectangle(p[16], p[23]))				
				
				"""
				trie.cube_left_front_lower = n3GeoTree(self.limit, n3AxisParallelRectangle([p[0], p[1], p[8], p[7], p[9], p[10], p[17], p[16]]))
				trie.cube_left_back_lower  = n3GeoTree(self.limit, n3AxisParallelRectangle([p[1], p[2], p[3], p[8], p[10], p[11], p[12], p[17]]))
				trie.cube_right_back_lower = n3GeoTree(self.limit, n3AxisParallelRectangle([p[8], p[3], p[4], p[5], p[17], p[12], p[13], p[14]]))
				trie.cube_right_front_lower= n3GeoTree(self.limit, n3AxisParallelRectangle([p[7], p[8], p[5], p[6], p[16], p[17], p[14], p[15]]))
				trie.cube_left_front_top    = n3GeoTree(self.limit, n3AxisParallelRectangle([p[9], p[10], p[17], p[16], p[18], p[19], p[26], p[25]]))
				trie.cube_left_back_top     = n3GeoTree(self.limit, n3AxisParallelRectangle([p[10], p[11], p[12], p[17], p[19], p[20], p[21], p[26]]))
				trie.cube_right_back_top    = n3GeoTree(self.limit, n3AxisParallelRectangle([p[17], p[12], p[13], p[14], p[26], p[21], p[22], p[23]]))
				trie.cube_right_front_top   = n3GeoTree(self.limit, n3AxisParallelRectangle([p[16], p[17], p[14], p[15], p[25], p[26], p[23], p[24]]))								
				"""				
				
				# assume the data is nearly linearly distributed, so we don't have to check whether our newly created container gets full
				for el in trie.container:					
					if el.x <= trie.middle_x:
						#left
						if el.y <= trie.middle_y:
							#front
							if el.z <= trie.middle_z:
								#bottom
								trie.cube_left_front_lower.container.append(el)
							else:
								#top
								trie.cube_left_front_top.container.append(el)								
						else:
							#back
							if el.z <= trie.middle_z:
								#bottom
								trie.cube_left_back_lower.container.append(el)
							else:
								#top
								trie.cube_left_back_top.container.append(el)
					else:
						#right
						if el.y <= trie.middle_y:
							#front
							if el.z <= trie.middle_z:
								#bottom
								trie.cube_right_front_lower.container.append(el)
							else:
								#top
								trie.cube_right_front_top.container.append(el)								
						else:
							#back
							if el.z <= trie.middle_z:
								#bottom
								trie.cube_right_back_lower.container.append(el)
							else:
								#top
								trie.cube_right_back_top.container.append(el)
				trie.is_leaf = False
				trie.container = []	
			
	def get_sourrounding_points (self, coords):
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		#return deepcopy(trie.container)
		return trie.container

	def get_points_within_radius (self, coord, radius):
		result = []		
		
		coords = n3Point(coord.x - radius, coord.y, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		result.append(trie)
		
		coords = n3Point(coord.x + radius, coord.y, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y - radius, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y + radius, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y, coord.z - radius)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y, coord.z + radius)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_front_lower
					else:
						#top
						trie = trie.cube_left_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_left_back_lower
					else:
						#top
						trie = trie.cube_left_back_top
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_front_lower
					else:
						#top
						trie = trie.cube_right_front_top
						
				else:
					#back
					if coords.z <= trie.middle_z:
						#bottom
						trie = trie.cube_right_back_lower
					else:
						#top
						trie = trie.cube_right_back_top
		result.append(trie)
		
		resultset = set(result)
		result = []
		for t in resultset:
			result[len(result):] = t.container
		return result

class n2GeoTree (object):
	"""Index Structure for accessing spatial information in log(n) time.
		It uses 4 (and splitting) axis parallel rectangles as way to access
		spatial data. It works like a tree on those rectangles
	"""
	
	def __init__ (self, limit, rectangle):
		self.limit = limit
		self.container = []
		self.rectangle = rectangle
		self.is_leaf = True
		self.middle_x = rectangle.get_middle().x
		self.middle_y = rectangle.get_middle().y

	def add_element (self, element):				
		trie = self
		while not trie.is_leaf:
			if element.x <= trie.middle_x:
				#left
				if element.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if element.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		trie.container.append(element)
		if len(trie.container) > trie.limit:
				#split				
				p = trie.rectangle.get_points()				
				
				p.insert(1, p[0] + (p[1] - p[0])/2)
				p.insert(3, p[2] + (p[3] - p[2])/2)
				p.insert(5, p[4] + (p[5] - p[4])/2)
				p.insert(7, p[6] + (p[0] - p[6])/2)
				p.insert(8, p[0] + (p[4] - p[0])/2)

				trie.cube_left_front = n2GeoTree(self.limit, n2AxisParallelRectangle(p[0], p[8]))
				trie.cube_left_back  = n2GeoTree(self.limit, n2AxisParallelRectangle(p[1], p[3]))
				trie.cube_right_back = n2GeoTree(self.limit, n2AxisParallelRectangle(p[8], p[4]))
				trie.cube_right_front= n2GeoTree(self.limit, n2AxisParallelRectangle(p[7], p[5]))				
				
				# assume the data is nearly linearly distributed, so we don't have to check whether our newly created container gets full
				for el in trie.container:					
					if element.x <= trie.middle_x:
						#left
						if element.y <= trie.middle_y:
							#front					
							trie = trie.rect_left_front.container.append(el)
						else:
							#back
							trie = trie.rect_left_back.container.append(el)
					else:
						#right
						if element.y <= trie.middle_y:
							#front
							trie = trie.rect_right_front.container.append(el)
						else:
							#back
							trie = trie.rect_right_back.container.append(el)
				trie.is_leaf = False
				trie.container = []	
			
	def get_sourrounding_points (self, coords):
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		#return deepcopy(trie.container)
		return trie.container

	def get_points_within_radius (self, coord, radius):
		result = []		
		
		coords = n3Point(coord.x - radius, coord.y, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		result.append(trie)
		
		coords = n3Point(coord.x + radius, coord.y, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y - radius, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y + radius, coord.z)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y, coord.z - radius)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		result.append(trie)
		
		coords = n3Point(coord.x, coord.y, coord.z + radius)
		trie = self
		while not trie.is_leaf:
			if coords.x <= trie.middle_x:
				#left
				if coords.y <= trie.middle_y:
					#front					
					trie = trie.rect_left_front
				else:
					#back
					trie = trie.rect_left_back
			else:
				#right
				if coords.y <= trie.middle_y:
					#front
					trie = trie.rect_right_front						
				else:
					#back
					trie = trie.rect_right_back
		result.append(trie)
		
		resultset = set(result)
		result = []
		for t in resultset:
			result[len(result):] = t.container
		return result


class TwoDimEuclidMetric (object):
	""" Simple 2 dimensional euclidean metric to compute the distance between 
	to points
	"""
	
	def compute_distance (self, p1, p2):
		a = p1.x - p2.x
		b = p1.y - p2.y
		return sqrt (a*a + b*b)
		
		
class ThreeDimEuclidMetric (object):
	""" Simple 3 dimensional euclidean metric to compute the distance between 
	to points
	"""
	
	def compute_distance (self, p1, p2):
		a = p1.x - p2.x
		b = p1.y - p2.y
		c = p1.z - p2.z 		
		return sqrt (a*a + b*b + c*c)
		
	
class Noise ():
	""" "Static class" with can be used to noise data
	"""	
	
	@staticmethod
	def add_gauss_noise (mean, variance):
		""" Returning point is gaussion distributed i.r.t. mean and variance
			If result is below zero, zero will be returned instead
		"""
		return max (0, int (random.gauss(mean, variance)))
		
	@staticmethod
	def add_linear_noise(lower_limit, upper_limit):
		""" Returns a point linearly distributed between both limits
			If result is below zero, zero will be returned instead
		"""
		return max (0, int (random.randint(lower_limit, upper_limit)))

"""
	The following lines define functions, which can be used for time distributions
	for the generators. Every class (should) define a get_value method, which
	computes a y (function value) for the given x.
"""

class GaussAlikeFunction ():
	""" This is originally a gauss function. But the area sums not to 1 anymore.
		But is scaled so, that the mean has the function value "highest_point".
		Means a function with very high variance behaves like a constant 
		function: for every x the function value is highest_point
	"""
	
	def __init__ (self, mean, standard_deviation, highest_point, lower_bound, upper_bound ):
		self.mean = mean
		self.highest_point = highest_point
		self.variance = 	standard_deviation * standard_deviation		
		self.lower_bound = lower_bound
		self.upper_bound = upper_bound	
	
	def get_value (self, x):
		""" Returns the function value for a given x
		"""
		if x < self.lower_bound or x > self.upper_bound:
			return 0
		else:	
			# proudly programmed by Katja			
			#gauss_param = 1 / sqrt(pi * self.variance*self.variance * 2)  * pow (e, - 1 * (pow((x - self.mean) / sqrt(2) / self.variance, 2)))
			#gauss_param = 1 / pow (e, pow (x - self.mean, 2) / 2 / self.variance) / sqrt (2 * pi * self.variance)			
			return self.highest_point / pow (e, pow (x - self.mean, 2) / 2 / self.variance)


class FuzzyFunction ():
	""" This function looks like a pyramid. The mean is the top of the pyramid,
		with lower and upper bound defining the base.
	
	"""	
	
	def __init__ (self, mean, highest_point, lower_bound, upper_bound):
		self.mean = mean
		self.highest_point = highest_point
		self.lower_bound = lower_bound
		self.upper_bound = upper_bound	
	
	def get_value (self, x):
		""" Returns the function value for a given x
		"""
		if x < self.lower_bound or x > self.upper_bound:
			return 0
		elif x <= self.mean :
			result = self.highest_point * (x - self.lower_bound) / (self.mean - self.lower_bound)
		else :
			result = self.highest_point * (1 - float(x - self.mean) / (self.upper_bound - self.mean))
		return result


class QuadraticFunction ():
	""" It needs 3 Points and then will compute the underlying square function
		of it. 	
	"""
		
	def __init__ (self, p1, p2, p3, lower_bound, upper_bound):		
		""" The computation of the function is according to this website (20.02.2015)
			http://www.arndt-bruenner.de/mathe/10/parabeldurchdreipunkte.htm
		"""
		liste = [p1, p2, p3]		
		liste.sort()
		x1 = liste[0].x
		x2 = liste[1].x
		x3 = liste[2].x
		y1 = liste[0].y
		y2 = liste[1].y
		y3 = liste[2].y
		self.a = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)
		self.a /= (x1 - x2) * (x1 - x3) * (x3 - x2)
		self.b = x1 * x1 * (y2 - y3) + x2 * x2 * (y3 - y1)	 + x3 * x3 * (y1 - y2)
		self.b /= (x1 - x2) * (x1 - x3) * (x2 - x3)
		self.c = x1 * x1 * (x2 * y3 - x3 * y2) + x1 * (x3 * x3 * y2 - x2 * x2 * y3) + x2 * x3 * y1 * (x2 - x3)
		self.c /= (x1 - x2) * (x1 - x3) * (x2 - x3)
		self.lower_bound = lower_bound
		self.upper_bound = upper_bound
		
	def get_value (self, x):
		""" Returns the function value for a given x. It is always positiv.
			When below zero, zero is instead returned
		"""
		if x < self.lower_bound or x > self.upper_bound:
			return 0
		else:	
			return max(0, self.a * x * x + self.b * x + self.c)

class Scipy1DInterpolation ():
	""" This function uses a Scipy interpolation function to define our function.
		Attention: the range is strictly limited to the given points.
		If you receive:
		ValueError: A value in x_new is above/below the interpolation range.
		it is because of this.
	"""
	
	def __init__ (self, *points):
		from scipy.interpolate import interp1d
		x = []
		y = []
		for p in points:
			x.append(p.x)
			y.append(p.y)
		self.f = interp1d(x, y, kind='cubic')
		
	def get_value (self, x):
		""" Returns the function value for a given x
		"""
		return self.f(x) * 1
		
	
class Distances ():
	""" Class which holds the distance matrix and the according functions like
		adding new distances, getting distances and doing the statistics.
		The elements given in add and get have to have the property dist_index
		which donates the index in this matrix.
	"""
	
	def __init__ (self, init_container_size = 10, growth_factor = 1.8):
		""" Simple constructor, just defining variables.
		"""
		self.growth_factor = growth_factor		
		self.container = np.zeros((init_container_size,init_container_size), dtype = float)
			
	def add (self, element1, element2, distance_function):
		""" Adds a distance to the matrix. If the distance is already there,
			no computation is done. If not the distance function is called
			with the 2 given arguments and then written to the matrix.
			If the container can't hold the new elements, the container will
			be extended.
			Finally True will be returned if the connection was computed and
			added and False otherwise. (For example, if distance was already
			computed False will be returned)
		"""
		container_size = len(self.container)
		
		if element1.dist_index > container_size or element2.dist_index > container_size:
			# grow			
			new_container_size = len(self.container) * self.growth_factor
			new_container = np.zeros((new_container_size,new_container_size), dtype = float)			
			for i in xrange(container_size):
				new_container[i][:container_size] = self.container[i]			
			self.container = new_container
			return self.add(element1, element2, distance_function)
		elif	self.container[element1.dist_index-1][element2.dist_index-1] == 0:
			self.container[element1.dist_index-1][element2.dist_index-1] = distance_function(element1, element2)
			return True
		else:
			return False

	def get (self, element1, element2):
		""" Returns the distance between those 2 elements
		"""
		return self.container[element1.dist_index-1][element2.dist_index-1]

	def compute_statistics (self, upper_limit, partition_bins = 10):
		""" Computes the statistics for our model. Means make a bar diagramm
			for the distribution of the distances.
			The maximal distance for the partitioning is also needed.		
		"""
		result = np.zeros((partition_bins), dtype = int)
		for i in xrange(len(self.container)):
			for j in xrange(len(self.container[i])):
				if self.container[i][j] != 0:
					result[int(self.container[i][j]/upper_limit * partition_bins)] += 1
		return result
		
	def compute_arithmetic_mean_length (self):
		""" Computes the arithmetric mean of non-zero entries and returns it.
		"""
		result = 0
		number_of_elements = 0
		for i in xrange(len(self.container)):
			for j in xrange(len(self.container[i])):
				if self.container[i][j] != 0:
					number_of_elements += 1
					result += self.container[i][j]
		return result / float(number_of_elements)
		
	def get_min_length (self):	
		""" Returns the minimal length in this matrix. Zero entries don't count in. 
		"""
		smallest_element = 10000
		for i in xrange(len(self.container)):
			for j in xrange(len(self.container[i])):
				if self.container[i][j] != 0  and self.container[i][j] < smallest_element:
					smallest_element = self.container[i][j]
		return smallest_element
		
	def get_max_length (self):
		""" Returns the maximal length in this matrix
		"""
		return np.amax(self.container)

	def get_median_length (self):
		""" Returns the median length of the distance matrix. Zero entries 
			don't count in.
		"""
		result =[]
		for i in xrange(len(self.container)):
			for j in xrange(len(self.container[i])):
				if self.container[i][j] != 0:
					result.append(self.container[i][j])	
		return np.median(result)

	def get_matrix (self):
		""" Returns the distance matrix
		"""
		return deepcopy(self.container)
		
def create_chemical_gradients_field (area, sources, filepath, inhibitors=0):
		bounding_box = area.get_bounding_box()
		chemical_gradients = np.zeros(((bounding_box.left_front_lower - bounding_box.right_front_lower).length, \
		(bounding_box.left_front_lower - bounding_box.left_back_lower).length, \
		(bounding_box.left_front_lower - bounding_box.left_front_top).length, len(sources) ))
		diagonal_length = (bounding_box.right_back_top - bounding_box.left_front_lower).length
		print chemical_gradients.shape
		p = n3Point(0, 0, 0)
		lx = bounding_box.left_front_lower.x
		ly = bounding_box.left_front_lower.y
		lz = bounding_box.left_front_lower.z
		for source in xrange(len(sources)):
			print "Generator: "
			print source
			s = sources[source]
			for x in xrange(chemical_gradients.shape[0]):					
				for y in xrange(chemical_gradients.shape[1]):
					for z in xrange(chemical_gradients.shape[2]):
						p.x, p.y, p.z = x + lx, y + ly, z + lz
						if area.lies_inside(p):
							if inhibitors:									
								for inhibitor in inhibitors:
									chemical_gradients[x, y, z, source] -= max (0, ( 1 - inhibitor.middle.distance_to(p) / inhibitor.radius)  )
							chemical_gradients[x, y, z, source] += (1 - p.distance_to(s) / diagonal_length)
							chemical_gradients[x, y, z, source] = max (0, chemical_gradients[x, y, z, source])
		np.save(filepath, chemical_gradients)
		
def debug_chem_grad_field (field, gen, thr):
	for x in xrange(field.shape[0]):
		for y in xrange(field.shape[1]):
			for z in xrange(field.shape[2]):
				if field[x, y, z, gen] > thr:
					print "%i, %i, %i ist %.4f" %(x, y, z, field[x, y, z, gen])


class Simulation ():
	""" This is our "Central processing unit". It manages the whole simulation
		including distance matrix, geotree, neuron generators and neurons.
	"""
	
	def __init__ (self, neuron_generators):
		""" Normal constructor. Just defining some variables we need. And 
			constructing our super area, which contains all other areas.
		"""
		self.generators = neuron_generators		
		self.upper_limit_container_capacity = 0.55
		self.tree = 0
		self.simulation_area = 0
		self.verbose = 0

	def set_container_capacity_heuristic_limit (self, float_number):
		""" Sets the percentage when a container will be considered full and
			therefore no new neurons will be placed in there.
		"""
		self.upper_limit_container_capacity = float_number

	def set_verbosity (self, integer):
		self.verbose = integer

	def set_bounding_area (self, area):
		self.simulation_area = area
		self.bounding_box = area.get_bounding_box()
		self.max_distance = self.simulation_area.longest_in_area_line_length
		
	def set_chemical_gradient_field (self, gradient_field):
		self.chemical_gradients = gradient_field

	def use_geotree (self, boolean, capacity = 50, search_radius = 5):
		""" Tells the simulation whether a indexing structure should be used,
			to speed the computation up. With 300 neurons the factor is 3
			(theoretically and practically about 2).
			And with more neurons the factor increases.
		"""		
		if boolean:
			self.tree = n3GeoTree(capacity, self.simulation_area)
			self.tree_search_radius = search_radius
		else :
			self.tree = 0

	def simulation_step (self):
		""" This method simulates a simulation step. DO NOT CALL DIRECTLY.
		"""		
		self.simulation_step_counter += 1
		added_neurons = 0
		for generator in self.generators:		
			for area in generator.areas:
				if area.sim_active:
					neurons_to_add = generator.new_neurons(self.simulation_step_counter, area)					
					degree_of_capacity = area.occupied_volume / area.volume
					if degree_of_capacity > self.upper_limit_container_capacity:
						if self.verbose:
							print "Error: %s has reached quite its capacity. Skipping this container for the rest of the simulation." %(area.id)
						area.sim_active = False
						neurons_to_add = 0					
					elif degree_of_capacity > 0.35 and self.verbose:
						print "Warning: %s is going into capacity problems." %(area.id) 					
					added_neurons += neurons_to_add
					for i in xrange(neurons_to_add):				
						nneuron = generator.next_neuron(self.simulation_step_counter, area)
						while not self.is_free(nneuron, nneuron.cellbody_radius):
							nneuron = generator.next_neuron(self.simulation_step_counter, area)
						nneuron.axon = n3Line(nneuron.position)
						nneuron.active = True
						nneuron.dist_index = self.dist_counter
						area.ocu_space += 4 / 3 * pi * nneuron.cellbody_radius * nneuron.cellbody_radius * nneuron.cellbody_radius
						self.dist_counter += 1
						if self.tree:
							self.tree.add_element(nneuron)
						self.neurons.append(nneuron)
						self.neuron_path.append([nneuron.position])
		for neuron in self.neurons:			
			if neuron.active:
				# grow 
				if neuron.type == "long" :
					x = round(neuron.axon.head.x)
					y = round(neuron.axon.head.y)
					z = round(neuron.axon.head.z)
					target_area = neuron.target_area
					gradient_field = [self.chemical_gradients[x - 1][y - 1][z - 1][target_area], self.chemical_gradients[x - 1][y][z - 1][target_area], \
					self.chemical_gradients[x - 1][y + 1][z - 1][target_area], self.chemical_gradients[x][y + 1][z - 1][target_area], \
					self.chemical_gradients[x + 1][y + 1][z - 1][target_area], self.chemical_gradients[x + 1][y][z - 1][target_area], \
					self.chemical_gradients[x + 1][y - 1][z - 1][target_area], self.chemical_gradients[x][y - 1][z - 1][target_area], \
					self.chemical_gradients[x][y][z - 1][target_area], \
					self.chemical_gradients[x - 1][y - 1][z][target_area], self.chemical_gradients[x - 1][y][z][target_area], \
					self.chemical_gradients[x - 1][y + 1][z][target_area], self.chemical_gradients[x][y + 1][z][target_area], \
					self.chemical_gradients[x + 1][y + 1][z][target_area], self.chemical_gradients[x + 1][y][z][target_area], \
					self.chemical_gradients[x + 1][y - 1][z][target_area], self.chemical_gradients[x][y - 1][z][target_area], \
					self.chemical_gradients[x][y][z][target_area], \
					self.chemical_gradients[x - 1][y - 1][z + 1][target_area], self.chemical_gradients[x - 1][y][z + 1][target_area], \
					self.chemical_gradients[x - 1][y + 1][z + 1][target_area], self.chemical_gradients[x][y + 1][z + 1][target_area], \
					self.chemical_gradients[x + 1][y + 1][z + 1][target_area], self.chemical_gradients[x + 1][y][z + 1][target_area], \
					self.chemical_gradients[x + 1][y - 1][z + 1][target_area], self.chemical_gradients[x][y - 1][z + 1][target_area], \
					self.chemical_gradients[x][y][z + 1][target_area]]						
					neuron.axon.push(neuron.grow(gradient_field))
				elif neuron.type == "short":					
					neuron.axon.push_delta(neuron.grow())
				self.neuron_path[neuron.dist_index].append(neuron.axon.head)
				# make connections
				if neuron.can_put_connection():
					if self.tree :						
						center = neuron.axon.middle
						radius = neuron.axon.get_length + self.tree_search_radius
						neurons = self.tree.get_points_within_radius(center, radius)
					else :
						neurons = self.neurons
					for ntest in neurons:
						if ntest is not neuron and ntest.can_receive_connection(neuron) and \
							neuron.axon.compute_distance(ntest.position) < ntest.dendrite_radius and \
							self.dmatrix.add(neuron, ntest, self.generators[0].metric.compute_distance):
							neuron.put_connection(ntest)
							ntest.receive_connection(neuron, neuron.axon)
				neuron.active = neuron.can_put_connection() and self.simulation_area.lies_inside(neuron.axon.head)
		if self.verbose:
			print "Step %i: Added %i new neurons." %(self.simulation_step_counter, added_neurons)

	def simulate (self):
		""" Main method to start the simulation.
		"""
		self.neurons = []
		self.neuron_path = []	
		self.dmatrix = Distances(100)
		self.simulation_step_counter = 0
		self.dist_counter = 0

		id_generator = 1
		for generator in self.generators:	
			for area in generator.areas:			
				area.occupied_volume = 0.0
				area.sim_active = True
				if not hasattr(area, 'id'):
					area.id = 'Area%i' %(id_generator)
					id_generator += 1
		if not self.simulation_area:	
			min_x, min_y, min_z, max_x, max_y, max_z = 100000, 100000, 100000, 0, 0, 0		
			for generator in self.generators:				
				for area in generator.areas:
					bounding_box = area.get_bounding_box()
					min_x = min( min_x, bounding_box.left_front_lower.x)
					min_y = min( min_y, bounding_box.left_front_lower.y)
					min_z = min( min_y, bounding_box.left_front_lower.z)
					max_x = max( max_x, bounding_box.right_back_top.x)				
					max_y = max( max_y, bounding_box.right_back_top.y)
					max_z = max( max_z, bounding_box.right_back_top.z)
			self.simulation_area = n3AxisParallelRectangle(n3Point(min_x, min_y, min_z), n3Point(max_x, max_y, max_z))
			self.bounding_box = self.simulation_area
			self.max_distance = self.simulation_area.longest_in_area_line_length
		if self.tree :
			self.tree = n3GeoTree(self.tree.limit, self.super_area)
		if self.verbose:
			print "Adding Neurons and growing Axons"
		min_iterations = 0
		for gen in self.generators:
			min_iterations = max(min_iterations, gen.minimun_iterations)
		for i in xrange(min_iterations):
			self.simulation_step()
		if self.verbose:
			print "Finishing Growth"
		while not self.finished():			
			self.simulation_step()	
		
	def print_simulation_meta_data (self):
		""" Prints some statistics onto the console.
		"""
		print "\nNeurons: ---------------------------------------------------\n"
		print "%i Neurons were added and they made %i connections.\n" %(len(self.neurons), self.dmatrix.compute_statistics(self.max_distance, 1)[0])		
		print "Connections: -----------------------------------------------\n"		
		print "The arithmetric mean of all connection lengths is %f." %(self.dmatrix.compute_arithmetic_mean_length())
		print "The median connection has a length of %f."%(self.dmatrix.get_median_length())
		print "The shortest connection has a length of %f and the longest %f.\n" %(self.dmatrix.get_min_length(), self.dmatrix.get_max_length())		
		print "Containers: ------------------------------------------------\n"		
		for gen in self.generators:
			for area in gen.areas:
				if area.sim_active :
					print 'Container "%s" reached %i%s of its capacity during the simulation.' %(area.id, int (area.ocu_space / area.max_space * 100), "%")
				else :
					print 'Container "%s" hit the given limit and was, at some point, excluded from further neuronal placement. Please check the settings for this area.' %(area.id)
				
		print "A capacity of 70 % means the container is full.\n"
	
	def finished (self):
		""" Determines whether our simulation has finished.
		"""
		for neuron in self.neurons:
			if neuron.active:
				return False
		return True
				
	def is_free(self, point, radius):
		""" Computes whether the given area is already (partly) occupied.
		"""
		if self.tree:
			neurons = self.tree.get_elements_within_radius(point, radius)
		else :
			neurons = self.neurons
		for neuron in neurons:			
			if self.generators[0].metric.compute_distance(point.position, neuron.position) < (radius + neuron.cellbody_radius):
				return False		
		return True

	def get_neurons (self):
		""" Returns the list containg all neurons. Please use with care.
		"""
		return self.neurons

	def get_statistics (self, partition_bins = 10):
		""" Computes the statistics for our model. Means make a bar diagramm
			for the distribution of the distances.
		"""		
		return self.dmatrix.compute_statistics(self.max_distance, partition_bins)

	def get_distance_matrix (self):
		""" Returns the distance matrix for this model. The indexes are in
			order in which neurons was placed. Means: First added neuron is 
			index 0, second added neuron is index 1 and so on.
		"""
		return self.dmatrix

	def saveModelNeuronsAsFile (self, path):
		""" Saves our model (neurons) as coordinates in
			a text file
		"""
		with open(path,'w') as f: 		
			for n in self.neurons:
				s = n.get_position().__str__().strip('()') + "\n" 
				f.write(s)

	def saveModelNeuronsWithAxonsAsFile (self, path):
		""" Saves our model (neurons and respecting axons) as coordinates in
			a text file
		"""
		with open(path,'w') as f:
			for n in self.neurons:
				s = n.get_position().__str__().strip('()') + ", " + n.axon.get_head().__str__().strip('()') +  "\n"
				f.write(s)


class ShortDistanceNeuron (object):
		
	def __init__ (self, position, bounding_area, area_type, rotation_angle):
		""" Trivial constructor
		"""
		self.position = position
		self.type = "short"
		self.cellbody_radius = 0.2
		self.dendrite_radius = 2
		self.bounding_area = bounding_area
		self.rotation_angle = rotation_angle
		self.area_type = area_type
		self.incoming_connections = 0
		self.outgoing_connections = 0
		self.direction_x = random.random() - 0.5
		self.direction_y = random.random() - 0.5
		self.direction_z = max(-0.5, min (0.5, random.normalvariate(0, 0.2)))
		self.direction = n3Point(self.direction_x, self.direction_y, self.direction_z)
		self.direction = self.direction / self.direction.length
		self.direction = self.direction.rotate_by_degree(rotation_angle[0], rotation_angle[1], rotation_angle[2])
		self.growth_speed = 0.3
		self.axon_flexibility = 0.02
	

	def grow (self, *i):
		""" Is called everytime this neuron shall grow. The return value
			is a direction vector.
		"""
		self.direction.x = random.normalvariate(self.direction.x, self.axon_flexibility) # flexibility of axon
		self.direction.y = random.normalvariate(self.direction.y, self.axon_flexibility)
		self.direction.z = random.normalvariate(self.direction.z, self.axon_flexibility)
		self.direction = self.direction / self.direction.length * self.growth_speed # growthspeed		
		return self.direction
		
	def can_put_connection(self):
		""" Is called, when the axon from this neuron found an neuron and 
			wants to establish a connection OR whether the axon to this 
			neuron can still grow.			
		"""
		return self.outgoing_connections < 10 and self.bounding_area.lies_inside(self.axon.head)
	
	def put_connection(self, target_neuron):
		""" Is called, when the outgoing connection is made.
		"""
		self.outgoing_connections += 1
		
	def can_receive_connection(self):
		""" Is called, when another neuron tries to make a connection.
			The return value states whether this is okay or not
		"""
		return self.incoming_connections < 10
		
	def receive_connection(self, origin_neuron, axon):
		""" Is called, when another neuron makes a connection.
		"""
		self.incoming_connections += 1
	
	
class LongDistanceNeuron (object):
		
	def __init__ (self, position, target_area):
		""" Trivial constructor
		"""
						
		self.type = "long"
		self.position = position
		self.cellbody_radius = 0.2
		self.dendrite_radius = 3		
		self.bounding_area = 0
		self.incoming_connections = 0
		self.outgoing_connections = 0
		self.target_area = target_area		
		self.grow_speed_constant = 7
		self.axon_flexibility = 0.02
		self.growth_speed = 0.3
		self.area_type = 0
		
	def grow (self, gradient_field):
		""" modified glücksrad auswahl
		"""
		inverse_distance_to_area = max(gradient_field)
		direction = n3Point(0, 0, 0)
		if inverse_distance_to_area > 0.85 :
			# "white matter"	
			liste = map(lambda x: (max(0, x-0.5))**3, gradient_field)
			summe = random.random() * sum(liste)
			growth_speed = (1.1 - inverse_distance_to_area) * self.grow_speed_constant
			i = 0
			while summe > 0:
				summe += liste[i]
				i += 1
			z = i / 9 -1
			i = i - z * 9		
			if i == 0 or i == 1 or i == 2:
				x = -1
			elif i == 3 or i == 7 or i == 8:
				x = 0
			else:
				x = 1
			if i == 0 or i ==6 or i == 7:
				y = -1 
			elif i == 1 or i == 5 or i == 8:
				y = 0
			else:
				y = 1
			direction.x = x * growth_speed + self.axon.head.x
			direction.y = y * growth_speed + self.axon.head.y
			direction.z = z * growth_speed + self.axon.head.z
		else:
			# "grey matter"
			direction.x = random.normalvariate(self.axon.direction_vector_norm.x , self.axon_flexibility) # flexibility of axon
			direction.y = random.normalvariate(self.axon.direction_vector_norm.y, self.axon_flexibility)
			direction.z = random.normalvariate(self.axon.direction_vector_norm.z, self.axon_flexibility)
			direction = direction / direction.length * self.growth_speed 
		return direction
			
		
	def can_put_connection(self, target_neuron):
		""" Is called, when the axon from this neuron found an neuron and 
			wants to establish a connection OR whether the axon to this 
			neuron can still grow.			
		"""		
		return self.outgoing_connections < 20 and target_neuron.area_type == self.target_area
	
	def put_connection(self, target_neuron):
		""" Is called, when the outgoing connection is made.
		"""
		self.outgoing_connections += 1
		
	def can_receive_connection(self):
		""" Is called, when another neuron tries to make a connection.
			The return value states whether this is okay or not
		"""
		return self.incoming_connections < 20
		
	def receive_connection(self, origin_neuron, axon):
		""" Is called, when another neuron makes a connection.
		"""
		self.incoming_connections += 1
		
	
class LayersNeuronGenerator (object):
	""" Very trivial generator. Every method implemented is used by the simulation.
		Every generator class must either directly or indirectly inherit from this.
	"""
	
	def __init__(self, area, area_type, area_target, rotation_angle):
		""" Trivial constructor. Metric and areas should generally exist.
		"""
		self.metric = ThreeDimEuclidMetric()
		self.minimun_iterations = 2
		self.area_type = area_type
		self.area_target = area_target
		self.chemical_gradient_source = area.middle
		points = area.get_ordered_points()
		self.bounding_box = area
		self.rotation_angle = rotation_angle
		bottom = points[:4]
		top = points [4:]
		layer6 = n3Rectangle([bottom[0], bottom[1], bottom[2], bottom[3], \
			bottom[0] + (top[0]-bottom[0])/6.0, bottom[1] + (top[1]-bottom[1])/6.0, \
			bottom[2] + (top[2]-bottom[2])/6.0, bottom[3] + (top[3]-bottom[3])/6.0])
		layer6.id = "Layer6"
		layer5 = n3Rectangle([bottom[0] + (top[0]-bottom[0])/6.0, bottom[1] + (top[1]-bottom[1])/6.0, \
			bottom[2] + (top[2]-bottom[2])/6.0, bottom[3] + (top[3]-bottom[3])/6.0, \
			bottom[0] + (top[0]-bottom[0])/3.0, bottom[1] + (top[1]-bottom[1])/3.0, \
			bottom[2] + (top[2]-bottom[2])/3.0, bottom[3] + (top[3]-bottom[3])/3.0])
		layer5.id = "Layer5"
		layer4 = n3Rectangle([bottom[0] + (top[0]-bottom[0])/3.0, bottom[1] + (top[1]-bottom[1])/3.0, \
			bottom[2] + (top[2]-bottom[2])/3.0, bottom[3] + (top[3]-bottom[3])/3.0, \
			bottom[0] + (top[0]-bottom[0])/2.0, bottom[1] + (top[1]-bottom[1])/2.0, \
			bottom[2] + (top[2]-bottom[2])/2.0, bottom[3] + (top[3]-bottom[3])/2.0])
		layer4.id = "Layer4"
		layer3 = n3Rectangle([bottom[0] + (top[0]-bottom[0])/2.0, bottom[1] + (top[1]-bottom[1])/2.0, \
			bottom[2] + (top[2]-bottom[2])/2.0, bottom[3] + (top[3]-bottom[3])/2.0, \
			bottom[0] + (top[0]-bottom[0])/1.5, bottom[1] + (top[1]-bottom[1])/1.5, \
			bottom[2] + (top[2]-bottom[2])/1.5, bottom[3] + (top[3]-bottom[3])/1.5])
		layer3.id = "Layer3"
		layer2 = n3Rectangle([bottom[0] + (top[0]-bottom[0])/1.5, bottom[1] + (top[1]-bottom[1])/1.5, \
			bottom[2] + (top[2]-bottom[2])/1.5, bottom[3] + (top[3]-bottom[3])/1.5, \
			bottom[0] + (top[0]-bottom[0])/1.2, bottom[1] + (top[1]-bottom[1])/1.2, \
			bottom[2] + (top[2]-bottom[2])/1.2, bottom[3] + (top[3]-bottom[3])/1.2])
		layer2.id = "Layer2"
		layer1 = n3Rectangle([bottom[0] + (top[0]-bottom[0])/1.2, bottom[1] + (top[1]-bottom[1])/1.2, \
			bottom[2] + (top[2]-bottom[2])/1.2, bottom[3] + (top[3]-bottom[3])/1.2, \
			top[0], top[1], top[2], top[3]])
		layer1.id = "Layer1"
		self.areas = [layer1, layer2, layer3, layer4, layer5, layer6]
		
	def new_neurons(self, simulation_step, area):
		""" Tells the simulation how many neurons will be added in this step
			for the given area
		"""
		if simulation_step < 50 :
			if   area.id=="Layer1":
				return 1
			elif area.id == "Layer2":
				return 1
			elif area.id == "Layer3":
				return 1
			elif area.id == "Layer4":
				return 1
			elif area.id == "Layer5":
				return 1
			else : #Layer6
				return 1
		else :
			return 0

	def next_neuron(self, simulation_step, area):
		""" Returning a new instance of neuron. This will be placed on the 
			given area. If it collides, it won't be placed, at all.
		"""
		x = random.random() * (area.right_front_lower - area.left_front_lower).x + area.left_front_lower.x
		y = random.random() * (area.left_back_lower - area.left_front_lower).y + area.left_front_lower.y	
		z = random.random() * (area.right_front_top - area.left_front_lower).z + area.left_front_lower.z			
		if   area.id=="Layer1":
			return ShortDistanceNeuron(n3Point(x, y, z), self.bounding_box, self.area_type, self.rotation_angle)
		elif area.id == "Layer2":
			return ShortDistanceNeuron(n3Point(x, y, z), self.bounding_box, self.area_type, self.rotation_angle)
		elif area.id == "Layer3":
			return ShortDistanceNeuron(n3Point(x, y, z), self.bounding_box, self.area_type, self.rotation_angle)
		elif area.id == "Layer4":
			return LongDistanceNeuron(n3Point(x, y, z), self.area_target)
		elif area.id == "Layer5":
			return ShortDistanceNeuron(n3Point(x, y, z), self.bounding_box, self.area_type, self.rotation_angle)
		else : #Layer6
			return ShortDistanceNeuron(n3Point(x, y, z), self.bounding_box, self.area_type, self.rotation_angle)
	
	
p1 = n3Point(0, 0, 0)
p2 = n3Point(0, 50, 0)
p3 = n3Point(50, 50, 0)
p4 = n3Point(50, 0, 0)
p5 = n3Point(0, 0, 30)
p6 = n3Point(0, 50, 30)
p7 = n3Point(50, 50, 30)
p8 = n3Point(50, 0, 30)

r1 = n3Rectangle([p1, p2, p3, p4, p5, p6, p7, p8])
r1.translate(r1.middle * -1 + n3Point(0, 0, 75))
lg1 = LayersNeuronGenerator(r1, "A1", "A2", (0, 0, 0))

r2 = n3Rectangle([p1, p2, p3, p4, p5, p6, p7, p8])
r2.translate(r2.middle * -1 + n3Point(75, 0, 0))
r2.rotate_by_self(0, 90, 0)
lg2 = LayersNeuronGenerator(r2, "A2", "A3", (0, 90, 0))

r3 = n3Rectangle([p1, p2, p3, p4, p5, p6, p7, p8])
r3.translate(r3.middle * -1 + n3Point(-75, 0, 0))
r3.rotate_by_self(0, -90, 0)
lg3 = LayersNeuronGenerator(r3, "A3", "A4", (0, -90, 0))

r4 = n3Rectangle([p1, p2, p3, p4, p5, p6, p7, p8])
r4.translate(r4.middle * -1 + n3Point(0, 0, -75))
r4.rotate_by_self(0, 180, 0)
lg4 = LayersNeuronGenerator(r4, "A4", "A5", (0, 180, 0))

r5 = n3Rectangle([p1, p2, p3, p4, p5, p6, p7, p8])
r5.translate(r5.middle * -1 + n3Point(0, 75, 0))
r5.rotate_by_self(-90, 0, 0)
lg5 = LayersNeuronGenerator(r5, "A5", "A6", (-90, 0, 0))

r6 = n3Rectangle([p1, p2, p3, p4, p5, p6, p7, p8])
r6.translate(r6.middle * -1 + n3Point(0, -75, 0))
r6.rotate_by_self(90, 0, 0)
lg6 = LayersNeuronGenerator(r6, "A6", "A1", (90, 0, 0))
	
#inhibitors
inh =  [n3Sphere(n3Point(58, 58, 58), 75), n3Sphere(n3Point(-58, 58, 58), 75), n3Sphere(n3Point(58, -58, 58), 75), n3Sphere(n3Point(-58, -58, 58), 75), \
n3Sphere(n3Point(58, 58, -58), 75), n3Sphere(n3Point(-58, 58, -58), 75), n3Sphere(n3Point(58, -58, -58), 75), n3Sphere(n3Point(-58, -58, -58), 75)]	
	
s = Simulation([lg1, lg2, lg3, lg4, lg5, lg6])
s.set_bounding_area = n3Sphere(n3Point(0, 0, 0), 100)
s.set_chemical_gradient_field(np.load("chemical_gradient_field.npy"))
s.set_verbosity(1)
#s.simulate()
#s.print_simulation_meta_data()
#print s.get_statistics()
#d = s.get_distance_matrix()
