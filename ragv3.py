# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 17:28:16 2015
Last Update on Fri Feb 20 16:02:53 2015

This is the third version of a python script modeling and simulating the 
random axon growth model.
See following papers:
#TODO: papers

@author: Jan Zelmer
@mail: jzelmer@cbs.mpg.de

Future possible improvements:
Add bounding box for geometric classes
Adding new geometric objects and its subroutines
Adding get_middle for geometric objects
Discuss and decide the use of different metrics
"""

# coding: utf8

#zest.releaser
#@watman

import random
from math import sqrt
from math import pi
from math import cos
from math import sin
from math import e
from copy import deepcopy
import numpy as np

# initializing the generator of random numbers
random.seed()


class n2Point():
	"""Actually this is a location vector. And therefore some vector operations are supported.
	"""		

	def __init__ (self, x, y):
		""" Simple constructor.
		"""
		self.x = x
		self.y = y
		
		
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

	def length (self):
		""" Returns the length of this vector.
		"""
		return sqrt ( self.x * self.x + self.y * self.y)

	def distance_to (self, p):
		""" Returns the distance between this and the given point by computing
			the length of the resulting vector.			
		"""
		return (self - p).length()

	def get_position(self):
		""" Returns its position which is basically itself.
		"""
		return self


class n2Line ():
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
		self.direction_vector_norm = self.direction_vector / self.direction_vector.length()

	def push_delta (self, p):
		""" Pushes the line by adding the coordinates form the point p to the
			actual head as the new head. Finally setting the tail to the old
			head.
		"""
		self.tail = self.head
		self.head = self.head + p
		self.direction_vector = self.head - self.tail
		self.direction_vector_norm = self.direction_vector / self.direction_vector.length()

	def get_middle (self):
		""" Returns the point representing the middle.
		"""
		return self.tail + (self.head - self.tail) / 2 

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
		r = self.direction_vector.scalar_product(v) / self.direction_vector.length()
		if r < 0:
			return self.tail.distance_to(p)
		elif r > self.direction_vector.length():
			return self.head.distance_to(p)
		else:
			return p.distance_to(self.direction_vector_norm * r + self.tail)

	def get_length (self):
		""" Returns the length of this line
		"""
		return (self.head - self.tail).length()


class n2Rectangle ():
	""" This class implements a rectangle and is implemented based on n2Point"""	
	
	def __init__ (self, p1, p2): 
		""" Simple constructor with following ordering of the points.
		"""		
		p3 = n2Point(p1.x, p2.y)
		p4 = n2Point(p2.x, p1.y)
		self.lbPoint, self.luPoint, self.ruPoint, self.rbPoint = self.order(p1,p2,p3,p4)
				
	def get_points (self): 
		""" Simple getter method. Returns the points of the rectangle in a list.
			The first element is the point in the left bottom, the second
			is in the left upper corner, the next in the right upper corner
			and finally the last is in the right bottom.		
		"""
		return [self.lbPoint, self.luPoint, self.ruPoint, self.rbPoint,]
		
	def order(self, p1,p2,p3,p4):
		"""Orders points according to: Most left Bottom Point is first, then
			left upper, then right upper and finally right bottom. 	The 
			points are determined via the distance to coordinate origin
			and a point on the x axis lying to the right of the rectangle.
		"""
		parray = [p1,p2,p3,p4]
		l1 = (p1.distance_to (n2Point(0,0)),0) 
		l2 = (p2.distance_to (n2Point(0,0)),1)
		l3 = (p3.distance_to (n2Point(0,0)),2)
		l4 = (p4.distance_to (n2Point(0,0)),3)
		larray = [l1,l2,l3,l4]
		larray.sort(lambda x,y: -1 if x < y else 1) 
		lb = parray[larray[0][1]] 
		ru = parray[larray[3][1]] 
		rarray = [p1,p2,p3,p4]
		rarray.sort(lambda x,y: -1 if x.x < y.x else 1) 
		maxx = rarray[3].x 
		l1 = (p1.distance_to (n2Point(maxx+10,0)),0) 
		l2 = (p2.distance_to (n2Point(maxx+10,0)),1)
		l3 = (p3.distance_to (n2Point(maxx+10,0)),2)
		l4 = (p4.distance_to (n2Point(maxx+10,0)),3)
		larray = [l1,l2,l3,l4]
		larray.sort(lambda x,y: -1 if x < y else 1)
		rb = parray[larray[0][1]]
		lu = parray[larray[3][1]]
		return  [lb,lu,ru,rb]
	
	def get_middle (self):
		""" Returns the point in the middle of the rectangle
		"""
		return self.lbPoint + (self.ruPoint - self.lbPoint) / 2


class n2AxisParallelRectangle (n2Rectangle):
	""" The lines of the rectangle are parallel to the axis.
	"""	
	
	def is_inside (self, p1):
		return 0 <= p1.x <= self.ruPoint.x and 0 <= p1.y <= self.ruPoint.y
		

class GeoTree ():
	"""Index Structure for accessing spatial information in log(n) time.
		It uses 4 (and splitting) axis parallel rectangles as way to access
		spatial data. It works like a tree on those rectangles
	"""
	
	def __init__(self, limit, rect, minimum_area_size):
		""" Simple constructor. Just defining some variables needed.
		"""
		self.limit = limit
		self.container = []
		self.rect = rect 
		self.is_leaf = True
		self.minimum_area_size = minimum_area_size
		
	@property
	def middle (self):
		""" Returns the middle of the underlying rectangle as property
		"""
		return self.rect.get_middle()

	@property
	def size_of_area (self):
		""" Returns the length of one side of the underlying rectangle = size
		"""
		return (self.rect.get_points()[1] - self.rect.get_points()[0]).length()

	def add_element (self, element):
		""" Add a new element. If this isn't a leaf, the tree will traversed
			a step and then tried to add it again. 
			If this container get's full and is still twice above the minimum
			size, it is splitted and looses its leaf property. The content
			is divided on the different subtrees.
		"""				
		if self.is_leaf:
			self.container.append(element)
			if len(self.container) > self.limit and self.minimum_area_size < ( self.size_of_area / 2):
				#split, if container is full AND minimum size isn't reached
				self.is_leaf = False
				corner_points = self.rect.get_points()
				lower_x = corner_points[0].x
				middle_x = self.middle.x
				higher_x = corner_points[2].x
				lower_y = corner_points[0].y
				middle_y = self.middle.y
				higher_y = corner_points[2].y			
				self.q1 = GeoTree(self.limit, n2Rectangle(n2Point(lower_x, lower_y),   n2Point(middle_x, middle_y)), self.minimum_area_size)
				self.q2 = GeoTree(self.limit, n2Rectangle(n2Point(lower_x, middle_y),  n2Point(middle_x, higher_y)), self.minimum_area_size)
				self.q3 = GeoTree(self.limit, n2Rectangle(n2Point(middle_x, middle_y), n2Point(higher_x, higher_y)), self.minimum_area_size)
				self.q4 = GeoTree(self.limit, n2Rectangle(n2Point(middle_x, lower_y),  n2Point(higher_x, middle_y)), self.minimum_area_size)
				for el in self.container:					
					self.get_quadrant(el).add_element(el)					
				self.container = []
		else :			
			self.get_quadrant(element).add_element(element)
				
	def get_elements(self, element):
		""" Get all elements of this leaf or traverse a step into the tree.
		"""
		if self.is_leaf :
			return deepcopy(self.container)
		else :
			return self.get_quadrant(element).get_elements(element)
				
	def get_elements_within_radius(self, element, radius):
		""" Method which returns all elements, which lay inside the radius of 
			the element and more. The corner points of the bounding box of
			the circle are used for getting all elements in the respectively
			leafs.
		"""
		if self.is_leaf :
			return deepcopy(self.container)
		else :
			lower_x = element.get_position().x - radius
			higher_x = element.get_position().x + radius
			lower_y = element.get_position().y - radius
			higher_y = element.get_position().y + radius
			result               = self.get_elements(n2Point(lower_x, lower_y))
			result[len(result):] = self.get_elements(n2Point(higher_x, lower_y))
			result[len(result):] = self.get_elements(n2Point(higher_x, higher_y))
			result[len(result):] = self.get_elements(n2Point(lower_x, higher_y))
			# return list(set(result))
			return result
			
	def get_quadrant (self, element):
		""" Method which gets the right quadrant according to the element
		"""
		x = element.get_position().x
		y = element.get_position().y			
		if x <= self.middle.x:
			if y <= self.middle.y:
				return self.q1
			else:
				return self.q2
		else:
			if y <= self.middle.y:
				return self.q4
			else:
				return self.q3

	
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
	

class Simulation ():
	""" This is our "Central processing unit". It manages the whole simulation
		including distance matrix, geotree, neuron generators and neurons.
	"""
	
	def __init__ (self, neuron_generators):
		""" Normal constructor. Just defining some variables we need. And 
			constructing our super area, which contains all other areas.
		"""
		self.generators = neuron_generators		
		self.neurons = []
		self.dmatrix = Distances(100)
		self.max_radius = 0
		min_x, min_y, max_x, max_y = 100000, 100000, 0, 0
		for generator in neuron_generators:
			self.areas = generator.get_areas()
			for area in self.areas:
				min_x = min( min_x, area.get_points()[0].x)
				min_y = min( min_y, area.get_points()[0].y)				
				max_x = max( max_x, area.get_points()[2].x)				
				max_y = max( max_y, area.get_points()[2].y)
		self.super_area = n2AxisParallelRectangle(n2Point(min_x, min_y), n2Point(max_x, max_y))
		self.tree = GeoTree(15, self.super_area, 2.0)		
		self.max_distance = sqrt( pow(max_x - min_x, 2) + pow(max_y - min_y, 2))
		self.dist_counter = 0
		self.simulation_step_counter = 0

	def simulation_step (self):
		""" This method simulates a simulation step. DO NOT CALL DIRECTLY.
		"""		
		self.simulation_step_counter += 1
		added_neurons = 0
		for generator in self.generators:		
			for area in generator.get_areas():
				neurons_to_add = generator.new_neurons(self.simulation_step_counter, area)
				added_neurons += neurons_to_add
				for i in xrange(neurons_to_add):				
					nneuron = generator.next_neuron(self.simulation_step_counter, area)
					while not self.is_free(nneuron, nneuron.radius):
						nneuron = generator.next_neuron(self.simulation_step_counter, area)			
					self.max_radius = max(nneuron.get_radius(), self.max_radius)
					nneuron.axon = n2Line(nneuron.get_position())
					nneuron.active = True
					nneuron.dist_index = self.dist_counter
					self.dist_counter += 1
					self.tree.add_element(nneuron)
					self.neurons.append(nneuron)		
		for neuron in self.neurons:			
			if neuron.active:
				# grow 
				neuron.axon.push_delta(neuron.grow())
				# make connections
				if neuron.can_put_connection():					
					center = neuron.axon.get_middle()
					radius = neuron.axon.get_length() + self.max_radius
					for ntest in self.tree.get_elements_within_radius(center, radius):
						if ntest is not neuron and 	ntest.can_receive_connection(neuron) and \
							neuron.axon.compute_distance(ntest.get_position()) < ntest.get_radius() and \
							self.dmatrix.add(neuron, ntest, self.generators[0].metric.compute_distance):
							neuron.put_connection(ntest)
							ntest.receive_connection(neuron, neuron.axon)
				neuron.active = neuron.can_put_connection() and self.super_area.is_inside(neuron.axon.get_head())
		print "Step %i: Added %i new neurons." %(self.simulation_step_counter, added_neurons)

	def simulate (self):
		""" Main method to start the simulation.
		"""
		print "Adding Neurons and growing Axons"
		min_iterations = 0
		for gen in self.generators:
			min_iterations = max(min_iterations, gen.get_minimum_iterations())
		for i in xrange(min_iterations):
			self.simulation_step()
		print "Finishing Growth"
		while not self.finished():			
			self.simulation_step()
		self.print_statistics()
		
	def print_statistics (self):
		""" Prints some statistics onto the console.
		"""
		print "%i Neurons were added and they made %i connections" %(len(self.neurons), self.dmatrix.compute_statistics(self.max_distance, 1)[0])
		print "The arithmetric mean of all connections is %f" %(self.dmatrix.compute_arithmetic_mean_length())
		print "The median connection is %f long"%(self.dmatrix.get_median_length())
		print "The shortest connection is %f short and the longest is %f long" %(self.dmatrix.get_min_length(), self.dmatrix.get_max_length())
	
	
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
		neurons = self.tree.get_elements_within_radius(point, radius)
		for neuron in neurons:			
			if self.generators[0].metric.compute_distance(point.get_position(), neuron.get_position()) < radius:
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
		if partition_bins > 0:			
			return self.dmatrix.compute_statistics(self.max_distance, partition_bins)
		else :
			return []

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

	def saveModelNeuronsWithStraightAxonsAsFile (self, path):
		""" Saves our model (neurons and respecting axons) as coordinates in
			a text file
		"""
		with open(path,'w') as f:
			for n in self.neurons:
				s = n.get_position().__str__().strip('()') + ", " + n.axon.line.get_head().__str__().strip('()') +  "\n"
				f.write(s)


class TwoDimEuclidMetric ():
	""" Simple 2 dimensional euclidean metric to compute the distance between 
	to points
	"""
	
	def compute_distance (self, p1, p2):
		a = p1.x - p2.x
		b = p1.y - p2.y
		return sqrt (a*a + b*b)

class TrivialNeuronGenerator ():
	""" Very trivial generator. Every method implemented is used by the simulation.
		Every generator class must either directly or indirectly inherit from this.
	"""
	
	def __init__(self):
		""" Trivial constructor. Metric and areas should generally exist.
		"""
		self.i = 0
		self.metric = TwoDimEuclidMetric()
		self.areas = [n2Rectangle(n2Point(0,0),n2Point(50,50))]
		
	def new_neurons(self, simulation_step, area):
		""" Tells the simulation how many neurons will be added in this step
			for the given area
		"""
		return 2

	def get_areas (self):
		""" Returns the areas for this generator
		"""
		return self.areas
		
	def get_metric (self):
		""" Returns the metric for distance computation this generator defines
		"""
		return self.metric

	def next_neuron(self, simulation_step, area):
		""" Returning a new instance of neuron. This will be placed on the 
			given area. If it collides, it won't be placed, at all.
		"""
		self.i += 1
		return TrivialNeuron(self.i - 1)

	def get_minimum_iterations (self):
		""" As stated below this tells the simulation how much minimum simulation
			steps have to be simulated.
		"""
		return 5

class TrivialNeuron ():
	""" Very trivial neuron. Used to either define methods and properties we need
		to run the simulation and to be inherited from.
	"""
	
	def __init__ (self, i):
		""" Trivial constructor
		"""
		self.position = n2Point(i,i)
		self.radius = 1.0

	@property
	def x (self):
		""" Returns the x-coordinate of the neuron
		"""
		return self.position.x
	
	@property
	def y (self):
		""" Returns the y-coordinate of the neuron
		"""
		return self.position.y
	
	def get_radius(self):
		""" Returns the radius of the neuron
		"""
		return self.radius
		
	def get_position(self):
		""" Returns the position of the neuron
		"""
		return self.position

	def grow (self):
		""" Is called everytime this neuron shall grow. The return value
			is a direction vector.
		"""
		return n2Point(0.2,0.3)
		
	def can_put_connection(self):
		""" Is called, when the axon from this neuron found an neuron and 
			wants to establish a connection OR whether the axon to this 
			neuron can still grow.			
		"""
		return True
	
	def put_connection(self, target_neuron):
		""" Is called, when the outgoing connection is made.
		"""
		pass
		
	def can_receive_connection(self):
		""" Is called, when another neuron tries to make a connection.
			The return value states whether this is okay or not
		"""
		return True
		
	def receive_connection(self, origin_neuron, axon):
		""" Is called, when another neuron makes a connection.
		"""
		pass


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
		
""" 	-----------------------------------------------------------------------------
	-----------------------------------------------------------------------------
	From below this point, the code could and should be altered to fit your needs
	
	Explanation: A neuron generator produces neurons and give it to the simulation
	instance. The simulation instance cares about all technical details, like 
	computing distances, do not allow overlapping neurons, make connections 
	between axons and neurons, doing statistics and much more.
	
	The neuron generator holds all crucial information about placing neurons:
	when, where and especially which. So you can have many different neuron 
	classes acting differently, but you can only have one generator. This 
	generator then specifies, which class to use and when etc.
	
	Every new neuron class has to implement the methods defined in TrivialNeuron
	(see above). This is easily done by inherit from it. See is shown in the 
	definition of SimpleNeuron. As you can see TrivialNeuron is in the brackets
	and therefore SimpleNeuron inherits all methods and property from TrivialNeuron.
	You can now reimplement the methods you like to alter and even create new
	methods and properties.
	The same is true for the SimpleNeurongenerator.
	
	In general the simulation ends, when they aren't any active neurons on the
	field. A neuron is active, as soon as is placed and will be deactived if 
	one of the 2 cases hold: its axon reaches the border or it can't make a 
	connection any more. A once deactivated neuron can't be activated again.
	Though the simulation would therefor never start, you have to provde a 
	minimum iteration length, which will definitly be simulated. After that 
	the above criteria holds.
	
	As you can see in the last lines: The Simulation class is the main object
	which organizes all and calls the rest of the code. Thats why you create
	an instance of it providing your generator. After that you start the
	simulation by calling simulate(). And finally cou can call any of the 
	statistics methods provided.
	-----------------------------------------------------------------------------
	-----------------------------------------------------------------------------
"""


class SimpleNeuron (TrivialNeuron):
	""" A class representing a neuron in our simulation. This neuron is simple
		and inherits methods and properties from TrvialNeuron.
		Its axon growth is linear and static over the time course.
		It can have 10 outgoing and 10 incoming connections.
	"""		
	
	def __init__ (self, position):
		""" Normal constructor. All properties except radius is just for the class itself.
		"""
		angle = random.random() * pi * 2
		self.dx = cos(angle)
		self.dy = sin(angle)
		#self.dx = cos(angle) * growth_speed
		self.position = position
		self.outgoing_limit = 10
		self.outgoing = 0
		self.ingoing_limit = 10
		self.ingoing = 0
		self.radius = 1.0
	
	def grow (self):
		""" Is called everytime this neuron shall grow. The return value
			is a direction vector.
		"""
		return n2Point(self.dx, self.dy)
		
	def can_put_connection(self):
		""" Is called, when the axon from this neuron found an neuron and 
			wants to establish a connection OR whether the axon to this 
			neuron can still grow.			
		"""
		return self.outgoing_limit > self.outgoing

	def put_connection(self, target_neuron):
		""" Is called, when the outgoing connection is made.
		"""
		self.outgoing += 1
		
	def can_receive_connection(self, origin_neuron):
		""" Is called, when another neuron tries to make a connection.
			The return value states whether this is okay or not
		"""
		return self.ingoing_limit > self.ingoing
		
	def receive_connection(self, origin_neuron, axon):
		""" Is called, when another neuron makes a connection.
		"""
		self.ingoing += 1


class SimpleNeuronGenerator (TrivialNeuronGenerator):
	""" A simple neuron generator. It provides only SimpleNeuron and they are
		evenly distributed spatially and gaussion distributed over the time
		course.
	"""
	
	def __init__(self):
		""" Simple constructor. Metric and area "Must be" provided or the 
			simulation will fail. The rest of the variables are just for
			it self.
		"""
		self.metric = TwoDimEuclidMetric()
		self.limit_x = 50
		self.limit_y = 50
		area = n2Rectangle(n2Point(0,0),n2Point(self.limit_x,self.limit_y))
		area.id = "TheOneAndOnly!"
		self.areas = [area]
		self.gauss = GaussAlikeFunction(35,10,5,10,60)

	def get_minimum_iterations (self):
		""" As stated above this tells the simulation how much minimum simulation
			steps have to be simulated.
		"""
		return 60

	def new_neurons(self, simulation_step, area):
		""" Tells the simulation how many neurons will be added in this step
			for the given area
		"""
		if simulation_step >= 10 and simulation_step <= 61 and area.id == "TheOneAndOnly!":
			return Noise.add_gauss_noise(self.gauss.get_value(simulation_step),2)
		else:
			return 0		

	def next_neuron(self, simulation_step, area):
		""" Returning a new instance of neuron. This will be placed on the 
			given area. If it collides, it won't be placed, at all.
		"""
		return SimpleNeuron(n2Point(random.random()*self.limit_x, random.random()*self.limit_y))

""" Till this line everything was just defined. The lines below start the actual simulation
"""		

s = Simulation([SimpleNeuronGenerator()])
s.simulate()
print s.get_statistics()
#d = s.get_distance_matrix()
#s.saveModelNeuronsWithStraightAxonsAsFile("ragv3naxons.txt")