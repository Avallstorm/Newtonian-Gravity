from pygame.locals import *
import numpy as np
import pygame.draw, pygame.time, sys
from math import sqrt, log, sin, cos
from copy import deepcopy as copy

#===============================================================================
#Global variables 

ORIGINX = 0
ORIGINY = 0
"""int: origine points for drawing"""

GRAVEC = (6.67 * (10**-11)) * 1000
"""float: Gravitational constant used in physics calculations"""

VERBOSE = False

#===============================================================================
#Main Class

class gravpoint(object):
	"""Points to be drawn in 3D space with mass and radius that simulate
		gravitational interactions between them


	__init__

	Initializes a gravepoint instance with the following arguments

		Args:
			pos   (:obj:`array` of :obj:`float` len = 3): position of the point
				in 3D space
			mass  (float): Mass used to calculated movement and collisions
			vel   (:obj:`array` of :obj:`float` len = 3): velocity of the point
				use in move method
			imov  (Bool): Whether or not the point's move method modifies it
				AKA if the point is immovable or not
			color (:obj:`tuple` of :obj:`int` len = 3): color in RBG of the 
				point when printed

		Returns: (:obj: gravepoint)

	"""

	def __init__(self, pos, mass, vel, imov = False, color = (255,255,255)):

		self.mass = mass
		self.pos = pos
		self.vel = vel
		self.radius = int((mass)**(1/3)) # Radius used in drawing point
		self.imov = imov
		self.color = color

	def move(self):
		""" Moves the point's location based on it's velocity 

		Note: Does nothing if the point is immovable

		Side effect: Mutates self

		Returns: None

		"""

		if not self.imov:
			self.pos = self.pos + self.vel
			

	def update(self, force):
		""" Updates the velocity of the object based on the force applied

		Args: force (:obj:`array` of :obj:`float` len = 3): force applied to the
			point in 3 dimentions

		Side effect: Mutates self

		Returns: None

		"""

		self.vel = self.vel + force/self.mass

	def __add__(self, other):
		""" overrides the + operator to allow two gravepoints to be added

		Args: other (:obj:'gravepoint'): other point beind added to self

		Notes: adds all values pointwise except for radius which is calculated upon 
			the initialization of the new gravepoint

		Returns: (:obj: gravepoint)

		"""

		if isinstance(other, self.__class__):
			newmass = self.mass + other.mass
			selfin = (self.vel * self.mass)
			otherin = (other.vel * other.mass)
			combmass = (self.mass + other.mass)
			newvel = (selfin + otherin) / combmass
			newimov = self.imov or other.imov
			newpos = self.pos + other.pos/2
			if newimov:
				if self.imov:
					newpos = self.pos
				else:
					newpos = other.pos
			newcol = (	(self.color[0] + other.color[0]) / 2,
						(self.color[1] + other.color[1]) / 2,
						(self.color[2] + other.color[2]) / 2)
			new = gravpoint(newpos,newmass,newvel,
			 imov = newimov, color = newcol)
			return new

	def __eq__(self, other):
		""" overrides the == operator to allow comparison of gravepoints 

		Args: other (:obj:'gravepoint'): other point beind compared to self

		Note: Never brings up a type error but will return false when comparing
			a gravepoint to a non-gravepoint

		Returns: Bool

		"""

		if isinstance(other, self.__class__):
			req1 = self.pos.all() == other.pos.all()
			req2 = self.vel.all() == other.vel.all()
			req3 = self.mass == other.mass
			if req1 and req2 and req3:
				return True
		return False

	def __ne__(self, other):
		""" overrides the != operator to allow comparison of gravepoints

		Args: other (:obj:'gravepoint'): other point beind compared to self

		Returns: Bool

		"""

		return not self == other

	def __str__(self):
		return "pos:{0}, mass:{1}, vel:{2}".format(self.pos,self.mass,self.vel)

#===============================================================================
#helper functions		

def distance(point1, point2):
	""" calculates the distance between two points in 3D space 

	Args:
		point1 (:obj:`array` of :obj:`float` len = 3): position of the point
				in 3D space
		point2 (:obj:`array` of :obj:`float` len = 3): position of the point
				in 3D space

	Returns: Int

	"""
	sum = 0
	for i in range(3):
		sum += (point1[i] - point2[i])**2
	return sqrt(sum)

def normalize(point):
	""" shortens the length of point to 1 while keeping it's verctor direction

	Args:
		point (:obj:`array` of :obj:`float` len = 3): position of the point
				in 3D space

	Side effects: mutates point

	Returns: None

	"""
	if point.all() == 0:
		return point
	zero = [0,0,0]
	dist = distance(point,zero)
	for i in range(3):
		point[i] = point[i]/dist

def merge_detect(gpoints):
	""" detects when two gpoints are close enough, then merges them

	Args: gpoints:(:obj:'list' of :obj: gravepoint)

	Side effects: mutates gpoints according to the merges that are occuring

	Returns: None

	"""
	for i in range(len(gpoints)):
		for j in range(i + 1, len(gpoints)):
			try:
				sumradi = gpoints[i].radius + gpoints[j].radius -1
				rad = gpoints[i].radius + gpoints[j].radius
				dist = distance(gpoints[i].pos, gpoints[j].pos)
				if dist < rad:
					new = gpoints[i] + gpoints[j]
					if i > j:
						gpoints.remove(gpoints[i])
						gpoints[j] = new
					else:
						gpoints.remove(gpoints[j])
						gpoints[i] = new
			except IndexError:
				continue

def gravity(gpoints):
	""" Main engine for applying newtonian gravity

	Args: gpoints:(:obj:'list' of :obj: gravepoint)

	Side effects: mutates gpoints according to the merges that are occuring

	Returns: None

	"""
	forces = []
	for i in range(len(gpoints)):
		#partial forces

		par_forces = []
		mainpoint = gpoints[i]

		for j in range(len(gpoints)):
			otherpoint = gpoints[j]
			if i == j:
				continue

			#calculating force and distance between points i and j
			force = (GRAVEC * (mainpoint.mass*1000 * otherpoint.mass*1000))
			dist = distance(mainpoint.pos, otherpoint.pos)

			#Simulating firctional forces
			if dist == 0:
				force = 0
			else: 
				force = force/(dist**2)
				if dist <= 5:
					force = log(force)

			#creating force vector
			par_force = copy(otherpoint.pos)
			par_force = par_force - mainpoint.pos
			normalize(par_force)
			par_force = par_force*force

			#adding force vector to list of force vectors
			par_forces.append(par_force)

		#adding all the forces of a point together	
		force = np.array([0,0,0])
		for i in range(len(par_forces)):
			force = force + par_forces[i]

		forces.append(force)

	#applying the appropritate for to each point
	for i in range(len(gpoints)):
		gpoints[i].update(forces[i])


def updateall(gpoints):
	""" Applying both merge detecting and gravity funtions and moving points """
	merge_detect(gpoints)
	gravity(gpoints)
	for gpoint in gpoints:
		gpoint.move()

#===============================================================================
#Drawing tools

def draw_3dcircle(surface, colour, a, radius):
	"""Draw a 2D circle around a 3D point.""" 
	ax, ay = a[0]+(a[2]*0.3)+ORIGINX, a[1]+(a[2]*0.3)+ORIGINY
	pygame.draw.circle(surface, colour, (int(ax), int(ay)), radius)

def draw_points(surface, colour, points, radi):
	if type(1) == type(colour):
		for i in range(len(points)):
			draw_3dcircle(surface, colour, points[i], radi[i])
	else:
		for i in range(len(points)):
			draw_3dcircle(surface, colour[i], points[i], radi[i])

#===============================================================================
#Data

#Datapoints used in multiple simulation
gpoints = []
centerp = gravpoint(np.array([0,0,0]), 1500, np.array([0,0,0]), imov = True)

#Datapoints for Line Simulation command
if "square" in sys.argv:
	gp1 = gravpoint(np.array([0,-50,0]), 100, np.array([0,0,0]))
	gp2 = gravpoint(np.array([0,50,0]), 100, np.array([0,0,0]))
	gp3 = gravpoint(np.array([50,0,0]), 100, np.array([0,0,0]))
	gp4 = gravpoint(np.array([-50,0,0]), 100, np.array([0,0,0]))
	gp5 = gravpoint(np.array([0,0,0]), 200, np.array([0,0,1]))
	gpoints += [gp1, gp2, gp3, gp4, gp5]

#Datapoints for Circle Simulation command
if "circle" in sys.argv:
	gpoints += [centerp]
	for i in range(0,360,10):
		cosang = cos(i)
		sinang = sin(i)
		sin2ang  = sin(i + 60)
		cos2ang  = cos(i + 60)
		gpoints += [gravpoint(np.array([cosang*100,sinang*100,0]),
		 50, np.array([cos2ang*10,sin2ang*10,0]),
		  color = (abs(sinang)*255,abs(cosang)*255,abs(sin2ang)*255))]

#Datapoints for Line Simulation command
if "line" in sys.argv:
	gp1 = gravpoint(np.array([0,-50,0]), 100, np.array([0,0,0]))
	gp2 = gravpoint(np.array([0,50,0]), 100, np.array([0,0,0]))
	gp3 = gravpoint(np.array([0,100,0]), 100, np.array([0,0,0]))
	gp4 = gravpoint(np.array([0,-100,0]), 100, np.array([0,0,0]))
	gp5 = gravpoint(np.array([0,200,0]), 100, np.array([0,0,0]))
	gp6 = gravpoint(np.array([0,-200,0]), 100, np.array([0,0,0]))
	gpoints += [gp1, gp2, gp3, gp4, gp5, gp6]

#Datapoints for Stable orbit simulation command
if "orbit" in sys.argv:
	gp1 = gravpoint(np.array([0,0,0]), 2000, np.array([-1,0,0]))
	gp2 = gravpoint(np.array([0,100,0]), 200, np.array([10,0,0]))
	gpoints += [gp1, gp2]

#===============================================================================

def main():

	#Checking for recognized command
	try:
		gpoints
	except NameError:
		print("No command given or command unrecognized")
		sys.exit()

	#Pygame screen setup
	global ORIGINX, ORIGINY
	pygame.init()
	screen = pygame.display.set_mode((640,640))
	ORIGINX = screen.get_width()/2
	ORIGINY = screen.get_height()/2

	#Interation tracking variable
	count = 0

	#Main Pygame loop
	while 1:

		#Exit conditions 
		for event in pygame.event.get():
			if event.type == pygame.QUIT:
				pygame.quit()
				sys.exit()
			elif event.type == pygame.KEYDOWN:
				if event.key == pygame.K_ESCAPE:
					pygame.quit()
					sys.exit()


		if VERBOSE:
			print("Loop {0}".format(count))
			for point in gpoints: print(point)

		#Pulling information from the list of points	
		points = [gpoint.pos for gpoint in gpoints]
		radi = [gpoint.radius for gpoint in gpoints]
		colors = [gpoint.color for gpoint in gpoints]

		#Drawing points to screen 
		draw_points(screen, colors, points, radi)
		event = pygame.event.poll()

		#Resetting screen and updating points
		pygame.display.flip()
		pygame.time.delay(25)
		screen.fill((0,0,0))
		updateall(gpoints)
		count += 1

if __name__ == '__main__':
	main()


