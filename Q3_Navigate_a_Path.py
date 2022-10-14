# imports statements (signal for timing and defaultdict for storing graph)
import signal
from collections import defaultdict

# opens input file and outputs file
input_file = open("navigate.in")
output_file = open("navigate.out", "w")

# accepts and type-casts inputs
n, m, k, p = [int(i) for i in input_file.readline().split()]

# defines global variables
max_exec_time = 1.0
points = []
polys = []
steps = []

# represents coordinate grid using "0" for points
graph = [[0 for row in range(m)] for col in range(n)]

# iterates through next k lines for
# coordinates of points to hit
for i in range(0,k):

	# accepts, type-casts, and appends coordinates to points
	x, y = [int(j) for j in input_file.readline().split()]
	points.append((x, y))

# iterates through p lines for polygons
for i in range(0, p):

	# accepts and type-casts number of vertices in a polygon
	num_verts = int(input_file.readline())
	verts = []

	# iterates through next num_verts lines for vertices
	for j in range(0, num_verts):

		# accepts, type-casts, and appends vertices to verts
		x_poly, y_poly = [int(num) for num in input_file.readline().split()]
		verts.append((x_poly, y_poly))

	# appends created list to polys
	polys.append(verts)

# marks 1 in graph for every point that must be hit
def mark_points_to_hit():
	for x, y in points:
		graph[x - 1][y - 1] = 1

# marks 2 graph for every point in polygon
def mark_poly():
	global n
	for i in range(len(polys)):
		temp_graph = [[0 for row in range(m)] for col in range(n)]
		verts = list(polys[i])
		verts.append(polys[i][0])

		# marks edges of polygon
		for num in range(len(verts) - 1):
			x1, y1 = verts[num]
			x2, y2 = verts[num + 1]

			x1 -= 1
			y1 -= 1
			x2 -= 1
			y2 -= 1

			# marks vertical edges of polygon
			if x1 == x2:
				for y in range(min(y1, y2), max(y1, y2) + 1):
					temp_graph[x1][y] = 2

			# marks horizontal edges of polygon
			if y1 == y2:
				for x in range(min(x1, x2), max(x1, x2) + 1):
					temp_graph[x][y1] = 2

			# marks edges of polygon with slope 1 or -1
			if (x1 != x2) and (y1 != y2):
				delta_x = 1 if x1 < x2 else -1
				delta_y = 1 if y1 < y2 else -1
				y = y1
				x = x1
				for i in range(min(x1, x2), max(x1, x2) +1):
					temp_graph[x][y] = 2
					x += delta_x
					y += delta_y

		# marks points inside polygon traversing horizontally
		for i in range(m):
			b_in_poly = False
			b_start = -1
			for j in range(n):
				if temp_graph[j][i] == 2:
					if not b_in_poly:
						b_in_poly = True
						b_start = j
					else:
						b_in_poly = False
						for num in range(b_start, j):
							temp_graph[num][i] = max(1, temp_graph[num][i])

				# marks points inside polygon traversing vertically
				for i in range(n):
					b_in_poly = False
					b_start = -1
					for j in range(m):
						if temp_graph[i][j] == 2:
							if not b_in_poly:
								b_in_poly = True
								b_start = j
							else:
								b_in_poly = False
								for num in range(b_start, j):
									temp_graph[i][num] = min(2, temp_graph[i][num] + 1)

		# overlays temp_graph to graph
		for i in range(m):
			for j in range(n):
				if temp_graph[j][i] == 2:
					graph[j][i] = 2

# finds all possible steps
def find_all_steps():
	for i in range(m):
		for j in range(n):
			if graph[j][i] != 2:

				# checks if West (W) and appends next coordinate if traversable
				# in single-number format to steps
				if j - 1 >= 0:
					if graph[j - 1][i] != 2:
						steps.append((j * m + i, (j - 1) * m + (i)))

				# checks if North (N) and appends next coordinate if traversable
				# in single-number format to steps
				if i + 1 < m:
					if graph[j][i + 1] != 2:
						steps.append((j * m + i, (j) * m + (i + 1)))

				# checks if North-East (NE) and appends next coordinate if traversable
				# in single-number format to steps
				if j + 1 < n and i + 1 < m:
					if graph[j + 1][i + 1] != 2:
						steps.append((j * m + i, (j + 1) * m + (i + 1)))

				# checks if East (E) and appends next coordinate if traversable
				# in single-number format to steps
				if j + 1 < n:
					if graph[j + 1][i] != 2:
						steps.append((j * m + i, (j + 1) * m + (i)))

				# checks if South-East (SE) and appends next coordinate if traversable
				# in single-number format to steps
				if j + 1 < n and i - 1 >= 0:
					if graph[j + 1][i - 1] != 2:
						steps.append((j * m + i, (j + 1) * m + (i - 1)))

				# checks if South (S) and appends next coordinate  if traversable
				# in single-number format to steps
				if i - 1 >= 0:
					if graph[j][i - 1] != 2:
						steps.append((j * m + i, (j) * m + (i - 1)))

				# checks if South-West (SW) and appends next coordinate if traversable
				# in single-number format to steps
				if j - 1 >= 0 and i - 1 >= 0:
					if graph[j - 1][i - 1] != 2:
						steps.append((j * m + i, (j - 1) * m + (i - 1)))

				# checks if North-West (NW) and appends next coordinate if traversable
				# in single-number format to steps
				if j - 1 >= 0 and i + 1 < m:
					if graph[j - 1][i + 1] != 2:
						steps.append((j * m + i, (j - 1) * m + (i + 1)))


# represents directed graph using adjacency list representation
# modeled after code found in: https://www.geeksforgeeks.org/find-paths-given-source-destination/
class Graph:

	def __init__(self, vertices):
		# number of vertices
		self.V = vertices
		self.short_path = []
		self.len_short_path = self.V + 1
		self.counter = 0

		# default dictionary to store graph
		self.graph = defaultdict(list)

	# resets graph
	def reset(self):
		self.short_path = []
		self.len_short_path = self.V + 1
		self.counter = 0

	# adds edge to graph
	def add_edge(self, u, v):
		self.graph[u].append(v)

	# returns next point to go to in single-number format
	def next_step(self, x, y, id):
		if id == 0:  # A
			return (x) * m + (y + 1)
		if id == 1:  # RA
			return (x + 1) * m + (y + 1)
		if id == 2:  # LA
			return (x - 1) * m + (y + 1)
		if id == 3:  # R
			return (x + 1) * m + (y)
		if id == 4:  # L
			return (x - 1) * m + (y)
		if id == 5:  # RB
			return (x + 1) * m + (y - 1)
		if id == 6:  # LB
			return (x - 1) * m + (y - 1)
		if id == 7:  # B
			return (x) * m + (y - 1)
		return None

	# returns preferred order of steps given direction of next point to hit
	# relative to current point
	def get_optimal_order_of_steps(self, u, d, def_steps):
		steps = []
		x1, y1 = ((int(u) // int(m)), (u % m))
		x2, y2 = ((int(d) // int(m)), (d % m))

		# N --> 0, NE --> 1, NW --> 2, E --> 3,
		# W --> 4, SE --> 5, SW --> 6, S --> 7
		preferred_order = []

		# direction is N -->  preferred order is N, NE, NW, E, W, SE, SW, S
		if (x1 == x2) and (y1 < y2):
			preferred_order = [0, 1, 2, 3, 4, 5, 6, 7]

		# direction is S -->  preferred order is S, SW, SE, W, E, NW, NE, N
		if (x1 == x2) and (y1 > y2):
			preferred_order = [7, 6, 5, 4, 3, 2, 1, 0]

		# direction is E -->  preferred order is E, SE, NE, S, N, SW, NW, W
		if (y1 == y2) and (x1 < x2):
			preferred_order = [3, 5, 1, 7, 0, 6, 2, 4]

		# direction is W -->  preferred order is W, NW, SW, N, S, NE, SE, E
		if (y1 == y2) and (x1 > x2):
			preferred_order = [4, 2, 6, 0, 7, 1, 5, 3]

		# direction is NE -->  preferred order is NE, E, N, SE, NW, S, W, SW
		if (x1 < x2) and (y1 < y2):
			preferred_order = [1, 3, 0, 5, 2, 7, 4, 6]

		# direction is SW -->  preferred order is SW, W, S, NW, SE, N, E, NE
		if (x1 > x2) and (y1 > y2):
			preferred_order = [6, 4, 7, 2, 5, 0, 3, 1]

		# direction is SE -->  preferred order is SE, E, S, NE, SW, N, W, NW
		if (x1 < x2) and (y1 > y2):
			preferred_order = [5, 3, 7, 1, 6, 0, 4, 2]

		# direction is NW -->  preferred order is NW, W, N, SW, NE, S, E, SE
		if (x1 > x2) and (y1 < y2):  # Going Left Above; preferred order is LA, L, A, LB, RA, B, R, RB
			preferred_order = [2, 4, 0, 6, 1, 7, 3, 5]

		# finds preferred next point and appends to created list
		for num in preferred_order:
			next = int(self.next_step(x1, y1, num))
			if next in def_steps:
				steps.append(next)
		if len(def_steps) != len(steps):
			print("WARNING: something could be wrong in determining the optimal next step")

		return steps

	# recursive funstion to help find all paths from source to destination
	# visited keeps track of vertices in current path
	# path stores actual vertices and path_index is current index in path
	def find_all_paths_util(self, u, d, visited, path):

		# marks current node as visited and stores in path
		visited[u] = True
		path.append(u)

		# checks if current point is same as destination
		# (i.e. destination is reached)
		if u == d:
			self.counter += 1
			self.len_short_path = len(path)
			self.short_path = list(path)

		# recurs for all points adjacent to current point
		# if current point is not same as destination
		else:
			next_steps = list(self.get_optimal_order_of_steps(u, d, self.graph[u]))

			# recurs only if path is less than shortest path found
			# in order to reduce run-time and look for optimal paths only
			for i in next_steps:
				if not visited[i]:
					if len(path) < self.len_short_path - 1:
						self.find_all_paths_util(i, d, visited, path)

		# removes current vertex from path and marks as unvisited
		path.pop()
		visited[u] = False

	# finds all paths from source to destination
	def find_all_paths(self, s, d):

		# marks all vertices as not visited
		visited = [False] * (self.V)
		path = []

		# calls the recursive function to find all paths
		self.find_all_paths_util(s, d, visited, path)

	# helps if debugging is needed
	def debug(self):
		print("Graph Nodes: ")
		print(self.graph)


# helps implement time limit for run-time
# modelled after code found in: https://stackoverflow.com/questions/25027122/break-the-function-after-certain-time
class TimeoutException(Exception):
	pass

# raises TimeoutException if run-time takes too long
def timeout_handler(signum, frame):
	raise TimeoutException

# creates signal that implements timeout
signal.signal(signal.SIGALRM, timeout_handler)

# calls helper functions (mark_points_to_hit,
# mark_poly, find_all_steps, and add_edge) and
# creates new Graph i_Graph
mark_points_to_hit()
mark_poly()
find_all_steps()
i_Graph = Graph(m * n)
for step in steps:
	i_Graph.add_edge(step[0], step[1])

# helps remove multiple points when combining shortest paths
last_x = -1
last_y = -1

# iterates through points to hit
for i in range(len(points) - 1):
	i_Graph.reset()

	# defines x and y coordinates for points to hit and
	# modifies them to fit in a grid starting at (1, 1)
	x1, y1 = points[i]
	x2, y2 = points[i + 1]
	x1 -= 1
	y1 -= 1
	x2 -= 1
	y2 -= 1

	# sets execution time to either a total of max_exec_time or
	# 0.1 seconds each, whichever is longer
	exec_time = max(max_exec_time / len(points), 0.1)
	signal.setitimer(signal.ITIMER_REAL, exec_time, 0.01)

	# tries to find shortest path, but if not possible in
	# given execution time, implements TimeoutException
	try:
		i_Graph.find_all_paths(x1 * m + y1, x2 * m + y2)
	except TimeoutException:
		pass
	else:
		signal.alarm(0)

	# appends found path to short_path
	finally:
		short_path = []
		for node in i_Graph.short_path:
			short_path.append(((node // m) + 1, (node % m) + 1))

	# checks if point at start of new path is repeated
	# (uses variables last_x and last_y defined above)
	for x, y in short_path:
		if last_x == x and last_y == y:
			continue

		# writes coordinates in output file
		output_file.write(str(x) + " " + str(y) + "\n")

		# changes last_x and last_y to check if they are repeated in next path
		last_x = x
		last_y = y






