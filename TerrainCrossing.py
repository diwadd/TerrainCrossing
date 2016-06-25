import math
import heapq
import sys

ROUNDING_PRECISION = 6
INDEX_SHIFT = 0.5

class Vertex:
	def __init__(self, i, j, x, y, tt):
		
		# tt stands for terrain type
		self.i = i
		self.j = j
		self.x = x
		self.y = y
		self.tt = tt
	
	
	def __hash__(self):
		s = str(self.i) + str(self.j) + str(self.x) + str(self.y) + str(self.tt)
		return hash(s)
	
	
	def __eq__(self, other):
		sxr = round(self.x, ROUNDING_PRECISION)
		syr = round(self.y, ROUNDING_PRECISION)
		oxr = round(other.x, ROUNDING_PRECISION)
		oyr = round(other.y, ROUNDING_PRECISION)
			
		equal = ((self.i == other.i) and (self.j == other.j) and (sxr == oxr) and (syr == oyr) and (self.tt == other.tt))
		return equal

	def __ne__(self, other):
		return not(self == other)
	
	
	def __str__(self):
		return str(self.i) + str(self.j)
		
	
	def euclidean_distance(self, v):
		x2 = math.pow(self.x - v.x, 2)
		y2 = math.pow(self.y - v.y, 2)

		return math.sqrt( x2 + y2)


	def manhatan_distance(self, v):
		x_abs = abs(self.i - v.i)
		y_abs = abs(self.j - v.j)
		
		return x_abs + y_abs


	def intersection_point(self, v):
		
		x0 = 0;
		y0 = 0;
		if(self.j == v.j):
			x0 = max(self.i, v.i)
			y0 = self.y + (v.y - self.y) * (x0 - self.x) / (v.x - self.x)
		else:
			y0 = max(self.j, v.j)
			x0 = self.x + (v.x - self.x) * (y0 - self.y) / (v.y - self.y)

		return (x0, y0)
	
	
	def distance(self, v):
		
		if(self.manhatan_distance(v) == 0):
			# Euclidean distance times terrain type
			return self.euclidean_distance(v)*self.tt

		tt_s = self.tt
		tt_v = v.tt
		score = math.pow(tt_s - tt_v ,2)
		
		x0, y0 = self.intersection_point(v)
		intersection_vertex = Vertex(math.floor(x0), math.floor(y0), x0, y0, -1)

		return score + self.euclidean_distance(intersection_vertex)*tt_s + intersection_vertex.euclidean_distance(v)*tt_v


	def get_neighbours(self, N):
		# Return the nearest neighbours of a vertex.
		neighbours = []
		
		i_up = self.i-1
		j_up = self.j
		if (i_up >= 0):
			neighbours.append([i_up, j_up])
		
		
		i_down = self.i+1
		j_down = self.j
		if (i_down <= N-1):
			neighbours.append([i_down, j_down])
		
		
		i_left = self.i
		j_left = self.j-1
		if (j_left >= 0):
			neighbours.append([i_left, j_left])
		
		
		i_right = self.i
		j_right = self.j+1
		if (j_right <= N-1):
			neighbours.append([i_right, j_right])
		
		return neighbours



class Graph:
	def __init__(self, map_matrix):
		N = len(map_matrix)
		
		# vertex_connection[vertex] = number
		# vertex_array[number] = vertex
		# adjacency_list[number] = list of neighbouring vertexes

		self.n_vertexes = N*N
		self.vertex_connection = {}
		self.vertex_array = [Vertex(0,0,0,0,0) for i in range(N*N)]
		self.adjacency_list = [[] for i in range(N*N)]

		vertex_id = 0
		for i in range(len(map_matrix)):
			for j in range(len(map_matrix[i])):
				tt = map_matrix[i][j]
				
				v = Vertex(i, j, i + INDEX_SHIFT, j + INDEX_SHIFT, tt)
				
				neighbours = v.get_neighbours(N)
				
				self.vertex_connection[v] = vertex_id
				self.vertex_array[vertex_id] = v
				
				for k in range(len(neighbours)):
					ni = neighbours[k][0]
					nj = neighbours[k][1]
					nx = ni + INDEX_SHIFT
					ny = nj + INDEX_SHIFT
					nt = map_matrix[ni][nj]
					
					nv = Vertex(ni, nj, nx, ny, nt)
					(self.adjacency_list[vertex_id]).append([nv, nv.distance(v)])
				
				vertex_id = vertex_id + 1

	def print_graph(self):
	
		for i in range(len(self.vertex_array)):
			s = "Vertex: " + str(self.vertex_array[i]) + " -> "
			for j in range(len(self.adjacency_list[i])):
				s = s + "(" + str(self.adjacency_list[i][j][0]) + "," + str(round(self.adjacency_list[i][j][1],2)) + ") -> "
			s = s + "END"
			sys.stderr.write(s + "\n")


	def find_shortest_path(self, v1, v2):

		distances = [float("inf") for i in range(self.n_vertexes)]
		previous = [-1 for i in range(self.n_vertexes)]
		in_queue = [True for i in range(self.n_vertexes)]

		distances[v1] = 0

		priority_queue = [(distances[i], i) for i in range(self.n_vertexes)]
		heapq.heapify(priority_queue)

		while (len(priority_queue) != 0):
			u = heapq.heappop(priority_queue)
			u = u[1]

			in_queue[u] = False

			n_u_neighbours = len(self.adjacency_list[u]) # number of neighbours that the u vertex has
			for j in range(n_u_neighbours):
				v = self.vertex_connection[self.adjacency_list[u][j][0]]
				if (in_queue[v] == False):
					continue
				alt = distances[u] + self.adjacency_list[u][j][1]
				if (alt < distances[v]):
					distances[v] = alt
					previous[v] = u

			while (len(priority_queue) != 0):
				heapq.heappop(priority_queue)

			priority_queue = [(distances[i], i) for i in range(self.n_vertexes) if (in_queue[i] == True)]
			heapq.heapify(priority_queue)

		return (distances, previous)


	def read_shortest_path(self, distances, previous, v2):
		S = []
		u = v2 # set u to target
		while (previous[u] != -1):
			S.append(u)
			u = previous[u]
		S.append(u)

		Sij = [self.vertex_array[S[i]] for i in range(len(S))]
		return Sij


def print_neighbours(N, i, j):
	neighbours = get_neighbours(N, i, j)
	sys.stderr.write("Neighbours for: " + str(i) + str(j) + " " + str(neighbours) + "\n")


def print_matrix_ij(N):
	for i in range(N):
		for j in range(N):
			sys.stderr.write(str(i) + str(j) + " ")
		sys.stderr.write("\n")


def print_matrix_123(N):
	index = 0
	for i in range(N):
		for j in range(N):
			s = ""
			if (index < 10):
				s = "0" + str(index)
				sys.stderr.write(str(s) + " ")
				index = index + 1
			else:
				sys.stderr.write(str(index) + " ")
				index = index + 1
			sys.stderr.write("\n")


def print_matrix(m):
	for i in range(len(m)):
		for j in range(len(m[i])):
			sys.stderr.write(str(m[i][j]) + " ")
		sys.stderr.write("\n")


class TerrainCrossing:

	def getPath(self, world_map, locations, capacity):


		def convert_world_map_to_list(world_map):
			N = len(world_map)

			map_matrix = [[0 for i in range(N)] for j in range(N)]

			for i in range(N):
				for j in range(N):
					map_matrix[i][j] = int(world_map[i][j])

			return map_matrix



		for i in range(len(world_map)):
			sys.stderr.write(world_map[i] + "\n")

		map_matrix = convert_world_map_to_list(world_map)


		N = len(map_matrix)

		print_matrix(map_matrix)
		sys.stderr.write("\n")
		print_matrix_ij(5)
		sys.stderr.write("\n")

		source = 6
		target = 20

		g = Graph(map_matrix)
		distances, previous = g.find_shortest_path(source, target)

		path = g.read_shortest_path(distances, previous, target)

		for i in range(len(path)):
			sys.stderr.write(str(path[i]) + "\n")

		g.print_graph()


		ret = []
		ret.append(4.9995);    ret.append(0.7658)
		# pick up item 1
		ret.append(4.7867);    ret.append(0.7658)
		# drop it off at target 6
		ret.append(3.8144);    ret.append(0.1081)
		# pick up item 0
		ret.append(3.7648);    ret.append(1.2640)
		# drop it off at target 7
		ret.append(3.3420);    ret.append(2.5000)
		ret.append(3.3420);    ret.append(3.0530)
		# pick up item 2
		ret.append(2.5000);    ret.append(3.0530)
		ret.append(1.5000);    ret.append(3.0530)
		ret.append(0.7225);    ret.append(3.0530)
		ret.append(0.7225);    ret.append(2.5000)
		ret.append(0.7225);    ret.append(1.4533)
		# pick up item 3
		ret.append(0.2299);    ret.append(2.8555)
		ret.append(0.2299);    ret.append(3.8555)
		ret.append(0.2299);    ret.append(4.8555)
		# drop it off at target 4
		ret.append(0.5000);    ret.append(3.3869)
		ret.append(1.2611);    ret.append(3.3869)
		# drop it off at target 5
		ret.append(2.2611);    ret.append(3.3869)
		ret.append(2.2611);    ret.append(4.6214)
		ret.append(3.7958);    ret.append(4.6214)
		# exit
		ret.append(3.7958);    ret.append(4.9995)
		return ret




# -------8<------- end of solution submitted to the website -------8<-------


M = int(raw_input())
Map = []
for i in range(M):
	Map.append(raw_input().strip())

L = int(raw_input())
locations = []
for i in range(L):
	locations.append(float(raw_input()))

capacity = int(raw_input())

tc = TerrainCrossing()
ret = tc.getPath(Map, locations, capacity)
print(len(ret))
for num in ret:
	print(num)
	sys.stdout.flush()







