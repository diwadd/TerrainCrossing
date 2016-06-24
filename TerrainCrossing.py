import math

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
			print(s)
	
	
	


	

def print_neighbours(N,i,j):
	
	neighbours = get_neighbours(N,i,j)
	print("Neighbours for: " + str(i) + str(j) + " " + str(neighbours))
	

def print_matrix_ij(N):
	for i in range(N):
		for j in range(N):
			print(str(i) + str(j) + " ", end="")
		print()


def print_matrix(m):
	for i in range(len(m)):
		for j in range(len(m[i])):
			print(str(m[i][j]) + " ", end="")
		print()


map_matrix = [[9,9,9,9,1],[9,0,9,9,1],[9,0,0,9,1],[9,9,0,9,1],[0,1,2,3,4]]
N = len(map_matrix)

print_matrix(map_matrix)
print_matrix_ij(5)

g = Graph(map_matrix)
g.print_graph()


