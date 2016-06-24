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



class Graph:
	def __init__(self, map_matrix):
		N = len(map_matrix)
		self.vertex_connection = {}
		self.adjacency_list = [[] for i in range(N*N)]

		vertex_id = 0
		for i in range(len(map_matrix)):
			for j in range(len(map_matrix[i])):
				tt = map_matrix[i][j]
				
				neighbours = get_neighbours(N,i,j)
				
				v = Vertex(i, j, i + INDEX_SHIFT, j + INDEX_SHIFT, tt)
				
			

def get_neighbours(N,i,j):
	neighbours = []
	
	i_up = i-1
	j_up = j
	if (i_up >= 0):
		neighbours.append([i_up, j_up])
	
	
	i_down = i+1
	j_down = j
	if (i_down <= N-1):
		neighbours.append([i_down, j_down])
	
	
	i_left = i
	j_left = j-1
	if (j_left >= 0):
		neighbours.append([i_left, j_left])
	
	
	i_right = i
	j_right = j+1
	if (j_right <= N-1):
		neighbours.append([i_right, j_right])
	
	return neighbours
	

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


map_matrix = [[9,9,9,9],[9,0,9,9],[9,0,0,9],[9,9,0,9]]
N = len(map_matrix)

M = 5
print_matrix_ij(M)

print_neighbours(M, 0, 0)
print_neighbours(M, 1, 1)
print_neighbours(M, 2, 3)
print_neighbours(M, 0, 0)
print_neighbours(M, 4, 4)


v1 = Vertex(0, 0, 0.5, 0.5, 9)
v2 = Vertex(1, 1, 1.5, 1.5, 9)

v1 = Vertex(1, 1, 1.5, 1.5, 1)
v2 = Vertex(0, 0, 0.5, 0.5, 9)

print(v1.euclidean_distance(v2))
print(v1.intersection_point(v2))
print("Distance: " + str(v1.distance(v2)))

