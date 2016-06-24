import math

ROUNDING_PRECISION = 6

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
	
	
	def euclidean_distance(v1, v2):
		x2 = math.pow(v1.x - v2.x, 2)
		y2 = math.pow(v1.y - v2.y, 2)

		return math.sqrt( x2 + y2)


	def manhatan_distance(v1, v2):
		x_abs = abs(v1.i - v2.i)
		y_abs = abs(v1.j - v2.j)
		
		return x_abs + y_abs


	def intersection_point(v1, v2):
		
		x0 = 0;
		y0 = 0;
		if(v1.j == v2.j):
			x0 = max(v1.i, v2.i)
			y0 = v1.y + (v2.y - v1.y) * (x0 - v1.x) / (v2.x - v1.x)
		else:
			y0 = max(v1.j, v2.j)
			x0 = v1.x + (v2.x - v1.x) * (y0 - v1.y) / (v2.y - v1.y)

		return (x0, y0)






class Graph:
	def __init__(self, map_matrix):
		N = len(map_matrix)
		self.vertex_connection = [[0, 0] for i in range(N*N)]
		self.adjacency_list = []

		for i in range(len(map_matrix)):
			for j in range(len(map_matrix[i])):
				get_neighbours(i,j)
			

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


print(v1.euclidean_distance(v2))
print(v1.intersection_point(v2))


