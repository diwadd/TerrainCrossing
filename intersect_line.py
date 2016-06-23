import math

class Location:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.xi = int(math.floor(x))
        self.yi = int(math.floor(y))


def distance_euclidean(l1, l2):
    x2 = math.pow(l1.x - l2.x,2)
    y2 = math.pow(l1.y - l2.y,2)

    return math.sqrt(x2 + y2)

def distance_manhatan(l1, l2):
    x_abs = abs(l1.xi - l2.xi)
    y_abs = abs(l1.yi - l2.yi)
    
    return x_abs + y_abs
    

def interection_point(l1, l2):
    
    if(distance_manhatan(l1,l2) == 0):
        return distance_euclidean(l1,l2)

    x0 = 0;
    y0 = 0;
    if(l1.yi == l2.yi):
        x0 = max(l1.xi, l2.xi)
        y0 = l1.y + (l2.y - l1.y) * (x0 - l1.x) / (l2.x - l1.x)
    else:
        y0 = max(l1.yi, l2.yi)
        x0 = l1.x + (l2.x - l1.x) * (y0 - l1.y) / (l2.y - l1.y)

    return (x0, y0)

# l1 = Location(1.25, 1.25)
# l2 = Location(2.25, 1.75)

l1 = Location(0, 0)
l2 = Location(50, 50)

de = distance_euclidean(l1, l2)
dm = distance_manhatan(l1, l2)

print("de: " + str(de))
print("dm: " + str(dm))

ip = interection_point(l1, l2)

print("ip: " + str(ip))


