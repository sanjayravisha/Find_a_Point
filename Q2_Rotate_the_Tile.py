# imports statements (math for math operations)
import math

# opens input file and output file
input_file = open("rotate.in")
output_file = open("rotate.out", "w")

# accepts and type-casts inputs
h, v, l, x, y = [float(i) for i in input_file.readline().split()]

# calculates altitude and width of image using formulas
a = l / (2 * math.tan(math.radians(v / 2)))
w = 2 * a * (math.tan(math.radians(h / 2)))

# calculates pitch and roll using formulas
# positive pitch if positive x and negative pitch if negative x
# negative roll if positive y and positive roll if negative y
pitch = math.degrees(math.atan(abs(x) / a)) if x >= 0 else -math.degrees(math.atan(abs(x) / a))
roll = -math.degrees(math.atan(abs(y) / a)) if y >= 0 else math.degrees(math.atan(abs(y) / a))

# writes altitude, width, roll, and pitch in output file
output_file.write(str(round(a, 1)) + " " + str(round(w, 1)) + " " +
                  str(round(roll, 1)) + " " + str(round(pitch, 1)))