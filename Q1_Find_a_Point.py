# imports statements (math for math operations)
import math

# opens input file and output file
input_file = open("find.in" )
output_file = open("find.out", "w")

# accepts and type-casts inputs
n, m, dx, dy, k = [float(i) for i in input_file.readline().split()]

# iterates through next k lines for coordinates
for i in range(0, int(k)):

	# accept and type-cast coordinates
	x, y = [float(j) for j in input_file.readline().split()]

	# uses smaller if distance between integer values is same
	# otherwise, rounds to nearest integer
	res_x = math.floor(x / dx) if (abs(math.ceil(x / dx) - x / dx)) == (abs(math.floor(x / dx) - x / dx)) \
		else round(x / dx)
	res_y = math.floor(y / dy) if (abs(math.ceil(y / dy) - y / dy)) == (abs(math.floor(y / dy) - y / dy)) \
		else round(y / dy)

	# uses farthest point on grid if off the grid
	if res_x >= n: res_x = n - 1
	if res_y >= m: res_y = m - 1

	# writes integer coordinates in output file
	output_file.write(str(int(res_x)) + " " + str(int(res_y)) + "\n")
