import mcschematic
import os

schem = mcschematic.MCSchematic()

# read a txt file containing x, y, z coordinates and create a schematic with stone blocks at those coordinates

# open file "output.txt" for reading
f = open("output.txt", "r")

# read every line in the file
for line in f:
    # coordinates are separated by a comma and a space with \n at the end of each line
    # split the line into a list of strings
    coords = line.split(", ")
    # convert the strings to integers
    x = int(coords[0])
    y = int(coords[1])
    z = int(coords[2])
    # set the block at the coordinates to stone
    schem.setBlock((x, y, z), "minecraft:emerald_block")

# close the file
f.close()

if not os.path.exists("myschems"):
    os.mkdir("myschems")

# save the schematic
schem.save("myschems", "univ-teatru", mcschematic.Version.JE_1_20_1)
