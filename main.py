#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
from meshparty import meshwork, trimesh_vtk

# skel_dir = 'data'
skel_dir = "C:\\Users\\hoda.aboueich\\Downloads\\hoda_skeletons\\hoda_skeletons\\data" # data file

# The files in the directory `data` contain skeletons of all of the basket cells in a central column of primary visual cortex.
# The number is the so-called `root id` of a given cell at a given state of proofreading, and the suffix `.h5` indicates that it's an "hdf5" file, which is a generic file format that can hold all sorts of data. In this case, these h5 files are formatted as a "meshwork" object, which holds information about the morphology and connectivity of neurons.

import glob

# Glob lets you do shell-like file lists in python
# test_str1 = f"{skel_dir}/*.h5"
# print(test_str1)
files = glob.glob(f"{skel_dir}/*.h5")

# This command loads a neuron file into python

# Define constants
NUMFILES = int(59) # This is the number of neuron files
SUM1TONUMFILES = int(NUMFILES*(NUMFILES-1)/2) # This is the sum of numbers of 1 and the number of files; this represents the number of unique dendrite-to-dendrite interactions
NUMCOLUMNS = int(9) # This is the number of data columns output by the program
NUMROWS = SUM1TONUMFILES

r, c = NUMCOLUMNS, NUMROWS
DataMatrix = [[0 for x in range(r)] for y in range(c)]

print("NUMFILES = ", NUMFILES, " SUM1TONUMFILES = ", SUM1TONUMFILES, " NUMCOLUMNS = ", NUMCOLUMNS, " NUMROWS = ", NUMROWS)

z = 0

j = 0
while j < NUMFILES-1:
    DataMatrix[z][0] = j

    nrn = meshwork.load_meshwork(files[j])

# Since we are interested in dendrite-dendrite interactions, the first thing we want to do is limit our analysis to dendrites.
# To show you how this goes, we're going to plot the cell along the way. We've been using vtk as a quick way to plot cells in 3d.

# The following two lines of code will bring up a 3d plot of the neuron, with the dendrite in red and the axon in blue.
# You'll have to close the window to continue.

    ###ska = trimesh_vtk.skeleton_actor(
       ### nrn.skeleton, line_width=2, vertex_data=nrn.anno.is_axon.skel_mask
    ###)

    ###trimesh_vtk.render_actors([ska])

# You can see all of the vertices of the object like so:

    nrn.skeleton.vertices

# How many vertices are there?

### len(nrn.skeleton.vertices)

    DataMatrix[z][2] = len(nrn.skeleton.vertices)

    nlength = DataMatrix[z][2]

# But now we want to mask out the axon and leave only the dendrite.
# When I made these skeletons to begin with, I already labeled the axon and dendrite (which is how we did the colors before), so we can just use that.
# However, these objects contain both skeletons (which are trees) and "meshes", which are somewhat more fleshed out.
# In this case, these aren't true meshes, but they still capture a bit more structure of the cells.
#
# Meshwork objects have meshes, skeletons, and annotations ( like axon/dendrite labels or like synapses) that are all linked together into a core object.
# One key operation is masking, which limits an object to a given subset.
# The function `apply_mask` takes a boolean array on the mesh, which we can get similar to the colors above.

    nrn.apply_mask(~nrn.anno.is_axon.mesh_mask)

# After applying a mask, you can see the effect by plotting the object again.

    ###ska = trimesh_vtk.skeleton_actor(
        ###nrn.skeleton, line_width=2, vertex_data=nrn.anno.is_axon.skel_mask.astype(int)
    ###)

    ###trimesh_vtk.render_actors([ska])

# This time, you should see just a red part, and the skeleton will be much smaller.
#
# You can reset the mask to the original state with `nrn.reset_mask()`
#
# Now how many vertices are there?

### len(nrn.skeleton.vertices)

    DataMatrix[z][3] = len(nrn.skeleton.vertices)

    nlength2 = DataMatrix[z][3]

# How do we measure the distance between dendrites?
#
# As a first approximation, we can measure the closest distances between the dendrite vertices.
#
# Scipy has an object called a K-d tree that does this very efficiently.
# Let's use it here.
#
# We're going to use the mesh vertices instead of the skeleton vertices to get a slightly richer representation.

    i = j + 1
    while i < NUMFILES:
        print("i = ", i, " z = ", z)
        nrn_other = meshwork.load_meshwork(files[i])
        nrn_other.apply_mask(~nrn_other.anno.is_axon.mesh_mask)

        DataMatrix[z][0] = j
        DataMatrix[z][1] = i
        DataMatrix[z][2] = nlength
        DataMatrix[z][3] = nlength2
        i += 1

        from scipy import spatial

        kdt = spatial.KDTree(nrn_other.mesh.vertices)
        ds, inds = kdt.query(
        nrn.mesh.vertices)  # This returns the distance and closest index in the tree object for each of the vertices queried.
        len(ds)

# What is the closest distance between the two dendrites?

        # print(np.min(ds) / 1000)  # The vertex units are in nanometers, so let's divide by 1000 to bring into microns
        DataMatrix[z][4] = (np.min(ds) / 1000)
        DataMatrix[z][5] = np.argmin(ds)  # Returns the index of the minimum value along an axis

        DataMatrix[z][6] = nrn.seg_id

        DataMatrix[z][7] = nrn_other.seg_id

        DataMatrix[z][8] = nrn.mesh.vertices[np.argmin(ds)] / [4, 4, 40]

        z += 1
    j += 1
# Now, how would you go about getting this minimum distance for all pairs of cells?

# Copy DataMatrix into Pandas DataFrame
df = pd.DataFrame(DataMatrix, columns=['Column_A', 'Column_B', 'Column_C', 'Column_D', 'Column_E', 'Column_F', 'Column_G', 'Column_H', 'Column_I'])
pd.options.display.max_rows = None
pd.options.display.max_columns = None


print(df)
print(type(df))

# Print DataMatrix rows 1-59
k = 0
while k < SUM1TONUMFILES:
    print(DataMatrix[k][0], " ", DataMatrix[k][1], " ", DataMatrix[k][2], " ", DataMatrix[k][3], " ", DataMatrix[k][4], " ", DataMatrix[k][5], " ", DataMatrix[k][6], " ", DataMatrix[k][7], " ", DataMatrix[k][8])
    k += 1

a = 0
count0_1 = 0
count1_2 = 0
count2_3 = 0
count3_4 = 0
count4_5 = 0
count_greater_than_5 = 0

# Create sublists
r, c = NUMCOLUMNS, NUMROWS
DataMatrix0_1 = [[0 for x in range(r)] for y in range(c)]
DataMatrix1_2 = [[0 for x in range(r)] for y in range(c)]
DataMatrix2_3 = [[0 for x in range(r)] for y in range(c)]
DataMatrix3_4 = [[0 for x in range(r)] for y in range(c)]
DataMatrix4_5 = [[0 for x in range(r)] for y in range(c)]
DataMatrix_greater_than_5 = [[0 for x in range(r)] for y in range(c)]

b = 0
c = 0
d = 0
e = 0
f = 0
g = 0

while a < SUM1TONUMFILES:
    if DataMatrix[a][4] > 0 and DataMatrix[a][4] <= 1:
        count0_1 += 1
        DataMatrix0_1[b] = DataMatrix[a]
        b += 1
    elif DataMatrix[a][4] <= 2:
        count1_2 += 1
        DataMatrix1_2[c] = DataMatrix[a]
        c += 1
    elif DataMatrix[a][4] <= 3:
        count2_3 += 1
        DataMatrix2_3[d] = DataMatrix[a]
        d += 1
    elif DataMatrix[a][4] <= 4:
        count3_4 += 1
        DataMatrix3_4[e] = DataMatrix[a]
        e += 1
    elif DataMatrix[a][4] <= 5:
        count4_5 += 1
        DataMatrix4_5[f] = DataMatrix[a]
        f += 1
    elif DataMatrix[a][4] > 5:
        count_greater_than_5 += 1
        DataMatrix_greater_than_5[g] = DataMatrix[a]
        g += 1

    a += 1

print("Between 0 and 1 ", count0_1)
print("Between 1 and 2 ", count1_2)
print("Between 2 and 3 ", count2_3)
print("Between 3 and 4 ", count3_4)
print("Between 4 and 5 ", count4_5)
print("Greater than 5", count_greater_than_5)

# Plot the data
import matplotlib.pyplot as plt

x = ['0-1', '1-2', '2-3', '3-4', '4-5', '>5']
y = [count0_1, count1_2, count2_3, count3_4, count4_5, count_greater_than_5]

plt.scatter(x, y)
plt.show()

print("Printing DataMatrix0_1, distances between 0 and 1")
k = 0
while k < b:
    print(DataMatrix0_1[k][0], " ", DataMatrix0_1[k][1], " ", DataMatrix0_1[k][2], " ", DataMatrix0_1[k][3], " ", DataMatrix0_1[k][4], " ", DataMatrix0_1[k][5], " ", DataMatrix0_1[k][6], " ", DataMatrix0_1[k][7], " ", DataMatrix0_1[k][8])
    k += 1
print("Printing DataMatrix1_2, distances between 1 and 2")
k = 0
while k < c:
    print(DataMatrix1_2[k][0], " ", DataMatrix1_2[k][1], " ", DataMatrix1_2[k][2], " ", DataMatrix1_2[k][3], " ", DataMatrix1_2[k][4], " ", DataMatrix1_2[k][5], " ", DataMatrix1_2[k][6], " ", DataMatrix1_2[k][7], " ", DataMatrix1_2[k][8])
    k += 1
print("Printing DataMatrix2_3, distances between 2 and 3")
k = 0
while k < d:
    print(DataMatrix2_3[k][0], " ", DataMatrix2_3[k][1], " ", DataMatrix2_3[k][2], " ", DataMatrix2_3[k][3], " ", DataMatrix2_3[k][4], " ", DataMatrix2_3[k][5], " ", DataMatrix2_3[k][6], " ", DataMatrix2_3[k][7], " ", DataMatrix2_3[k][8])
    k += 1
print("Printing DataMatrix3_4, distances between 3 and 4")
k = 0
while k < e:
    print(DataMatrix3_4[k][0], " ", DataMatrix3_4[k][1], " ", DataMatrix3_4[k][2], " ", DataMatrix3_4[k][3], " ", DataMatrix3_4[k][4], " ", DataMatrix3_4[k][5], " ", DataMatrix3_4[k][6], " ", DataMatrix3_4[k][7], " ", DataMatrix3_4[k][8])
    k += 1
print("Printing DataMatrix4_5, distances between 4 and 5")
k = 0
while k < f:
    print(DataMatrix4_5[k][0], " ", DataMatrix4_5[k][1], " ", DataMatrix4_5[k][2], " ", DataMatrix4_5[k][3], " ", DataMatrix4_5[k][4], " ", DataMatrix4_5[k][5], " ", DataMatrix4_5[k][6], " ", DataMatrix4_5[k][7], " ", DataMatrix4_5[k][8])
    k += 1
print("Printing DataMatrix_greater_than_5, distances greater than 5")
k = 0
while k < g:
    print(DataMatrix_greater_than_5[k][0], " ", DataMatrix_greater_than_5[k][1], " ", DataMatrix_greater_than_5[k][2], " ", DataMatrix_greater_than_5[k][3], " ", DataMatrix_greater_than_5[k][4], " ", DataMatrix_greater_than_5[k][5], " ", DataMatrix_greater_than_5[k][6], " ", DataMatrix_greater_than_5[k][7], " ", DataMatrix_greater_than_5[k][8])
    k += 1

df0_1 = pd.DataFrame(DataMatrix0_1, columns=['Column_A', 'Column_B', 'Column_C', 'Column_D', 'Column_E', 'Column_F', 'Column_G', 'Column_H', 'Column_I'])
df1_2 = pd.DataFrame(DataMatrix1_2, columns=['Column_A', 'Column_B', 'Column_C', 'Column_D', 'Column_E', 'Column_F', 'Column_G', 'Column_H', 'Column_I'])
df2_3 = pd.DataFrame(DataMatrix2_3, columns=['Column_A', 'Column_B', 'Column_C', 'Column_D', 'Column_E', 'Column_F', 'Column_G', 'Column_H', 'Column_I'])
df3_4 = pd.DataFrame(DataMatrix3_4, columns=['Column_A', 'Column_B', 'Column_C', 'Column_D', 'Column_E', 'Column_F', 'Column_G', 'Column_H', 'Column_I'])
df4_5 = pd.DataFrame(DataMatrix4_5, columns=['Column_A', 'Column_B', 'Column_C', 'Column_D', 'Column_E', 'Column_F', 'Column_G', 'Column_H', 'Column_I'])
df_greater_than_5 = pd.DataFrame(DataMatrix_greater_than_5, columns=['Column_A', 'Column_B', 'Column_C', 'Column_D', 'Column_E', 'Column_F', 'Column_G', 'Column_H', 'Column_I'])
