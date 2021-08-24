import matplotlib as mt
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from numpy.ma.core import sqrt


#to apply boundary condition to numerator of temperature of a cell
def bc_n(boundary,i,j,wt,x,y):
    n = 0 
    if boundary['type'] == 1:
        n =  wt[i,j,boundary['wt_no']]*boundary['t']*2
    elif boundary['type'] == 2:
        n = int(0)
    elif boundary['type'] == 3:
        dx = x[j + 1] - x[j]
        dy = y[i + 1] - y[i]
        if boundary['name'] in ['left','right']:
            d = dx
        else: 
            d = dy
        n =  wt[i,j,boundary['wt_no']]*d*boundary['t']
    elif boundary['type'] == 4:
        dx = x[j + 1] - x[j]
        dy = y[i + 1] - y[i]
        c = boundary['c']
        To = boundary['t']
        if boundary['name'] in ['left','right']:
            d = dx
        else: 
            d = dy
        phi = (c * To * (-1)) / ((1 / d) - (c / 2))
        n = wt[i,j,boundary['wt_no']]*phi
    return n


#to apply boundary condition to denominator of temperature of a cell
def bc_d(boundary,i,j,wt,x,y):
    r = 0
    if boundary['type'] == 1:
        r = wt[i,j,boundary['wt_no']]*(-1)
    elif boundary['type'] == 2:
        r = wt[i,j,boundary['wt_no']]
    elif boundary['type'] == 3:
        r = wt[i,j,boundary['wt_no']]
    elif boundary['type'] == 4:
        dx = x[j + 1] - x[j]
        dy = y[i + 1] - y[i]
        c = boundary['c']
        if boundary['name'] in ['left','right']:
            d = dx
        else: 
            d = dy
        epsilon = ((1 / d) + (c / 2)) / ((1 / d) - (c / 2))
        r = wt[i,j,boundary['wt_no']]*epsilon
    return r


# Base parameter required for the simulation
l = float(input("length: "))
w = float(input("width: "))
rx = float(input("successive ratio in x direction: "))
ry = float(input("successive ratio in y direction: "))

# Number of cells
cx = int(input("cells in x direc:"))
cy = int(input("cells in y direc:"))

# calculating the length of first cell using GP
if rx == 1:
    fclx = l/cx
else:
    fclx = l * abs((1 - rx)/(1 - rx**cx))

if ry == 1:
    fcly = w/cy
else:
    fcly = w * abs((1 - ry)/(1 - ry**cy))

# Determining x coordinates and y coordiantes of the TempM including boundary
x_coo = np.zeros(1)
t = 0
for i in range(cx):
    t += fclx*(rx**i)
    x_coo = np.append(x_coo, t)
x_cells = cx

y_coo = np.zeros(1) 
t = 0
for i in range(cy):
    t += fcly*(ry**i)
    y_coo = np.append(y_coo, t)
y_cells = cy

# plotting the TempM using x and y coordinates
for i in x_coo:
    plt.vlines(i, ymin=0, ymax= w)
for j in y_coo:
    plt.hlines(j, xmin=0, xmax= l)
plt.show()

# Matrix of the mesh
Mesh = np.empty(shape=(y_cells + 1, x_cells + 1))

# boundary conditions
print("Enter 1 for dirichiletch \n Enter 2 for Neumann Homogenous \n Enter 3 for Neumann Heterogenous \n Enter 4 for Robins")
lb = {
    'name' : 'left',
    'cells' : [],
    'wt_no' : 3
}
rb = {
    'name' : 'right',
    'cells' : [],
    'wt_no' : 2
}
tb = {
    'name' : 'top',
    'cells' : [],
    'wt_no' : 0
}
bb = {
    'name' : 'bottom',
    'cells' : [],
    'wt_no' : 1
}
lb['type'] = int(input("left boundary type:"))
rb['type'] = int(input("right boundary type:"))
tb['type'] = int(input("top boundary type:"))
bb['type'] = int(input("bottom boundary type:"))
for i in [lb, rb, tb, bb]:
    if i['type'] == 1:
        i['t'] = float(input("Enter temperature for " + i['name'] + " boundary: "))
    elif i['type'] == 2:
        pass
    elif i['type'] == 3:
        i['t'] = float(input("Please enter value of heat flux for " + i['name'] + " boundary: "))
    elif i['type'] == 4:
        print("equation : dT/dx = c(T- To)")
        i['c'] = float(input("Please enter value of c for " + i['name'] + " boundary: "))
        i['t'] = float(input("Please enter value of To for " + i['name'] + " boundary: "))
    else:
        print("The input for " + i['name'] + ' boundary is invalid.')
        exit()
        
# real boundary
rows = [i for i in range(0, y_cells)]
columns = [i for i in range(0, x_cells)]

# Initialising a cell center temperaure matrix
TempM = np.empty(shape=(y_cells, x_cells))

# Initialising the value of cell centers in the matrix
for i in rows:
    for j in columns:
        TempM[i, j] = 0

for i in rows:
    for j in columns:
        if j == x_cells - 1 and i in range(0, y_cells): #east
            rb['cells'].append((i,j))
        if i == 0 and j in range(0, x_cells): # south
            bb['cells'].append((i,j))
        if i == y_cells - 1 and j in range(0, x_cells): # north
            tb['cells'].append((i,j))
        if j == 0 and i in range(0, y_cells): # west
            lb['cells'].append((i,j))

#Calculating weights
wt = np.zeros((y_cells, x_cells, 5))

for i in rows:
    for j in columns:
        if (i >= 1 and i <= y_cells - 2) and (j >= 1 and j <= x_cells - 2):
         # matrix excluding boundary cells
            dxw = (x_coo[j + 1] - x_coo[j - 1])/2
            dxe = (x_coo[j + 2] - x_coo[j])/2
            dyn = (y_coo[i + 2] - y_coo[i])/2
            dys = (y_coo[i + 1] - y_coo[i - 1])/2
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
# boundary cells
        elif j == x_cells - 1 and i in range(1, y_cells - 2):#east
            dxw = (x_coo[j + 1] - x_coo[j - 1])/2
            dxe = (x_coo[j + 1] - x_coo[j])
            dyn = (y_coo[i + 2] - y_coo[i])/2
            dys = (y_coo[i + 1] - y_coo[i - 1])/2
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        elif i == 0 and j in range(1, x_cells - 2):
            dxw = (x_coo[j + 1] - x_coo[j - 1])/2
            dxe = (x_coo[j + 2] - x_coo[j])/2
            dyn = (y_coo[i + 2] - y_coo[i])/2
            dys = (y_coo[i + 1] - y_coo[i])
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        elif i == y_cells-1 and j in range(1, x_cells - 2):
            dxw = (x_coo[j + 1] - x_coo[j - 1])/2
            dxe = (x_coo[j + 2] - x_coo[j])/2
            dyn = (y_coo[i + 1] - y_coo[i])
            dys = (y_coo[i + 1] - y_coo[i - 1])/2
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        elif j == 0 and i in range(1, y_cells - 1):
            dxw = (x_coo[j + 1] - x_coo[j])
            dxe = (x_coo[j + 2] - x_coo[j])/2
            dyn = (y_coo[i + 2] - y_coo[i])/2
            dys = (y_coo[i + 1] - y_coo[i - 1])/2
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
# corners
        elif (i, j) == (0, 0):
            dxw = (x_coo[j + 1] - x_coo[j])
            dxe = (x_coo[j + 2] - x_coo[j])/2
            dyn = (y_coo[i + 2] - y_coo[i])/2
            dys = (y_coo[i + 1] - y_coo[i])
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        elif (i, j) == (y_cells - 1, x_cells - 1):
            dxw = (x_coo[j + 1] - x_coo[j - 1])/2
            dxe = (x_coo[j + 1] - x_coo[j])
            dyn = (y_coo[i + 1] - y_coo[i])
            dys = (y_coo[i + 1] - y_coo[i - 1])/2
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        elif (i, j) == (0, x_cells-1):
            dxw = (x_coo[j + 1] - x_coo[j - 1])/2
            dxe = (x_coo[j + 1] - x_coo[j])
            dyn = (y_coo[i + 2] - y_coo[i])/2
            dys = (y_coo[i + 1] - y_coo[i])
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        elif (i, j) == (y_cells-1, 0):
            dxw = (x_coo[j + 1] - x_coo[j])
            dxe = (x_coo[j + 2] - x_coo[j])/2
            dyn = (y_coo[i + 1] - y_coo[i])
            dys = (y_coo[i + 1] - y_coo[i - 1])/2
            dx = x_coo[j + 1] - x_coo[j]
            dy = y_coo[i + 1] - y_coo[i]
            a_e = dy/dxe
            a_n = dx/dyn
            a_w = dy/dxw
            a_s = dy/dys
            a_p = (a_e + a_w + a_n + a_s) * (-1)
        # wt are saved in the format North-South-East-West-Summation
        wt[i, j, 0] = a_n
        wt[i, j, 1] = a_s
        wt[i, j, 2] = a_e
        wt[i, j, 3] = a_w
        wt[i, j, 4] = a_p

#iterations
iterations = int(input("iterations: "))
for k in range(iterations):
    for i in rows:
        for j in columns:
            n = float(0)
            d = float(wt[i,j,4])
            #west
            if (i,j) in lb['cells']:
                n += float(bc_n(lb,i,j,wt,x_coo,y_coo))
                d += float(bc_d(lb,i,j,wt,x_coo,y_coo))
            else:
                n += float(wt[i,j,3]*TempM[i,j-1])
            #east
            if (i,j) in rb['cells']:
                n += float(bc_n(rb,i,j,wt,x_coo,y_coo))
                d += float(bc_d(rb,i,j,wt,x_coo,y_coo))
            else: 
                n += float(wt[i,j,2]*TempM[i,j+1])
            #north
            if (i,j) in tb['cells']:
                n += float(bc_n(tb,i,j,wt,x_coo,y_coo))
                d += float(bc_d(tb,i,j,wt,x_coo,y_coo))
            else:
                n += float(wt[i, j, 0]*TempM[i+1,j])
            #south
            if (i,j) in bb['cells']:
                n += float(bc_n(bb,i,j,wt,x_coo,y_coo))
                d += float(bc_d(bb,i,j,wt,x_coo,y_coo))
            else:
                n += float(wt[i,j,1]*TempM[i-1,j])

            TempM[i, j] = (-1)*n/d

            if k == iterations - 2 :
                TempM_old = np.copy(TempM)

print(TempM, wt)

TempM_old = np.multiply(TempM_old,TempM_old)
TempM_new = np.multiply(TempM, TempM)
sum_new = np.sum(TempM_new)/float((i*j))
sum_old = np.sum(TempM_old)/float((i*j))
rms_old = float(sqrt(sum_old))
rms_new = float(sqrt(sum_new))
residual = float(rms_new - rms_old)
print(residual)
##---------------CFD CODE WORKS FOR SQUARE/RECTANGLE UNIFORM/NON-UNIFORM MESHES------------------##

#Variation of temperature in y direction for given value of x
tpx = np.empty(shape=y_cells)
y = np.empty(shape=y_cells)
cx = int(input("cell number in x direction to plot a temperature graph of: "))
for i in range(y_cells):
    tpx[i] = TempM[i][cx]
    y[i] = (y_coo[i + 1] + y_coo[i])/2
plt.plot(y,tpx)
plt.show()

#Variation of temperature in y direction for given value of x
tpy = np.empty(shape=x_cells)
x = np.empty(shape=x_cells)
cy = int(input("cell number in y direction to plot a temperature graph of: "))
for i in range(x_cells):
    tpy[i] = TempM[cy][i]
    x[i] = (x_coo[i + 1] + x_coo[i])/2
plt.plot(x,tpy)
plt.show()

##---------------3d graph and contour works for squares only------------------## 
#Plotting a contour
xcc = np.empty(shape=(y_cells, x_cells))
ycc = np.empty(shape=(y_cells, x_cells))

for i in range(x_cells):
    for j in range(y_cells):
        xcc[i][j] = (x_coo[i+1] + x_coo[i])/2
for i in range(x_cells):
    for j in range(y_cells):
        ycc[i][j] = (y_coo[j+1] + y_coo[j])/2

plt.contourf(ycc, xcc,TempM,colours = 'black')
plt.contourf(ycc, xcc,TempM,40,cmap=plt.cm.jet)
plt.show()

#plotting 3d surface graph
fig = plt.figure()
axes = fig.gca(projection ='3d')
axes.plot_surface(ycc, xcc, TempM, cmap=plt.cm.jet)
  
plt.show()

###  IMPORTANT :check code for robins condition again  ###