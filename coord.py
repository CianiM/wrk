import numpy as np


Nx,Ny=69,49
Np=Nx*Ny
HEAD=2
f=open('grid_coord.dat','r')
wx=open('x_coord.dat','w')
wy=open('y_coord.dat','w')

dat=f.readlines()[HEAD:]
list_coord=[word for line in dat for word in line.split()]
#coord=[float(i) for i in list_coord]
#x=np.array(coord)
x_c=list_coord[:Nx]
for c in x_c:
    x_coord=str(c)
    wx.write(x_coord+'\n')
y_c=list_coord[Np:]
for j in np.arange(0,Ny):
    y_coord=str(y_c[j*Nx])
    wy.write(y_coord+'\n')


#for j in np.arange(0,Ny):
#    n=j*Np
#    y_coord=str(y_c[n])
#    wy.write(y_coord+'\n')
print('done')
