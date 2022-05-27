import numpy as np 
from scipy.interpolate import interp1d 
import matplotlib.pyplot as plt
import sys
R=287.05
gamma=1.4
T=298
a = np.sqrt(R*T*gamma)
U = 34.6
d99 = 0.08
lchar = 0.4*d99
Re = 9.36e5
#x0,Nx,Ny = 0, 103, 28
def main(x0,Nx,Ny):
    with open('Input_u.dat','r') as f:
        dat = f.readlines()[3:]
        list_dat =  [ word for line in dat for word in line.split()]
        yu = [ float(s) for s in list_dat]
        yu = np.array(yu)
        yu = yu.reshape((101,2))
        yu[:,1]=yu[:,1]*a/U
    #    fit = np.polyfit(yu[:,0],yu[:,1],1)
    #    print(fit)
    #line = np.poly1d(fit)
        x1 = yu[:,0]
    #f = interp1d(x, yu[0,1], kind = 'cubic')
    '''
    with open('Input_up.dat','r') as f:
        dat = f.readlines()[3:]
        list_dat =  [ word for line in dat for word in line.split()]
        yup = [ float(s) for s in list_dat]
        yup = np.array(yup)
        yup = yup.reshape((42,2))
        yup[:,1]=yup[:,1]*Re*lchar
        x2 = yup[:,0]
    '''
    with open('uprime_hw.dat','r') as f:
        dat = f.readlines()[1:]
        list_dat =  [ word for line in dat for word in line.split()]
        yup = [ float(s) for s in list_dat]
        yup = np.array(yup)
        dim = len(yup)/3
        yup = yup.reshape((int(dim),3))
        yup[:,2]=yup[:,2]*Re*lchar
        kk = np.zeros((int(dim)+1,3))
        kk[1:,:] = yup[:,:]
        yup = kk
        x2 = yup[:,0]/420

    with open('y_coord.dat','r') as f:
        y_coord = f.readlines()
        y_coord = [float(s) for s in y_coord]
        y_coord = np.array(y_coord)
        y_coord = y_coord.reshape((Ny,Nx))
        y_coord = y_coord.T
        yy = y_coord[x0:Nx,:]
        n1 = 0
        n2 = 0
        y_dat = yy[0,:]
        for i in y_dat:
            if i <= 0.15:
                #print(i)
                n1+=1
        for i in y_dat:
            if i <= 0.12:
                n2+=1
    #print(n)
        f1 = interp1d(x1,yu[:,1])
        f2 = interp1d(x2,yup[:,1])
        f2 = interp1d(x2,yup[:,2])
    interp1 = f1(y_dat[0:n1])
    interp2 = f2(y_dat[0:n2])

#plt.plot(x,yu[0:,1],'-r',y[0:n],f(y[0:n]),'--')
#plt.show()

    vel = np.zeros(Ny)
    nut = np.zeros(Ny)
    vel[0:n1] = interp1
    nut[0:n2] = interp2
    for i in np.arange(n1,Ny):
        vel[i] = interp1[-1]
    for i in np.arange(n2,Ny):    
        #c1 = interp2[-2]/(1-np.exp(-yy[x0,n2]+yy[x0,-1]))
        c1 = interp2[-2]/(np.tanh(-yy[x0,n2]+yy[x0,-1]))**10
        #nut[i] = interp2[-1] #- interp2[-1]/dy*yy[x0,i]
        #nut[i] = c1*(1-np.exp(-yy[x0,i]+yy[x0,-1]))
        nut[i] = c1*(np.tanh(-yy[x0,i]+yy[x0,-1]))**10  

#print(y[1])
#plt.plot(y,vel)
#plt.show()
    with open('x_coord.dat','r') as f:
        x_coord = f.readlines()
        x_coord = [float(s) for s in x_coord]
        x_coord = np.array(x_coord)
        x_coord = x_coord.reshape((Ny,Nx))
        x_coord = x_coord.T
    xx = x_coord[x0:Nx,:]
    #print(len(xx[:,0]))
    x_wall = xx[:,0]
    y_wall = yy[:,0]
    r_wall = x_wall**2.0+y_wall**2.0
    r      = np.zeros((Nx-x0,Ny))
    r[:,:] = xx[:,:]**2.0+yy[:,:]**2.0
    mindist = np.zeros((Nx-x0,Ny))
    temp = np.zeros(Nx-x0)
    for j in np.arange(1,Ny):
        for i in np.arange(0,Nx-x0):
            for k in np.arange(0,Nx-x0):
                temp[k] = np.sqrt(abs(r[i,j] - r_wall[k]))
            mindist[i,j] = min(temp)
    dist=np.zeros((Nx-x0,Ny-1))
    for i in np.arange(0,Nx-x0):
        for j in np.arange(0,int(yy.shape[1]-1)):
            dist[i,j] = yy[i,j+1]-yy[i,j]
    #print(dist.min())
    return xx,yy,vel,nut,mindist

#if __name__ == '__main__':
#    name,x0,Nx,Ny = sys.argv
#    xx,yy,vel,nut=main(np.int(x0),np.int(Nx),np.int(Ny))
    #plt.plot(x2,yup[:,1],'-or',yy[0,:],nut,'ok')
    #plt.plot(x1,yu[:,1],'-or',yy[0,:],vel,'ok')
    #plt.show()
xx,yy,vel,nut_dat,mindist = main(x0=24,Nx=409,Ny=109)
#print(xx[0,0])
#plt.plot(yy[0,:],nut_dat,'ok')
#plt.plot(yy[0,:],vel,'ok')
#plt.show()