# -----------------------------------------------------------------------------
#
# CASE: 2D URANS equations - Flat Plate 
#       -- with kernels: rhs.py genRhs.py
#
# -----------------------------------------------------------------------------

# ======================================================================== LOAD
from re import T
from tkinter.tix import Tree
import dnami as dn # dNami kernel

from dnami import np, sys # re-import non dnami modules
from owned_io import write_data
import os

# ===================================================================  FOLDERS

# -- restarts 
try :
    os.mkdir('./restarts/')
except FileExistsError:
    pass

# ========================================================================= ASK

# Parameters for the case ...
alpha = dn.cst(2.5)
Ma    = dn.cst(0.2)
Re    = dn.cst(5000000)
Pr    = dn.cst(0.71)
gamma = dn.cst(1.4)
R     = dn.cst(287.05)
Rho0  = dn.cst(1.0)
U0    = dn.cst(1.0)        #Reference value of Velocity 
V0    = dn.cst(0.0)
T0    = dn.cst(1.0)        #Reference value of Temperature
nut0  = dn.cst(0.000025)



dimPcr = False
if dimPcr:
	Cv = dn.cst(1.0/(gamma-1.0))
else:
	Cv = dn.cst(U0**2/(T0*Ma**2*gamma*(gamma-1.0)))     #Non dimensional specific heat 

P0 = dn.cst(U0**2/(gamma*Rho0*Ma**2))                   #Pressure at the inlet
#T0 = dn.cst(P0/(Rho0*(gamma-1)*Cv))
et0 = dn.cst(U0**2*( Cv*T0 + 0.5))

filtr_amp = dn.cst(0.5)    # filter amplitude

# ... in time ...
with_dt = dn.cst(1e-9)
#nitmax  = 50000 # for actual run
nitmax  = 1000000  # for test case 

# ... in space ...
#L = dn.cst(2.*np.pi) 
L = dn.cst( 2.0 + 1.0/3 )
H = dn.cst(1.0)
with_length = [L,H]      # domain length in each direction
#with_grid   = [nx-2*hlo,64]   # number of points in each direction

# ... as fast as possible!
with_proc     = [1,5] # mpi proc. topology

# ===================================================================== PREPARE

dtree = dn.create_tree()

# .. assign user-defined values
dtree['eqns']['coeff'][0][1] = dn.cst(1.0/Re)
dtree['eqns']['coeff'][1][1] = dn.cst(1.0/( (gamma-1.0)*Ma**2*Re*Pr ))
dtree['eqns']['coeff'][2][1] = dn.cst(gamma-1.)
dtree['eqns']['coeff'][3][1] = Cv
dtree['eqns']['coeff'][4][1] = U0                    #Uref
dtree['eqns']['coeff'][5][1] = dn.cst(0.1355)                 #Cb1
dtree['eqns']['coeff'][6][1] = dn.cst(0.622)                  #Cb2
dtree['eqns']['coeff'][7][1] = dn.cst(2.0/3)                  #sigma
dtree['eqns']['coeff'][8][1] = dn.cst(0.41)                   #k
dtree['eqns']['coeff'][9][1] = dn.cst(0.1355/0.41**2+(1+0.622)/2.0/3) #Cw1
dtree['eqns']['coeff'][10][1] = dn.cst(0.3)                   #Cw2
dtree['eqns']['coeff'][11][1] = dn.cst(2)                     #Cw3
dtree['eqns']['coeff'][12][1] = dn.cst(7.1)                   #Cv1
dtree['eqns']['coeff'][13][1] = dn.cst(1)                     #Ct1
dtree['eqns']['coeff'][14][1] = dn.cst(2)                     #Ct2
dtree['eqns']['coeff'][15][1] = dn.cst(1.1)                   #Ct3
dtree['eqns']['coeff'][16][1] = dn.cst(2)                     #Ct4
dtree['eqns']['coeff'][17][1] = dn.cst(1.0/2.0/3)             #sigmaI
dtree['eqns']['coeff'][18][1] = dn.cst(0.25)                  #esse
dtree['eqns']['coeff'][19][1] = dn.cst(1.0)                   #L_ref
dtree['eqns']['coeff'][20][1] = dn.cst(110.4)                 #suth constant
dtree['eqns']['coeff'][21][1] = dn.cst(L)                     #Reference Length
dtree['eqns']['coeff'][22][1] = dn.cst(1.4)                   #Gamma 
dtree['eqns']['coeff'][23][1] = Rho0                          #Gamma 
dtree['eqns']['coeff'][24][1] = T0                            #Gamma 
dtree['eqns']['coeff'][25][1] = P0                            #Gamma 
dtree['eqns']['coeff'][26][1] = nut0

# .. shortcut key handles
numerics = dtree['num']
grid     = dtree['grid']['size']
geom     = dtree['grid']['geom']
mpi      = dtree['mpi']['split']

hlo  = numerics['hlo']
with_grid = [69-2*hlo,49-2*hlo]

geom['Lx'] = 1.0 #with_length[0] 
geom['Ly'] = 1.0 #with_length[1] 
geom['Lz'] = 0.

mpi['nxpr'] = with_proc[0] 
mpi['nypr'] = with_proc[1] 
mpi['nzpr'] = 1

grid['nxgb'] = with_grid[0] 
grid['nygb'] = with_grid[1]
grid['nzgb'] = 1
# .. start the message passing interface
dtree = dn.start_mpi(dtree)
dMpi = dtree['mpi']['dMpi']

nxgb = dMpi.nxgb
nygb = dMpi.nygb
nx  = dMpi.nx
ny  = dMpi.ny

# .. create the *hlocomputational grid and write to file
dtree = dn.create_grid(dtree)
dn.dnami_io.hello_world(dtree)

# define useful aliases
xloc, yloc = geom['xloc'], geom['yloc']
Lx  , Ly   = geom['Lx']  , geom['Ly']
dx  , dy   = geom['dx']  , geom['dy']


numerics['tint']['tstep'] = with_dt 
dt = numerics['tint']['tstep']
numerics['filtr']['eps'] = filtr_amp 

# .. allocate tree
large = 10000 #no cache blocking in this example
dtree['libs']['cache blocking'] = [large,large,large]
dtree = dn.allocate(dtree)

# - Primitive variables
q  = dtree['eqns']['qvec']['views']['q'] 

# - Metrics
d  = dtree['eqns']['qvec']['views']['d']
ksi=dtree['eqns']['qvec']['views']['ksi']
eta=dtree['eqns']['qvec']['views']['eta']
#xx=dtree['eqns']['qvec']['views']['xx']

rho = dtree['eqns']['qvec']['views']['rho']
u = dtree['eqns']['qvec']['views']['u']
v = dtree['eqns']['qvec']['views']['v']
et = dtree['eqns']['qvec']['views']['et']
#nut = dtree['eqns']['qvec']['views']['nut']

#with open('grid_coord.dat','r') as f:
#    grid = f.readlines()[2:]
#    coord=[word for line in grid for word in line.split()] 
#gridNASA = np.array(coord,dtype=np.float128)
#gridNASA = gridNASA.reshape((2,273,193))


xx=np.zeros((nxgb+2*hlo,nygb+2*hlo))
yy=np.zeros((nxgb+2*hlo,nygb+2*hlo))
with open('x_coord.dat','r') as f:
    dat=f.readlines()
    dat=np.array(dat)
#    xx=dat
    for j in np.arange(0,nygb+2*hlo):
#        ksi[:,j]=dat
        xx[:,j]=dat
with open('y_coord.dat','r') as f:
    dat=f.readlines()
    nelem = len(dat)
    dat=np.array(dat).T
#    yy=dat
    for i in np.arange(0,nxgb+2*hlo):
#        eta[i,:]=dat
        yy[i,:]=dat

#ksi[:,:] = gridNASA[0,dMpi.ibeg-1:dMpi.iend+2*hlo,dMpi.jbeg-1:dMpi.jend+2*hlo]
#eta[:,:] = gridNASA[1,dMpi.ibeg-1:dMpi.iend+2*hlo,dMpi.jbeg-1:dMpi.jend+2*hlo]
ksi[:,:] = xx[dMpi.ibeg-1:dMpi.iend+2*hlo,dMpi.jbeg-1:dMpi.jend+2*hlo]
eta[:,:] = yy[dMpi.ibeg-1:dMpi.iend+2*hlo,dMpi.jbeg-1:dMpi.jend+2*hlo]



deltaxI = dtree['eqns']['qvec']['views']['deltaxI']
deltayI = dtree['eqns']['qvec']['views']['deltayI']
# - Store variables aliases if any
if 'qstored' in dtree['eqns']['qvec']['views'].keys():
	qstored  = dtree['eqns']['qvec']['views']['qstored'] 

# ================================================================== FUNCTIONS 

#def sound_speed():
#    e = et - .5*(ux*ux+uy*uy)
#    T = (1./alpha)*( e*Ma*Ma )
#    c = np.sqrt( T*(1.+1./alpha) )/Ma
#    return c	
def temperature(et,u,v,Cv):
    dim=u.shape
    temp=np.zeros((dim[0],dim[1]))
    for i in range(dim[0]):
        for j in range(dim[1]):
            temp[i,j]=1/Cv*(et[i,j]-0.5*(u[i,j]**2+v[i,j]**2))
    return temp

# ================================================================== INITIALISE

# initial clock
ti = dn.cst(0.0)
ni = 1
trstart = dn.cst(0.)

#init thermo
#T0   = dn.cst(1.0)
#P0   = dn.cst(1.0)/(Ma**2*gamma)
#Rho0 = dn.cst(1.0)#P0/T0*Ma**2*gamma

#numpy slice refering to the core of the domain
#dom = np.s_[0:273,0:193]

#rho[dom] = Rho0
#u[dom] = U0
#v[dom] = V0
#et[dom] = et0
#p[dom]  = P0
rho[:,:] = Rho0
u[:,:] = U0
v[:,:] = V0
et[:,:] = et0

# -- Swap 
#print('---swap---')
dMpi.swap(q,hlo,dtree)
#print('---swap---')


# -- Write the first restart
dn.dnami_io.write_restart(0,ti,0,dtree)
field = ['rho','u','v','et']
path= './output'
write_data(field,path,0,0,dtree)
# ========================================================================= RUN

intparam,fltparam,data = (dtree['libs']['fort']['integers'],
						  dtree['libs']['fort']['floats'],
						  dtree['libs']['fort']['data'])

if 'qstored' in dtree['eqns']['qvec']['views'].keys():
    dn.dnamiF.stored(intparam,fltparam,data,1)	
    dn.dnamiF.stored(intparam,fltparam,data,0)
    dMpi.swap(qstored,hlo,dtree)

mod_filter = 1
mod_output = 500000
mod_info   = 100000.


for n in range(1,nitmax+1):
    ti = ti + dt
#    print(eta[0,:])
#    te=temperature(et,u,v,Cv)
#    print(str(ti))
#    e = et - .5*(u*u)
#    p = eos_p(rho,e)
#    f=open(str(ti)+'Uvelocità.dat','w')
#    vv=open(str(ti)+'Vvelocità.dat','w')
#    pres=open(str(ti)+'Pressure.dat','w')
#    dens=open(str(ti)+'Density.dat','w')
#    temp=open(str(ti)+'Tempera.dat','w')
#    stem=open(str(ti)+'stemp.dat','w')
#    deltxI=open(str(ti)+'deltxI.dat','w')
#    deltyI=open(str(ti)+'deltyI.dat','w')
#    ksiout=open(str(ti)+'ksi.dat','w')
#    etaout=open(str(ti)+'eta.dat','w')
#    sacout=open(str(ti)+'dksidx.dat','w')
#    sacout2=open(str(ti)+'detady.dat','w')
#    sacout3=open(str(ti)+'dksidy.dat','w')
#    sacout4=open(str(ti)+'detadx.dat','w')
#    viscos=open(str(ti)+'viscosit.dat','w')
#    coeff1=open(str(ti)+'fw.dat','w')
#    coeff2=open(str(ti)+'gg.dat','w')
#    coeff3=open(str(ti)+'fv2.dat','w')
    udm=u.shape
#    print(udm)
#    print(ksi.shape)
#    print(eta.shape)
#    print(ksi[:,0])
#    print(ksi[0,:])
#    print(sacc[0,:])
#    deriv = np.zeros(len(ksi[1:,0]-1))
#    for i in np.arange(1,len(ksi[:,0])-1):
#        deriv[i]= (ksi[i+1,0]-2*ksi[i,0]+ksi[i-1,0])/dx
#    cic=(ksi[2,0]-2*ksi[1,0]+ksi[0,0])/dx
#    print('deriv=', deriv)
#    print(cic)
#    for i in np.arange(0,udm[0]):
#        for j in np.arange(0,udm[1]):
#            f.write('Uvel= '+ str(u[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            vv.write('Vvel= '+ str(v[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            pres.write('Pres= '+ str(p[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            dens.write('Rho= '+ str(rho[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            temp.write('T= '+ str(te[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            stem.write('stem= '+ str(stemp[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            deltxI.write('deltaxI= '+ str(deltaxI[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            deltyI.write('deltayI= '+ str(deltayI[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            ksiout.write('ksi= '+ str(ksi[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            etaout.write('eta= '+ str(eta[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            sacout.write('dksidx= '+ str(dksidx[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            sacout2.write('detady= '+ str(detady[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            sacout3.write('dksidy= '+ str(dksidy[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            sacout4.write('detadx= '+ str(detadx[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            viscos.write('visc= '+ str(viscosità[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            coeff1.write('fw= '+ str(fw[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            coeff2.write('gg= '+ str(gg[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#            coeff3.write('fv2= '+ str(fv2[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
#    f.close()
#    vv.close()
#    pres.close()
#    dens.close()
#    temp.close()
#    stem.close()
#    deltxI.close()
#    deltyI.close()
#    etaout.close()
#    ksiout.close()
#    sacout.close()
#    sacout2.close()
#    sacout3.close()
#    sacout4.close()
#    viscos.close()
#    coeff1.close()
#    coeff2.close()
#    coeff3.close()
#    print(d[:,0])
#    print(deltayI[0,:])
#    print(Cv)
#    print('----ksi----')
#    print(dtree['eqns']['qvec']['views']['ksi'][0,:])
#    print('----ksi----')
    # - RK loop
    for nrk in range(1,4):
        intparam[7] = nrk
        dMpi.swap(	q,hlo,dtree)

        if 'qstored' in dtree['eqns']['qvec']['views'].keys():
        	dn.dnamiF.stored(intparam,fltparam,data)
        	dMpi.swap(	qstored,hlo,dtree)

        dn.dnamiF.time_march(intparam,fltparam,data)	

    # - Filter
    if np.mod(n,mod_filter) == 0:

        dMpi.swapXc(q,hlo,dtree)
        dn.dnamiF.filter(1,intparam,fltparam,data)
        
        dn.dnamiF.applybc(intparam,fltparam,data)

        dMpi.swapYc(q,hlo,dtree)
        dn.dnamiF.filter(2,intparam,fltparam,data)
        
        dn.dnamiF.applybc(intparam,fltparam,data)


    # - Output restarts 
    if np.mod(n,mod_output) == 0:
            dn.dnami_io.write_restart(n,ti,0,dtree)
            write_data(field,path,n,ti,dtree)
    # - Output information
    if np.mod(n,mod_info) == 0:

        if dMpi.ioproc:
            print('____________________________________________________________')
            print('iteration',n,' with time t =',ti)
            sys.stdout.flush()
        dn.dnami_io.globalMinMax(dtree,rho,'r')
        dn.dnami_io.globalMinMax(dtree,u,'u')
        dn.dnami_io.globalMinMax(dtree,v,'v')
        dn.dnami_io.globalMinMax(dtree,et,'et')
        if dMpi.ioproc:
            print('convective CFL numbers')
            sys.stdout.flush()
        cfl = dt*np.abs(u)/dx
        dn.dnami_io.globalMax(dtree,cfl,'cfl-x')
        if dMpi.ioproc:
                print('acoustic CFL numbers')
                sys.stdout.flush()
#        cfl = dt*(np.abs(u[hlo:nx+hlo])+c)/dx
#        dn.dnami_io.globalMax(dtree,cfl,'cfl-x')
    
 
with open('u_sol.dat','w') as f:
     for i in np.arange(0,udm[0]):
        for j in np.arange(0,udm[1]):
           # f.write('Uvel= '+ str(u[i,j])+' ,i= '+ str(i)+' ,j= '+str(j)+'\n')
            f.write(str(u[i,j])+'\n')
with open('v_sol.dat','w') as f:
    for i in np.arange(0,udm[0]):
        for j in np.arange(0,udm[1]):
            f.write(str(v[i,j])+'\n')

with open('rho_sol.dat','w') as f:
    for i in np.arange(0,udm[0]):
        for j in np.arange(0,udm[1]):
            f.write(str(rho[i,j])+'\n')

with open('et_sol.dat','w') as f:
    for i in np.arange(0,udm[0]):
        for j in np.arange(0,udm[1]):
            f.write(str(et[i,j])+'\n')
print('done')
#   
# ----------------------------------------------------------------------------



