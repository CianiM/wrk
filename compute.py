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
from cancellami import xx,yy, vel,nut_dat, mindist
from dnami import np, sys # re-import non dnami modules
import os

from owned_io import write_data

# ===================================================================  FOLDERS

# -- restarts 
try :
    os.mkdir('./restarts/')
except FileExistsError:
    pass

# ========================================================================= ASK

# Parameters for the case ...
alpha = dn.cst(2.5)
Ma    = dn.cst(0.1)
Re    = dn.cst(936000)
Pr    = dn.cst(0.71)
gamma = dn.cst(1.4)
R     = dn.cst(287.05)
Rho0  = dn.cst(1.0)
U0    = dn.cst(1.0)        #Reference value of Velocity 
V0    = dn.cst(0.0)
# T0    = dn.cst(1.0)        #Reference value of Temperature
nut0  = dn.cst(3.0)



dimPcr = False

P0   = dn.cst(1.0)/(Ma**2*gamma)
#P0   = dn.cst(0.0)
Rho0 = dn.cst(1.0)#P0/T0*Ma**2*gamma

U0   = dn.cst(1.0)
dimPcr = False
if dimPcr:
    Cv = dn.cst(1.0/(gamma-1.0))
else:
    Cv = dn.cst(1.0/(Ma**2*gamma*(gamma-1.0)))

T0   = P0/(Rho0*(gamma-1.))/dn.cst(Cv)
et0 = dn.cst(U0**2*0.5 + Cv*T0)

filtr_amp = dn.cst(0.1)    # filter amplitude

# ... in time ...
with_dt = dn.cst(3e-7)
#nitmax  = 50000 # for actual run
 # for test case 

# ... in space ...
#L = dn.cst(2.*np.pi) 
L = dn.cst( 1.0 )
H = dn.cst( 0.909050 )
with_length = [L,H]      # domain length in each direction
#with_grid   = [nx-2*hlo,64]   # number of points in each direction

# ... as fast as possible!
with_proc     = [1,7] # mpi proc. topology

restart,nrst = False,1800000

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
dtree['eqns']['coeff'][7][1] = dn.cst(2.0/3.0)                  #sigma
dtree['eqns']['coeff'][8][1] = dn.cst(0.41)                   #k
dtree['eqns']['coeff'][9][1] = dn.cst(0.1355/0.41**2 + (1.0 + 0.622)*3.0/2.0) #Cw1
dtree['eqns']['coeff'][10][1] = dn.cst(0.3)                   #Cw2
dtree['eqns']['coeff'][11][1] = dn.cst(2.0)                     #Cw3
dtree['eqns']['coeff'][12][1] = dn.cst(7.1)                   #Cv1
dtree['eqns']['coeff'][13][1] = dn.cst(1)                     #Ct1
dtree['eqns']['coeff'][14][1] = dn.cst(2)                     #Ct2
dtree['eqns']['coeff'][15][1] = dn.cst(1.2)                   #Ct3
dtree['eqns']['coeff'][16][1] = dn.cst(0.5)                   #Ct4
dtree['eqns']['coeff'][17][1] = dn.cst(3.0/2.0)               #sigmaI
dtree['eqns']['coeff'][18][1] = dn.cst(0.25)                  #esse
dtree['eqns']['coeff'][19][1] = dn.cst(0.42)                   #L_ref
dtree['eqns']['coeff'][20][1] = dn.cst(110.4/300.0)                 #suth constant
dtree['eqns']['coeff'][21][1] = dn.cst(L)                     #Reference Length
dtree['eqns']['coeff'][22][1] = dn.cst(gamma)                   #Gamma 
dtree['eqns']['coeff'][23][1] = Rho0                          #Gamma 
dtree['eqns']['coeff'][24][1] = T0                            #Gamma 
dtree['eqns']['coeff'][25][1] = P0                            #Gamma 
dtree['eqns']['coeff'][26][1] = nut0
dtree['eqns']['coeff'][27][1] = dn.cst(L)
dtree['eqns']['coeff'][28][1] = dn.cst(H)
dtree['eqns']['coeff'][29][1] = dn.cst(3.58)
dtree['eqns']['coeff'][30][1] = dn.cst(3.58)
dtree['eqns']['coeff'][31][1] = dn.cst(3.58)
dtree['eqns']['coeff'][32][1] = alpha
dtree['eqns']['coeff'][33][1] = V0
dtree['eqns']['coeff'][34][1] = dn.cst(0.25)
dtree['eqns']['coeff'][35][1] = dn.cst(0.0)
dtree['eqns']['coeff'][36][1] = dn.cst(0.2)
dtree['eqns']['coeff'][37][1] = dn.cst(0.0)
dtree['eqns']['coeff'][38][1] = dn.cst(0.2)
#dtree['eqns']['coeff'][34][1] = dn.cst(0.7)
#dtree['eqns']['coeff'][34][1] = dn.cst(0.9)
# .. shortcut key handles
numerics = dtree['num']
grid     = dtree['grid']['size']
geom     = dtree['grid']['geom']
mpi      = dtree['mpi']['split']

hlo  = numerics['hlo']
hhh = 409 - 25
kkk = 109
with_grid = [hhh-2*hlo,kkk-2*hlo]

geom['Lx'] = dn.cst(1.0)#with_length[0] 
geom['Ly'] = dn.cst(1.0)#ith_length[1] 
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
dksidx=dtree['eqns']['qvec']['views']['dksidx']
detady=dtree['eqns']['qvec']['views']['detady']
dksidy=dtree['eqns']['qvec']['views']['dksidy']
detadx=dtree['eqns']['qvec']['views']['detadx']
J=dtree['eqns']['qvec']['views']['J']
Jm1=dtree['eqns']['qvec']['views']['Jm1']
rho = dtree['eqns']['qvec']['views']['rho']
u = dtree['eqns']['qvec']['views']['u']
v = dtree['eqns']['qvec']['views']['v']
et = dtree['eqns']['qvec']['views']['et']
Pressure = dtree['eqns']['qvec']['views']['Pressure']
U_phy = dtree['eqns']['qvec']['views']['U_phy']
V_phy = dtree['eqns']['qvec']['views']['V_phy']
nut = dtree['eqns']['qvec']['views']['nut']
U_inlet = dtree['eqns']['qvec']['views']['U_inlet']
nut_inlet = dtree['eqns']['qvec']['views']['nut_inlet']
tau_wall = dtree['eqns']['qvec']['views']['tau_wall']
'''
visc_SA = dtree['eqns']['qvec']['views']['visc_SA']
visc_turb = dtree['eqns']['qvec']['views']['visc_turb']
Pressure = dtree['eqns']['qvec']['views']['Pressure']
chi_coeff = dtree['eqns']['qvec']['views']['chi_coeff']
P = dtree['eqns']['qvec']['views']['P']
stemp = dtree['eqns']['qvec']['views']['stemp']
D = dtree['eqns']['qvec']['views']['D']
Dis = dtree['eqns']['qvec']['views']['Dis']
'''
geom['dx'] = geom['Lx']/(grid['nxgb']+2*hlo-1.0)
geom['dy'] = geom['Ly']/(grid['nygb']+2*hlo-1.0)

uu = np.zeros((hhh,kkk))
nutnut = np.zeros((hhh,kkk))
for i in np.arange(0,hhh):
    uu[i,:] = vel[:]
    nutnut[i,:] = nut_dat[:]


ksi[:,:] = xx[dMpi.ibeg-1:dMpi.iend+2*hlo,dMpi.jbeg-1:dMpi.jend+2*hlo]
eta[:,:] = yy[dMpi.ibeg-1:dMpi.iend+2*hlo,dMpi.jbeg-1:dMpi.jend+2*hlo]
d[:,:] = mindist[dMpi.ibeg-1:dMpi.iend+2*hlo,dMpi.jbeg-1:dMpi.jend+2*hlo]

#print(ksi.shape)
#for j in np.arange(dMpi.jbeg-1,dMpi.jend+2*hlo):
#    U_inlet[j] = uu[0,j]
#    nut_inlet[j] = nutnut[0,j]
U_inlet[:,0] = uu[0,dMpi.jbeg-1:dMpi.jend+2*hlo]
#U_inlet[:,0] = dn.cst(1.0)
nut_inlet[:,0] = nutnut[0,dMpi.jbeg-1:dMpi.jend+2*hlo]
#nut_inlet[:,0] = dn.cst(3.0)


write_data(['ksi','eta'],'./grid/',0,0,dtree)


# - Store variables aliases if any
if 'qstored' in dtree['eqns']['qvec']['views'].keys():
	qstored  = dtree['eqns']['qvec']['views']['qstored'] 

# ================================================================== FUNCTIONS 

def sound_speed():
    e = et - .5*(u*u+v*v)
    T = (1./alpha)*( e*Ma*Ma )
    c = np.sqrt( T*(1.+1./alpha) )/Ma
    return c    

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


rho[:,:]  = Rho0
u[:,:]    = uu[dMpi.ibeg-1:dMpi.iend+2*hlo,dMpi.jbeg-1:dMpi.jend+2*hlo]
#u[:,:] = dn.cst(1.0)
#u[:,:]    = dn.cst(0.0)
v[:,:]    = dn.cst(0.0)
#nut[:,:]  = nutnut[:,:]
nut[:,:] = nutnut[dMpi.ibeg-1:dMpi.iend+2*hlo,dMpi.jbeg-1:dMpi.jend+2*hlo]
#nut[:,:] = dn.cst(3.0)
#et[:,:] = et0
for j in np.arange(0,ny+2*hlo):
    for i in np.arange(0,nx+2*hlo):
        #p = dn.cst(0.0) 
        p = P0
        et[i,j] = p/rho[i,j]*dn.cst(1./(gamma-1.)) + dn.cst(0.5)*( u[i,j]**2 + v[i,j]**2 )


if restart:
    from genrst import genRst

    genRst(nrst)

    dn.dnami_io.read_restart(dtree,fname='restart')    
    ti = dtree['eqns']['time']
    ni = dtree['num']['tint']['itn'] + 1 

# -- Swap 
dMpi.swap(q,hlo,dtree)
field = ['rho','u','v','et','nut','Pressure','tau_wall','U_phy','V_phy']
#field = ['rho','u','v','et','nut','tau_wall','visc_SA','visc_turb','Pressure','chi_coeff']

#field = ['rho','u','v','et','nut','tau_wall','visc_SA','visc_turb','Pressure','chi_coeff','P','D','Dis','stemp']
#field = ['rho','u','v','et','Pressure','U_phy','V_phy']
# -- Write the first restart
dn.dnami_io.write_restart(0,ti,0,dtree)

path  = './output/'

write_data(field,path,0,0,dtree)

# ========================================================================= RUN



intparam,fltparam,data = (dtree['libs']['fort']['integers'],
						  dtree['libs']['fort']['floats'],
						  dtree['libs']['fort']['data'])

dn.dnamiF.applybc(intparam,fltparam,data)
if 'qstored' in dtree['eqns']['qvec']['views'].keys():
    dMpi.swap(  qstored,hlo,dtree)
    dn.dnamiF.stored(intparam,fltparam,data,1)
    dn.dnamiF.stored(intparam,fltparam,data,0)	
    dMpi.swap(qstored,hlo,dtree)



nitmax  = 50000000
mod_filter = 1
mod_output = 100000
mod_info   = 1000



dn.dnamiF.applybc(intparam,fltparam,data)


for n in range(ni,nitmax+ni):

    ti = ti + dt

    # - RK loop
    for nrk in range(1,4):

        intparam[7] = nrk
        dMpi.swap(	q,hlo,dtree)

        if 'qstored' in dtree['eqns']['qvec']['views'].keys():
        	dn.dnamiF.stored(intparam,fltparam,data)
        	dMpi.swap(	qstored,hlo,dtree)

        dn.dnamiF.time_march(intparam,fltparam,data)	
        #for i in np.arange(0,rho.shape[0]):
        #    for j in np.arange(0,rho.shape[1]):
        #        print(i,j,u[i,j])

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
            dMpi.swap(  q,hlo,dtree)
            write_data(field,path,n,ti,dtree)
            dn.dnami_io.write_restart(n,ti,0,dtree)
    # - Output information
    if np.mod(n,mod_info) == 0:
           # dn.dnami_io.write_data(['omg'],n,ti,dtree,'./out/omg/','omg')
            if dMpi.ioproc:
                print('____________________________________________________________')
                print('iteration',n,' with time t =',ti)
                sys.stdout.flush()

            dn.dnami_io.globalMinMax(dtree,rho,'rho')
            dn.dnami_io.globalMinMax(dtree,u,'u')
            dn.dnami_io.globalMinMax(dtree,v,'v')
            dn.dnami_io.globalMinMax(dtree,et,'et')
            dn.dnami_io.globalMinMax(dtree,nut,'nut')

            sys.stdout.flush()

            # Compute CFL:

            if dMpi.ioproc:
                print('convective CFL numbers')
                sys.stdout.flush()

            #cfl = dt*np.abs(u/dksidx)/(dx*1.0)
            #dn.dnami_io.globalMax(dtree,cfl,'cfl-x')
            #cfl = dt*np.abs(v/detady)/(dy*1.0)
            #dn.dnami_io.globalMax(dtree,cfl,'cfl-y')

            if dMpi.ioproc:
                print('acoustic CFL numbers')
                sys.stdout.flush()
            c = sound_speed()
        
            cfl = dt*(np.abs(u)+c)/dksidx/(dx*1.0)
            dn.dnami_io.globalMax(dtree,cfl,'cfl-x')
            cfl = dt*(np.abs(v)+c)/detady/(dy*1.0)
            dn.dnami_io.globalMax(dtree,cfl,'cfl-y')
            #U = np.sqrt(u*u+v*v); locMach = U/c
            #dn.dnami_io.globalMinMax(dtree,locMach,'M')
            if dMpi.ioproc:
                print('diffusive CFL numbers')
            #cfl = dt/detady/detady/(dy*dy*Re)
            dn.dnami_io.globalMax(dtree,cfl,'cfl-y')
            # print('dt/(Re*dy*dy)',cfl)
            sys.stdout.flush()
            if dMpi.ioproc: sys.stdout.flush()  
 
