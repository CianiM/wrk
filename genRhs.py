import sys
import re
import numpy as np
from genKer import  genrk3, genrk3update, genFilter, genBC, append_Rhs
import os


wp = 'float64'
  
from  equations_MA_SA import *

# Set install path (needed by genKer)
instpath = os.environ['INSTALLPATH']
incPATH  = instpath+'/src_for/includes/gen/'

hlo_glob = 2

def main():
      
    from genKer import rhs_info    
    rhs = rhs_info(dim,wp,hlo_glob,incPATH,varsolved,varname, 
                   consvar=consvar,varstored=varstored,varloc=varloc,varbc=varbc,
                   coefficients=coefficients)


# Generate LHS:
    genrk3(len(varsolved)      ,rhs=rhs) 
    genrk3update(len(varsolved),rhs=rhs)

# Generate RHS:
    append_Rhs(Src_conv, 5,4, rhsname, locname_conv, update=False,rhs=rhs,stored=True)
    append_Rhs(Src_dif , 3,2 , rhsname, locname_dif , update=True ,rhs=rhs)  

# Generate Filters (if required):      
    genFilter(5,4, len(varsolved),rhs=rhs)     #No filter required

# Generate BCs:
    genBC(Src_conv ,5,4,rhsname , locname_conv, update=False,rhs=rhs,stored=True)
    genBC(Src_dif  ,3,2 ,rhsname , locname_dif , update=True ,rhs=rhs)

    #--j1
    genBC(Src_BC_phy_j1  ,3,2,rhsname , locname_bc,     setbc=[True,{'Low_surf':{'j1':['q']}}]   , update=False, rhs=rhs)
    #genBC(Src_BC_rhs_j1  ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'Low_surf':{'j1':['rhs']}}] , update=False, rhs=rhs )
    genBC(Src_BC_Wall  ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'Low_surf':{'j1':['rhs']}}] , update=False, rhs=rhs )
    #--i1
    genBC(src_rhsBC['i1']  ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'Inlet_surf':{'i1':['rhs']}}] , update=False, rhs=rhs )
    genBC(Src_BC_phy_i1  ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'Inlet_surf':{'i1':['q']}}] , update=False, rhs=rhs )
    
    #--imax
    genBC(src_rhsBC['imax'],3,2,rhsname , locname_bc_rhs, setbc=[True,{'Out_surf':{'imax':['rhs']}}] , update=False, rhs=rhs )


    #--jmax
    #genBC(src_rhsBC['jmax'],3,2,rhsname , locname_bc_rhs, setbc=[True,{'Top_surf':{'jmax':['rhs']}}] , update=False, rhs=rhs )
    #genBC({'nut':'nut0'}   ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'Top_surf':{'jmax':['q']}}] , update=False, rhs=rhs )
    genBC(Src_BC_phy_jmax  ,3,2,rhsname , locname_bc,     setbc=[True,{'Top_surf':{'jmax':['q']}}]   , update=False, rhs=rhs)
    genBC(Src_BC_Symm  ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'Top_surf':{'jmax':['rhs']}}] , update=False, rhs=rhs )

    #--i1j1
    #genBC(src_rhsBC['i1']  ,3,2,rhsname , locname_bc,     setbc=[True,{'corner_inwall':{'i1j1':['rhs']}}]   , update=False, rhs=rhs)
    genBC(src_rhsBC['i1']  ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'corner_inwall':{'i1j1':['rhs']}}]   , update=False, rhs=rhs )
    genBC(Src_BC_phy_i1    ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'corner_inwall':{'i1j1':['q']}}] , update=False, rhs=rhs )

    #--imaxj1
    genBC(Src_BC_phy_j1  ,3,2,rhsname , locname_bc,     setbc=[True,{'corner_outwall':{'imaxj1':['q']}}]   , update=False, rhs=rhs)
    genBC(Src_BC_rhs_j1  ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'corner_outwall':{'imaxj1':['rhs']}}] , update=False, rhs=rhs )

    #--i1jmax
    # genBC(src_rhsBC['i1jmax'],3,2,rhsname , locname_bc_rhs, setbc=[True,{'corner_inout':{'i1jmax':['rhs']}}] , update=False, rhs=rhs )
    # genBC(Src_BC_phy_i1      ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'corner_inout':{'i1jmax':['q']}}] , update=False, rhs=rhs )
    genBC(src_rhsBC['i1']  ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'corner_inout':{'i1jmax':['rhs']}}] , update=False, rhs=rhs )
    genBC(Src_BC_phy_i1   ,3,2,rhsname , locname_bc,     setbc=[True,{'corner_inout':{'i1jmax':['q']}}]   , update=False, rhs=rhs)
    # genBC(Src_BC_phy_i1    ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'corner_inout':{'i1jmax':['q']}}] , update=False, rhs=rhs )
   
    #--imaxjmax
    genBC(Src_BC_Symm,3,2,rhsname , locname_bc_rhs, setbc=[True,{'corner_outout':{'imaxjmax':['rhs']}}] , update=False, rhs=rhs )
    genBC(Src_BC_phy_jmax  ,3,2,rhsname , locname_bc,     setbc=[True,{'corner_outout':{'imaxjmax':['q']}}]   , update=False, rhs=rhs)
        
# Extract RHS info:
    rhs.export()

if __name__ == '__main__':
    main()
