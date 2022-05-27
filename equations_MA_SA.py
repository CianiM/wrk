# =============================================================================
# 3D navier stokes equations  
# =============================================================================
from re import L
import sympy as sym
import numpy as np
import sys

# number of dimensions
dim = 2 

# coefficients ////////////////////////////////////////////////////////////////
coefficients = {'ReI'        : 1, 
                'kappa'      : 2,  
                'gamma_m1'   : 3,
                'Cv'         : 4,
                'U0'         : 5,
                'Cb1'        : 6,
                'Cb2'        : 7,
                'sigma'      : 8,
                'k'          : 9,
                'Cw1'        : 10,
                'Cw2'        : 11,
                'Cw3'        : 12,
                'Cv1'        : 13,
                'Ct1'        : 14,
                'Ct2'        : 15,
                'Ct3'        : 16,
                'Ct4'        : 17,
                'sigmaI'     : 18,
                'esse'       : 19,
                'L_ref'      : 20,
                'sut'        : 21,
                'L'          : 22,
                'gamma'      : 23,
                'Rho0'       : 24,
                'T0'         : 25,
                'P0'         : 26,
                'nut0'       : 27,
                'L_x'        : 28,
                'L_y'        : 29,
                'eta_1'      : 30,
                'eta_2'      : 31,
                'eta_3'      : 32,
                'Rgas'       : 33,
                'V0'         : 34,
                'sigmaBC'    : 35,
                'nut_filt_x' : 36,
                'nut_filt_y' : 37,
                'u_filt_x'   : 38,
                'u_filt_y'   : 39                 
                }

# unknowns to march in time ///////////////////////////////////////////////////
varname      = {'rho'  : 1,
                  'u'  : 2,  
                  'v'  : 3, 
                  'et' : 4,
                  'nut': 5}
                  
varsolved = ['rho','u','v','et','nut']

# of which
consvar      = [2,3,4,5] # are conservative variables         

# derived local variables /////////////////////////////////////////////////////

# -- Two divergence forms for inside and outside first derivative
#divops  = '  ( ( [u]_1x ) + ( [v]_1y ) ) '
#ddivops = '  ( ( {u}_1x ) + ( {v}_1y ) )'

dnutdx = ' J*(   ( [nut]_1x )*detady - ( [nut]_1y )*detadx ) '
dnutdy = ' J*( - ( [nut]_1y )*dksidy + ( [nut]_1y )*dksidx )'
ddnutdx= ' J*(   ( [nut]_1x )*detady - ( [nut]_1y )*detadx ) '
ddnutdy= ' J*( - ( [nut]_1y )*dksidy + ( [nut]_1y )*dksidx )'
drhodx = ' J*(   ( [rho]_1x )*detady - ( [rho]_1y )*detadx ) '
drhody = ' J*( - ( [rho]_1x )*dksidy - ( [rho]_1y )*dksidx ) '
# -- Srain tensor
#S =   ' J*( 0.5_wp*( dksidy*( {u}_1x ) - dksidx*( {u}_1y ) + detady*( {v}_1x ) - detadx*( {v}_1y ) ) '
S =   ' J*( 0.5_wp*( dksidy*( [u]_1x ) - dksidx*( [u]_1y ) + detady*( [v]_1x ) - detadx*( [v]_1y ) ) ) '
#for i in S.keys():
#    s = s + '('+ S[i] +')**2 + '
#s= '(( ( '+S['uv'] + ' )**2 + ' + '( '+S['vu'] +' )**2 )*2)**0.5'
#print(s)

varloc = {  'p'     : ' ( gamma_m1 )*rho*(e) ',
            'e'     : ' ( et-0.5_wp*( u*u+v*v ) ) ',
            'T'     : ' ( e )/( Cv ) ',
            'chi'   : ' (  nut*rho/visc ) ',
#                'visc': ' ( 1 + sut )/( T + sut )*T**1.5 ',                #Sutherland's law
            'visc'  : ' ( 1.0_wp )',
            'fv1'   : ' ( chi**3.0_wp /( chi**3.0_wp + Cv1**3.0_wp ) ) ',                                         #Adimensional turbulent viscosity
            'visc_t': ' ( max( ( visc + fv1*nut*rho )*ReI, ( visc )*ReI ) )',                    #Sum of adim dynamic viscosity and turbulent viscosity
#            'visc_t': ' ( visc*ReI )', 
            'fv2'   : ' ( 1.0_wp - chi/( 1.0_wp + chi*fv1) ) ',
            'ft2'   : ' ( Ct3*exp( -Ct4*chi**2.0_wp ) ) ',
                #'symm':'( ( sign( 1.0_wp, ksi) - 1.0_wp ) /(-2.0_wp) ) )' ,
            'c'     : ' ( gamma*p/rho )**0.5 ',
#            'Ssymm' : ' ( max( fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp), 0.0_wp)  + 0.00000000000001_wp ) ',
#            'rsymm' : ' ( min( ReI*( nut/(Ssymm*k**2.0_wp*d**2.0_wp ) ), 10.0_wp ) )',
 #           'gsymm' : ' ( rsymm + Cw2*( rsymm**6.0_wp - rsymm ) )',
 #           'fwsymm': ' ( gsymm*( ( 1.0_wp + Cw3**6.0_wp )/( gsymm**6.0_wp + Cw3**6.0_wp) )**( 1.0_wp/6.0_wp ) )',
#            'dist'  : ' ( ( eta**2.0_wp + ksi**2.0_wp)**0.5_wp )',
            'Utilde': ' (   detady*u - dksidy*v )',
            'Vtilde': ' ( - detadx*u + dksidx*v ) ',
            'x_xi'  : ' ( J*detady )',
            'x_eta' : ' ( - J*dksidy )',
            'y_xi'  : ' ( - J*detadx )',
            'y_eta' : ' ( J*dksidx )'}

varstored={ 'd'     :{'symb' : 'd' ,#( ( 1.0_wp - ( sign(1.0_wp, ksi) - 1.0_wp ) /(-2.0_wp) )*eta + ( ( sign(1.0_wp, ksi) - 1.0_wp ) /(-2.0_wp) )*dist )'  ,
                               'ind': 1 ,
                               'static': True}, # normal distance
            'eta'   :{'symb' : ' eta ' ,
                               'ind': 2 ,
                               'static': True},
            'ksi'   :{'symb' : ' ksi ' ,
                               'ind': 3,
                               'static': True},
            'detady':{'symb' : ' [ eta ]_1y ',
                               'ind': 4,
                               'static': True},
            'dksidy':{'symb' : ' [ ksi ]_1y ',
                               'ind': 5,
                               'static': True},
            'detadx':{'symb' : ' [ eta ]_1x ',
                               'ind': 6,
                               'static': True},
            'dksidx':{'symb' : ' ( [ ksi ]_1x ) ',
                               'ind': 7,
                               'static': True},
            'J'     :{'symb' : ' 1.0_wp / ( dksidx*detady - detadx*dksidy )',
                               'ind': 8,
                               'static': True},
            'Jm1'   :{'symb' : '( dksidx*detady - detadx*dksidy ) ',
                               'ind': 9,
                               'static': True},
            'Pressure':{'symb':'(gamma_m1)*rho*(e)',
                               'ind': 10,
                               'static': False},
            'Fnut'  :{'symb' : dnutdx,
                               'ind': 11,
                               'static': False},
            'Gnut'  :{'symb' : dnutdy,
                               'ind': 12,
                               'static': False},
            'Fnut2' :{'symb' : ddnutdx,
                               'ind': 13,
                               'static': False},
            'Gnut2' :{'symb' : ddnutdy,
                               'ind': 14,
                               'static': False},
            'stemp' :{'symb' : ' ( 2.0_wp*( dabs( '+S+') ) ) '  ,
                               'ind': 15,
                               'static': False},            
            'SS'    :{'symb' : '( max(stemp + fv2*ReI*nut/(k**2.0_wp*d**2.0_wp), 0.3_wp*stemp) + 0.000000000000001_wp ) ',
                               'ind': 16,
                               'static': False},
            'r'     :{'symb' : ' ( min( ReI*( nut/(SS*k**2.0_wp*eta**2.0_wp ) ), 10.0_wp ) ) ',
                               'ind': 17,
                               'static': False},
            'g'     :{'symb' : ' ( r + Cw2*( r**6.0_wp - r ) ) ',
                               'ind': 18,
                               'static': False},
            'fw'    :{'symb' : ' ( g*( ( 1.0_wp + Cw3**6.0_wp )/( g**6.0_wp + Cw3**6.0_wp) )**( 1.0_wp/6.0_wp ) ) ',
                               'ind': 19,
                               'static': False},
            'tau_wall':{'symb':' ( J*visc_t*( dksidx*( [u]_1y ) - dksidy*( [u]_1x ) ) )',
                               'ind': 20,
                               'static': False },
            'visc_SA':{'symb': 'nut*rho',
                               'ind': 21,
                               'static': False},
            'visc_turb':{'symb':'nut*rho*fv1',
                               'ind': 22,
                               'static': False},
            'drho_x':{'symb' : drhodx,
                               'ind': 23,
                               'static': False,},
            'drho_y':{'symb' : drhody,
                               'ind': 24,
                               'static': False},
            'U_phy' :{'symb' : ' J*(detady*u - dksidy*v) ',
                               'ind': 25,
                               'static': False},
            'V_phy' :{'symb' : ' J*( - detadx*u + dksidx*v ) ',
                               'ind': 26,
                               'static': False}}

# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = { 'rho': ' d(rho)/dt '  , 
            'u'  : ' d(rho u)/dt ',
            'v'  : ' d(rho v)/dt ', 
            'et' : ' d(rho et)/dt ',
            'nut': ' d(rho nut)/dt '}
       

# .. name tags to use for intermediate variables created by the constructor
locname_rhs = {'rho': ' rhs_rho ',
               'u'  : ' rhs_rhou ',
               'v'  : ' rhs_rhov ', 
               'et' : ' rhs_et  ',
               'nut': ' rhs_nut '}


locname_dif = {'rho': ' dif_rho ',
               'u'  : ' dif_rhou ',
               'v'  : ' dif_rhov ', 
               'et' : ' dif_et  ',
               'nut': ' dif_nut '}

locname_conv = {'rho': ' conv_rho ',
                'u'  : ' conv_rhou ',
                'v'  : ' conv_rhov ',
                'et' : ' conv_et  ',
                'nut': ' conv_nut '}

locname_bc  = {'rho': ' bc_rho ',
               'u'  : ' bc_u ',
               'v'  : ' bc_v ',
               'et' : ' bc_et  ',
               'nut': ' bc_nut '}

locname_bc_rhs= { 'rho' : ' rhs_rho ',
                  'u'   : ' rhs_u ' ,
                  'v'   : ' rhs_v ',
                  'et'  : ' rhs_et ',
                  'nut' : ' rhs_nut ' }
         
# RHS terms ///////////////////////////////////////////////////////////////////

# Euler 

Fx = {'rho' : ' ( rho*Utilde ) ',
      'u'   : ' ( rho*u*Utilde  + p*detady ) ',
      'v'   : ' ( rho*v*Utilde  - p*dksidy ) ',
      'et'  : ' ( ( rho*et + p )*Utilde ) ',
      'nut' : ' ( rho*nut*Utilde ) '}

Fy = {'rho' : ' ( rho*Vtilde ) ',
      'u'   : ' ( rho*u*Vtilde  - p*detadx ) ', 
      'v'   : ' ( rho*v*Vtilde  + p*dksidx ) ', 
      'et'  : ' ( ( rho*et + p )*Vtilde ) ',
      'nut' : ' ( rho*nut*Vtilde ) '} 

Src_conv={}

for key in Fx.keys():
    Src_conv[key]= 'J*( [ ' +Fx[key] + ' ]_1x ' + ' + ' + ' [ '+Fy[key]+ ' ]_1y ) '

######################################
#                                    #
# Navier-Stokes Diffusive terms only # 
#                                    #
######################################

#from NS_curvi2D import dNamiFlx_x, dNamiFlx_y

#Fx = { 'u' : ' -  Jm1*( '+dNamiFlx_x[1] + ')',}
'''
tau = { 'xx' :  ' visc_t*J*( 4.0_wp/3.0_wp*( [Utilde]_1x ) - 2.0_wp/3.0_wp*( [Vtilde]_1y ) + 2.0_wp/3.0_wp*dksidx*( - [v]_1x + [v]_1y )) ',
        'xy' :  ' visc_t*J*( [Utilde]_1y + [Vtilde]_1x ) ',
        'yx' :  ' visc_t*J*( [Vtilde]_1x + [Utilde]_1y ) ',
        'yy' :  ' visc_t*J*( 4.0_wp/3.0_wp*( 2.0_wp*( [v]_1y ) - [u]_1x ) '}
'''
E   = { 'xx' : ' J*( detady*( {u}_1x ) - detadx*( {u}_1y ) ) ',
        'xy' : ' J*( 0.5_wp*( {(detady*v-dksidy*u)}_1x + {(dksidx*u-detady*v)}_1y ) )' ,
#        'yx' : ' J*( 0.5_wp*( {(detady*v-dksidy*u)}_1x + {(dksidx*u-detady*v)}_1y ) )' ,
        'yy' : ' J*( - dksidy*( {v}_1x ) + dksidx*( {v}_1y ) ) '}

E_symm = { 'xx' : ' J*(   detady*( {u}_1x ) ) ',
           'xy' : ' J*( 0.5_wp*( {(detady*v-dksidy*u)}_1x ) ) ',
#           'yx' : ' J*( 0.5_wp*( {(detady*v-dksidy*u)}_1x ) )',
           'yy' : ' J*( - dksidy*( {v}_1x ) ) '}


divU   = ' J*( {Utilde}_1x + {Vtilde}_1y ) '
divU_x = ' J*( {Utilde}_1x ) '

tau_xy = { 'xx' : 'visc_t*( 2.0_wp*( '+E['xx']+' ) -1.0_wp/3.0_wp*'+divU+' ) ', 
           'xy' : 'visc_t*2.0_wp*( '+E['xy']+' ) ',
#           'yx' : 'visc_t*2.0_wp*( '+E['yx']+' ) ',
           'yy' : 'visc_t*( 2.0_wp*( '+E['yy']+' ) -1.0_wp/3.0_wp*'+divU+' ) '}

tau_xy_symm = { 'xx' : 'visc_t*( 2.0_wp*( '+E_symm['xx']+' ) -1.0_wp/3.0_wp*'+divU_x+' ) ', 
                'xy' : 'visc_t*2.0_wp*( '+E_symm['xy']+' ) ',
#                'yx' : 'visc_t*2.0_wp*( '+E_symm['yx']+' ) ',
                'yy' : 'visc_t*( 2.0_wp*( '+E_symm['yy']+' ) -1.0_wp/3.0_wp*'+divU_x+' ) '}

'''
tau = { 'xx' : ' visc_t*2.0_wp/3.0_wp*( 2.0_wp*( [u]_1x ) - [v]_1y ) ',
        'xy' : ' visc_t*( [u]_1y + [v]_1x ) ',
        'yx' : ' visc_t*( [v]_1x + [u]_1y ) ',
        'yy' : ' visc_t*2.0_wp/3.0_wp*( 2.0_wp*( [v]_1y ) - [u]_1x )' }
'''
q = { 'x': ' kappa*J*(   detady*( {T}_1x ) - detadx*( {T}_1y ) ) ',
      'y': ' kappa*J*( - dksidy*( {T}_1x ) + dksidx*( {T}_1y ) ) '}

q_symm = { 'x': '  kappa*J*detady*( {T}_1x ) ',
           'y': '- kappa*J*dksidy*( {T}_1x ) '}

f = { 'u'  : ' - ( ' + tau_xy['xx'] + ' ) ',
      'v'  : ' - ( ' + tau_xy['xy'] + ' ) ',
      'et' : ' - ( u*'+tau_xy['xx'] + '+ v*'+tau_xy['xy']+ ' + '+q['x'] + ' ) ',
      'nut': ' - ( sigmaI*ReI*(visc/rho + nut)*Fnut2 ) ' }

g = { 'u'  : ' - ( ' + tau_xy['xy'] + ' ) ',
      'v'  : ' - ( ' + tau_xy['yy'] + ' ) ',
      'et' : ' - ( u*'+tau_xy['xy'] + '+ v*'+tau_xy['yy']+ ' + '+q['y'] + ' ) ', 
      'nut': ' - ( sigmaI*ReI*(visc/rho + nut)*Gnut2 ) ' }

f_symm = { 'u'  : ' - ( ' + tau_xy_symm['xx'] + ' ) ',
           'et' : ' - ( u*'+tau_xy_symm['xx'] + '+ v*'+tau_xy_symm['xy']+ ' + '+q_symm['x'] + ' ) '}

g_symm = { 'u'  : ' - ( ' + tau_xy_symm['xy'] + ' ) ',
           'et' : ' - ( u*'+tau_xy_symm['xx'] + '+ v*'+tau_xy_symm['xy']+ ' + '+q_symm['y'] + ' ) '}

f_ad_wall = { 'et' : ' - ( u*'+tau_xy['xx'] + '+ v*'+tau_xy['xy']+ ' ) '}
g_ad_wall = { 'et' : ' - ( u*'+tau_xy['xy'] + '+ v*'+tau_xy['yy']+ ' ) '}

Fx = {}
Fy = {}
Fx_symm = {}
Fy_symm = {}
Fx_wall = {}
Fy_wall = {}

for key in f.keys():
      Fx[key] = ' (   detady*' + f[key] + ' - dksidy*' + g[key] + ' ) '
      Fy[key] = ' ( - detadx*' + f[key] + ' + dksidx*' + g[key] + ' ) '
for key in f_symm.keys():
    Fx_symm[key] = ' (   detady*' + f_symm[key] + ' - dksidy*' + g_symm[key] + ' ) '
    Fy_symm[key] = ' ( - detadx*' + f_symm[key] + ' + dksidx*' + g_symm[key] + ' ) '

Fx_wall = ' (   detady*' + f_ad_wall['et'] + ' - dksidy*' + g_ad_wall['et'] + ' ) '
Fy_wall = ' ( - detadx*' + f_ad_wall['et'] + ' + dksidx*' + g_ad_wall['et'] + ' ) '

'''

Fx = {'u'   : ' - ' + tau_x['u'],
      'v'   : ' - ' + tau_x['v'], 
      
      'et'  : ' - ' + q_x 
              ' - u*(visc_t*( 2.0_wp *( {u}_1x ) - 2.0_wp/3.0_wp * ( '+ ddivops +'  )))'
              ' - v*(visc_t*( ( {u}_1y ) + ( {v}_1x )))',
      'nut' : ' - ( ReI*( visc/rho + nut )*sigmaI*( { nut }_1x ) )'}

Fy = {'u'   : ' - ' + tau_y['u'],
      'v'   : ' - ' + tau_y['v'], 
      'et'  : ' - ' + q_y 
              ' - u*(visc_t*( ( {u}_1y ) + ( {v}_1x )))'
              ' - v*(visc_t*( 2.0_wp *( {v}_1y ) - 2.0_wp/3.0_wp * ( '+ ddivops +'  )))',
      'nut' : ' - ( ReI*( visc/rho + nut )*sigmaI*( { nut }_1y ) )'}     
# -- Divergence formulation



voc = ['u','v','et']
Fx = {}
Fy = {}
print(voc[0])
for i in np.arange(0,3):
      Fx[voc[i]] = dNamiFlx_x[i+1]
      Fy[voc[i]] = dNamiFlx_y[i+1]
#      print(voc[i])

'''
#Fx['nut'] = ' - ( ReI*( visc/rho + nut )*sigmaI*( '+ddnutdx+'  ) ) '
#Fy['nut'] = ' - ( ReI*( visc/rho + nut )*sigmaI*( '+ddnutdy+' ) ) '


Src_dif  = {}

#for i in np.arange(1,4):
#      Src_dif[voc[i]] = 'J*( [ Jm1*( ' + Fx[i] +') ]_1x )' + ' + ' + ' [ Jm1*('+ Fy[i]  +') ]_1y ) '

for key in Fx.keys():
    Src_dif[key]= ' J*( [ ' + Fx[key] +' ]_1x ' + ' + ' + ' [ '+ Fy[key]  + ' ]_1y ) '

#for key in Fx_symm.keys():
#      Src_dif_symm[key]= ' J*( [ ' + Fx_symm[key] +' ]_1x ) '

Src_dif['nut'] =   Src_dif['nut'] + ' \
                   - ReI*Cb2*sigmaI*( Fnut*Fnut + Gnut*Gnut ) \
                   - rho*Cb1*(1.0_wp-ft2)*SS*nut + rho*ReI*(Cw1*fw-Cb1/k**2.0_wp*ft2)*(nut/eta)**2.0_wp \
                   + sigmaI*ReI*( visc/rho + nut )*( ( drho_x )*Fnut + ( drho_y )*Gnut )' 




########################################################################################################## 
#----------------------------------------BOUNDARY CONDITIONS---------------------------------------------#
##########################################################################################################
print('start building boundary conditions')
#--Physical BC--
Src_BC_phy = {}
Src_BC_phy_i1={}
Src_BC_phy_j1={}
Src_BC_phy_imax={}
Src_BC_phy_jmax={}

#--RHS--
Src_BC_rhs ={}
Src_BC_rhs_i1={}
Src_BC_rhs_j1={}
Src_BC_rhs_imax={}
Src_BC_rhs_jmax={}

######################################
#                                    #
#----------Symmetry, j1--------------#  
#                                    #
######################################


Src_BC_Symm_conv ={}

Fx_symm_conv = {'rho' : ' ( rho*Utilde ) ',
                'u'   : ' ( rho*u*Utilde + p*detady ) ',
                #'v'  : ' rho*v*Utilde  - p*dksidy ',
                'et'  : ' ( ( rho*et+ p)*Utilde )',
                'nut' : ' ( rho*Utilde*nut ) '}


for key in Fx_symm.keys():
      Src_BC_Symm_conv[key] = 'J*( [ ' + Fx_symm_conv[key] + ' ]_1x ) '

#Src_BC_Symm_conv['nut'] = ' u*J*( ( [ nut ]_1x )*detady - ( [ v ]_1y )*detadx ) '

Src_BC_Symm_dif={}

for key in Fx_symm.keys():
      Src_BC_Symm_dif[key]= ' J*( [ ' + Fx_symm[key] +' ]_1x ) '

#Src_BC_Symm_dif['nut']= Src_BC_Symm_dif['nut'] + ' - ReI*Cb2*sigmaI*( ( [ nut ]_1x )**2.0_wp ) \
#                         - rho*( Cb1*(1.0_wp-ft2)*Ssymm*nut + ReI*(Cw1*fwsymm-Cb1/k**2.0_wp*ft2)*(nut/d)**2.0_wp ) '
#--Building Symmetry boundary conditions

Src_BC_Symm = {}
for key in Src_BC_Symm_conv.keys():
      if key in Src_BC_Symm_dif.keys():
            Src_BC_Symm[key] = '( '+ Src_BC_Symm_dif[key] + ' + ' + Src_BC_Symm_conv[key] + ' )' 
      else:
            Src_BC_Symm[key] = '( '+ Src_BC_Symm_conv[key]+ ' )'

######################################
#                                    #
#-------------Wall, j1---------------#  
#                                    #
######################################

Src_BC_Wall={}
Src_BC_Wall_conv = {}
Src_BC_Wall_conv['rho'] = Src_conv['rho']
Src_BC_Wall_conv['et']  = Src_conv['et']

#Src_BC_Wall_conv = { 'rho' : ' rho*( [v]_1y )',
#                     'et'  : ' ( rho*et +p )*( [v]_1y )'}
'''
Src_BC_Wall_dif = { 'et'  : ' - visc_t*( [u]_1y )*( [u]_1y ) - 4.0_wp/3.0_wp*visc_t*( [v]_1y )*( [v]_1y ) ' } #\
                             # - kappa*( [( {T}_1y )*deltayI ]_1y ) \   
                             # - kappa*( [( {T}_1x )*deltaxI ]_1x ) '} # Adiabatic wall
'''
Src_BC_Wall_dif = { 'et'  : ' J*( [ ' + Fx_wall +' ]_1x ' + ' + ' + ' [ '+ Fy_wall + ' ]_1y ) ' }


#--Building Wall boundary conditions
for key in Src_BC_Wall_conv.keys():
      if key in Src_BC_Wall_dif.keys():
            Src_BC_Wall[key]= Src_BC_Wall_conv[key] + ' + ' + Src_BC_Wall_dif[key] 
      else:
            Src_BC_Wall[key]= Src_BC_Wall_conv[key] 

'''
for key in Src_BC_Symm.keys():
      if key in Src_BC_Wall.keys():
            Src_BC_rhs_j1[key] = Src_BC_Symm[key] + ' + ' + Src_BC_Wall[key]
      else:
            Src_BC_rhs_j1[key] = Src_BC_Symm[key]
'''

##--Physical BC
Src_BC_phy_j1 = { 'u'  : ' 0.0_wp ',
                  'v'  : ' 0.0_wp ',
                  'nut': ' 0.0_wp '}


#################################################
#                                               #
#---- Characteristic Boundary Conditions -------#  
#                                               #
#################################################


##-- Charactertistics --##
from CharsForConsLaw import characteristics
from genNSBC import sympy2dNami

Char={}
Char=characteristics('Euler')


x      = sym.Symbol('x')
y      = sym.Symbol('y')
z      = sym.Symbol('z')
t      = sym.Symbol('t')

rho   = sym.Function('rho')(x,y,z,t)
u     = sym.Function('u')  (x,y,z,t)
v     = sym.Function('v')  (x,y,z,t)
w     = sym.Function('w')  (x,y,z,t)
et    = sym.Function('et') (x,y,z,t)
p     = sym.Function('p')  (x,y,z,t)

rho   = sym.Symbol('rho')
u     = sym.Symbol('u')  
v     = sym.Symbol('v')  
w     = sym.Symbol('w')  
et    = sym.Symbol('et') 
p     = sym.Symbol('p')  
gamma = sym.Symbol('gamma')

Q_CS = sym.Matrix([[rho],
                   [rho*u],
                   [rho*v],
                   [0.5*rho*(u*u+v*v)+1.0/(gamma-1)*p]])
                   # [rho*w],
                   # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])

Q    = sym.Matrix([[rho],
                   [u],
                   [v],
                   # [w],
                   [p]])


Pmatrix = Q_CS.jacobian(Q)


#############################
#
# i1 subsonic inflow :
#
#############################


Li_BC_i1sym = Char['xi'][0].copy()



Li_BC_i1 = []
for i in Li_BC_i1sym:
    Li_BC_i1.append(sympy2dNami(i))

for i in np.arange(0,len(Li_BC_i1)):
    Li_BC_i1[i]=' ( '+Li_BC_i1[i]+' )'

Li_BC_i1[3] = Li_BC_i1[2]
Li_BC_i1[1] = ' 0.5_wp*( gamma_m1 )*( ' + Li_BC_i1[3] + ' + ' +   Li_BC_i1[3] + ' ) '
Li_BC_i1[0] = ' 0.0_wp'  

Src_BC_phy_i1 = { 'u' : 'U_inlet',
                  'v' : 'V0',
                  #'et':  'P0/(rho*gamma_m1) + 0.5_wp*(U0*U0 + V0*V0)',
                  'et':  'Cv*T0 + 0.5_wp*( U_inlet*U_inlet + V0*V0 )',
                  'nut' : 'nut_inlet'}


#############################
#
# imax non-reflecive outflow:
#
#############################


Li_BC_imaxsym = Char['xi'][0].copy()

Li_BC_imax = []
for i in Li_BC_imaxsym:
    Li_BC_imax.append(sympy2dNami(i))

for i in np.arange(0,len(Li_BC_imax)):
    Li_BC_imax[i]=' ( '+Li_BC_imax[i]+' )'

Li_BC_imax[2] = 'sigmaBC*c*(1-M_imax*M_imax)/L_y*( p - P0)'


#############################
#
# jmax non-reflective outfow:
#
#############################

Mi_BC_jmaxsym = Char['eta'][0].copy()

Mi_BC_jmax = []
for i in Mi_BC_jmaxsym:
    Mi_BC_jmax.append(sympy2dNami(i))
for i in np.arange(0,len(Mi_BC_jmax)):
    Mi_BC_jmax[i]=' ( '+Mi_BC_jmax[i]+' )'

Mi_BC_jmax[2] = 'sigmaBC*c*(1-M_jmax*M_jmax)/L_x*( p - P0)'



Src_BC_phy_jmax = { 'v' : ' u*detadx/dksidx ',
                    'nut':' 0.0_wp '}

########################################
#
#
#   Generate RHS for conservative form
#
#
########################################

curvi = {'x':'xi','y':'eta'}

Pleft_inv = {}
for d in ['x','y']:
  Pleft_inv[d] = Char[curvi[d]][2].inv()

def setLiNames(bc,dir='x'):  
    
    pre = 'L'

    if dir == 'y': pre = 'M'
    if dir == 'z': pre = 'N'

    L0 = sym.Symbol(pre+'0_'+bc) # u-c
    L1 = sym.Symbol(pre+'1_'+bc) # u+c
    L2 = sym.Symbol(pre+'2_'+bc) # v
    L3 = sym.Symbol(pre+'3_'+bc) # u

    return [L0,L1,L2,L3]

src_rhsBC = {}
tdir = {'x':'y','y':'x'}

lodiRHSdNami = {}
for bc in ['1','max']:
  for dir in ['x','y']:

    idir = 'i'
    if dir == 'y':idir = 'j'
    if dir == 'z':idir = 'j'
    

    [L0,L1,L2,L3] = setLiNames(idir+bc,dir)

    lodiRHS       = Pmatrix*(Pleft_inv[dir].applyfunc(lambda x: sym.simplify(sym.factor(x)))*sym.Matrix([L0,L1,L2,L3]))
    # lodiRHS       = (Pleft_inv[dir].applyfunc(lambda x: sym.simplify(sym.factor(x)))*sym.Matrix([L0,L1,L2,L3]))
     
    # lodiRHS       = lodiRHS.applyfunc(lambda x:J*   (sym.expand(sym.simplify(x.subs(betaswp[dir][0],betaswp[dir][1]))))) #+ Ftild[tdir[dir]]

    lodiRHSdNami[idir+bc] = []
    for symb in lodiRHS:
        eqdNami = sympy2dNami(sym.simplify(symb)).replace('c','c_'+idir+bc)
        # eqdNami = eqdNami.replace('beta_'+dir,'beta_'+idir+bc)        
        lodiRHSdNami[idir+bc].append(eqdNami)

    src_rhsBC[idir+bc] = {'rho':lodiRHSdNami[idir+bc][0] ,#+ '+ ('+Ftild[tdir[dir]]['rho']+' )',
                          'u'  :lodiRHSdNami[idir+bc][1] ,#+ '+ ('+Ftild[tdir[dir]]['u'  ]+' )',
                          'v'  :lodiRHSdNami[idir+bc][2] ,#+ '+ ('+Ftild[tdir[dir]]['v'  ]+' )',
                          # 'w'  :,.replace('c','c_'+idir+bc)
                          'et' :lodiRHSdNami[idir+bc][3] }#+ '+ ('+Ftild[tdir[dir]]['et']+' )'}

'''
src_rhsBC['imax']['nut'] = ' v*( [ nut ]_1y )' \
                           ' - ( [ ( ReI*( visc/rho + nut )*sigmaI*( { nut }_1y ) ) ]_1y )'\
                           ' - ReI*Cb2*sigmaI*( [ nut ]_1y )**2.0_wp '\
                           ' - rho*( Cb1*(1.0_wp-ft2)*SS*nut + rho*ReI*(Cw1*fw-Cb1/k**2.0_wp*ft2)*(nut/d)**2.0_wp ) '
'''
src_rhsBC['imax']['nut'] = ' J*( ' + ' [ '+ Fy['nut']  + ' ]_1y ) \
                   - ReI*Cb2*sigmaI*( Gnut*Gnut ) \
                   - rho*Cb1*(1.0_wp-ft2)*SS*nut + rho*ReI*(Cw1*fw-Cb1/k**2.0_wp*ft2)*(nut/eta)**2.0_wp \
                   + sigmaI*ReI*( visc/rho + nut )*( ( drho_y )*Gnut )'  

print('bc for imax nut done')
#src_rhsBC['jmax']['nut'] = ' u*( [ nut ]_1x ) '\
#                           ' - ( [ ( ReI*( visc/rho + nut )*sigmaI*( { nut }_1x ) ) ]_1x )'\
#                           ' - ReI*Cb2*sigmaI*( [ nut ]_1x )**2.0_wp  '\
#                           ' - rho*( Cb1*(1.0_wp-ft2)*SS*nut + ReI*(Cw1*fw-Cb1/k**2.0_wp*ft2)*(nut/d)**2.0_wp ) '

print('bc for jmax nut done')
# src_rhsBC['i1j1'] = {'rho':lodiRHSdNami['i1'][0] + '+ ('+lodiRHSdNami['j1'][0]+' )',
#                      'u'  :lodiRHSdNami['i1'][1] + '+ ('+lodiRHSdNami['j1'][1]+' )',
#                      'v'  :lodiRHSdNami['i1'][2] + '+ ('+lodiRHSdNami['j1'][2]+' )',
#                      # 'w'  :,.replace('c','c_'+idir+bc)
#                      'et' :lodiRHSdNami['i1'][3] + '+ ('+lodiRHSdNami['j1'][3]+' )'}

src_rhsBC['i1jmax'] = {'rho':lodiRHSdNami['i1'][0] + '+ ('+lodiRHSdNami['jmax'][0]+' )'},
                     #'u'  :lodiRHSdNami['i1'][1] + '+ ('+lodiRHSdNami['jmax'][1]+' )',
                     #'v'  :lodiRHSdNami['i1'][2] + '+ ('+lodiRHSdNami['jmax'][2]+' )',
                     # 'w'  :,.replace('c','c_'+idir+bc)
                     # 'et' :lodiRHSdNami['i1'][3] + '+ ('+lodiRHSdNami['jmax'][3]+' )'}     


# src_rhsBC['imaxj1'] = {'rho':lodiRHSdNami['imax'][0] + '+ ('+lodiRHSdNami['j1'][0]+' )',
#                      'u'  :lodiRHSdNami['imax'][1] + '+ ('+lodiRHSdNami['j1'][1]+' )',
#                      'v'  :lodiRHSdNami['imax'][2] + '+ ('+lodiRHSdNami['j1'][2]+' )',
#                      # 'w'  :,.replace('c','c_'+idir+bc)
#                      'et' :lodiRHSdNami['imax'][3] + '+ ('+lodiRHSdNami['j1'][3]+' )'}

src_rhsBC['imaxjmax'] = {'rho':lodiRHSdNami['imax'][0] + '+ ('+lodiRHSdNami['jmax'][0]+' )',
                     'u'  :lodiRHSdNami['imax'][1] + '+ ('+lodiRHSdNami['jmax'][1]+' )',
                     'v'  :lodiRHSdNami['imax'][2] + '+ ('+lodiRHSdNami['jmax'][2]+' )',
                     # 'w'  :,.replace('c','c_'+idir+bc)
                     'et' :lodiRHSdNami['imax'][3] + '+ ('+lodiRHSdNami['jmax'][3]+' )'}






varbc = {'M_i1'   : {'symb'  : ' (u)/sqrt(gamma*Rgas*T) ', 
                    'ind'   : 1 ,
                    'static': False,
                    'face'  : 'i1'}, 
                      
         'M_imax' : {'symb'  : ' (u)/sqrt(gamma*Rgas*T) ', 
                    'ind'   : 2 ,
                    'static': False,
                    'face'  : 'imax'},
         'M_j1'   : {'symb'  : ' (v)/sqrt(gamma*Rgas*T) ', 
                    'ind'   : 1 ,
                    'static': False,
                    'face'  : 'j1'},                    
         'M_jmax' : {'symb'  : ' (v)/sqrt(gamma*Rgas*T) ', 
                    'ind'   : 2 ,
                    'static': False,
                    'face'  : 'jmax'},
         'c_i1'   : {'symb'   : 'sqrt(gamma*p/rho) ', 
                    'ind'   : 3 ,
                    'static': False,
                    'face'  : 'i1'},             
         'c_imax'   : {'symb'   : 'sqrt(gamma*p/rho) ', 
                    'ind'   : 4 ,
                    'static': False,
                    'face'  : 'imax'},
         'c_j1'   : {'symb'   : 'sqrt(gamma*p/rho) ', 
                    'ind'   : 3 ,
                    'static': False,
                    'face'  : 'j1'},                   
         'c_jmax'   : {'symb'   : 'sqrt(gamma*p/rho) ', 
                    'ind'   : 4 ,
                    'static': False,
                    'face'  : 'jmax'},


          # Li's definition i1: 
         'L0_i1' : {'symb'   : Li_BC_i1[0].replace('c','c_i1')+' ', 
                     'ind'   : 5,
                     'static': False,
                     'face'  : 'i1'}, 
            
         'L1_i1' : {'symb'   : Li_BC_i1[1].replace('c','c_i1')+' ', 
                     'ind'   : 6 ,
                     'static': False,
                     'face'  : 'i1'},
         'L2_i1' : {'symb'   : Li_BC_i1[2].replace('c','c_i1')+' ', 
                     'ind'   : 7 ,
                     'static': False,
                     'face'  : 'i1'},            
         'L3_i1' : {'symb'   : Li_BC_i1[3].replace('c','c_i1')+' ', 
                     'ind'   : 8 ,
                     'static': False,
                     'face'  : 'i1'},
           
          #Li's definition imax: 
         'L0_imax' : {'symb'   : Li_BC_imax[0].replace('c','c_imax')+' ', 
                     'ind'   : 9 ,
                     'static': False,
                     'face'  : 'imax'},            
         'L1_imax' : {'symb'   : Li_BC_imax[1].replace('c','c_imax')+' ', 
                     'ind'   : 10 ,
                     'static': False,
                     'face'  : 'imax'},
         'L2_imax' : {'symb'   : Li_BC_imax[2].replace('c','c_imax')+' ', 
                     'ind'   : 11  ,
                     'static': False,
                     'face'  : 'imax'},            
         'L3_imax' : {'symb'   : Li_BC_imax[3].replace('c','c_imax')+' ', 
                     'ind'   : 12 ,
                     'static': False,
                     'face'  : 'imax'},  
      'U_inlet' : {'symb': 'U_inlet',
                      'ind' : 13,
                      'static':True,
                      'face' : 'i1'},  
      'nut_inlet':{'symb' : 'nut_inlet',
                      'ind'  : 14,
                      'static':True,
                      'face' : 'i1'},  
          # Mi's definition jmax: 
         'M0_jmax' : {'symb'   : Mi_BC_jmax[0].replace('c','c_jmax')+' ', 
                     'ind'   : 5 ,
                     'static': False,
                     'face'  : 'jmax'},            
         'M1_jmax' : {'symb'   : Mi_BC_jmax[1].replace('c','c_jmax')+' ', 
                     'ind'   : 6 ,
                     'static': False,
                     'face'  : 'jmax'},
         'M2_jmax' : {'symb'   : Mi_BC_jmax[2].replace('c','c_jmax')+' ', 
                     'ind'   : 7 ,
                     'static': False,
                     'face'  : 'jmax'},            
         'M3_jmax' : {'symb'   : Mi_BC_jmax[3].replace('c','c_jmax')+' ', 
                     'ind'   : 8 ,
                     'static': False,
                     'face'  : 'jmax'}


          }           
