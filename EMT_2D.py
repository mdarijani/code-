import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio




# DEFINING SOME VARIABLES
# HOMOGENEOUS (BACKGROUND) COMPLIANCE MATRIX 
rho = 2285.0  
Vp = 2954.0
Vs = 1829.0

# SIMULATION PARAMETERS (COME FROM SEM)
first_timestep = 225
timestep = 225
total_timestep = 337500

total_timestep_rcv = 337500

f = 500000. # frequency

dt = 0.4e-9 # dt!

dx = 0.15 / 1000 # come from after interpolation
dz = dx

# Center of the interested area
xcc = 0.075
zcc = 0.075

wavelength = Vp / f
# X and Z size of the interested area
xi = 10. * wavelength
zi = 10. * wavelength



C_b_ij = np.zeros((3,3))
C_b_ij[0,0] = C_b_ij[1,1] = rho * Vp**2  
C_b_ij[2,2] = rho * Vs**2
C_b_ij[0,1] = C_b_ij[1,0] = C_b_ij[0,0] - 2.0 * C_b_ij[2,2]

S_b_ij = np.linalg.inv(C_b_ij)



dummy = np.zeros(((total_timestep-first_timestep)/timestep)) 
timeStep = dummy.copy()

P_vel , S_vel = dummy.copy(), dummy.copy()

S_ij = np.zeros(((total_timestep-first_timestep)/timestep,3,3))
C_ij = np.zeros(((total_timestep-first_timestep)/timestep,3,3))

#ave_sig_xx, ave_sig_zz, ave_sig_xz = dummy.copy(), dummy.copy(), dummy.copy()
#ave_eps_zz, ave_eps_xz = dummy.copy(), dummy.copy()

bndAve_sig_xx, bndAve_sig_zz, bndAve_sig_xz = dummy.copy(), dummy.copy(), dummy.copy()
#bndAve_eps_xx, bndAve_eps_zz, bndAve_eps_xz = dummy.copy(), dummy.copy(), dummy.copy()

Z22, Z11 = dummy.copy(), dummy.copy()
uz, ux = dummy.copy(), dummy.copy()
uz2, ux2, uz1, ux1 = dummy.copy(), dummy.copy(), dummy.copy(), dummy.copy()

sigma_mat = np.zeros(((total_timestep-first_timestep)/timestep,2,2))
disp_mat = np.zeros(((total_timestep-first_timestep)/timestep,2,1))
Z_mat = np.zeros(((total_timestep-first_timestep)/timestep,2,1))



# DISPLACEMENT FILES (2media and homogeneous)
with open("Uz_file_single_2m.bin", 'rb') as f:
    disp_z = np.fromfile(f, dtype=np.float32)
    U_z = np.reshape(disp_z, [total_timestep_rcv,np.size(disp_z)/total_timestep_rcv], order='F')
with open("Ux_file_single_2m.bin", 'rb') as f:
    disp_x = np.fromfile(f, dtype=np.float32)
    U_x = np.reshape(disp_x, [total_timestep_rcv,np.size(disp_x)/total_timestep_rcv], order='F')
    
with open("Uz_file_single_h.bin", 'rb') as f:
    disp_z_h = np.fromfile(f, dtype=np.float32)
    U_z_h = np.reshape(disp_z_h, [total_timestep_rcv,np.size(disp_z_h)/total_timestep_rcv], order='F')
with open("Ux_file_single_h.bin", 'rb') as f:
    disp_x_h = np.fromfile(f, dtype=np.float32)
    U_x_h = np.reshape(disp_x_h, [total_timestep_rcv,np.size(disp_x_h)/total_timestep_rcv], order='F')


# CRACK BOUNDARY FILES (TOP AND BOTTOM COORDINATES)
with open("xcoor_crackbnd_top.txt", 'r') as f:
    xcT = np.genfromtxt(f)
    
with open("xcoor_crackbnd_bot.txt", 'r') as f:
    xcB = np.genfromtxt(f)



wvlt = 0#0.40 * wavelength

# Min and maz value of X and Z
xx_right = xcc + (xi/2)
zz_top = zcc - (zi/2)
xx_left = xcc - (xi/2)
zz_bottom = zcc + (zi/2) - wvlt

# dist_bot = (zi/2) - wvlt

# Index of the area of interest
ind_x_right = int(round(xx_right / dx))
ind_z_top = int(round(zz_top / dz))
ind_x_left = int(round(xx_left / dx))
ind_z_bottom = int(round(zz_bottom / dz))


# Making an odd matrix
if (ind_x_right - ind_x_left) % 2 == 0 :
    ind_x_right = ind_x_right + 1
    
if (ind_z_bottom - ind_z_top) % 2 == 0 :
    ind_z_bottom = ind_z_bottom + 1


# Surface of area
S = xi * (zi - wvlt)
zcent = (zi - wvlt) /2



# x vector
mm = np.arange(0,dx*(ind_x_right - ind_x_left),dx)
kkkk = mm - np.mean(mm)
kkkk = -kkkk
# z vector
mm_z = np.arange(0,dx*(ind_z_bottom - ind_z_top),dx)
kkkk_z = mm_z - np.mean(mm_z)
kkkk_z = -kkkk_z



for it in xrange(first_timestep, total_timestep, timestep) :
    
    k = (it-first_timestep)/timestep
    
    timeStep[k] = it
    
# Reading files    
#     with open("epsilon_xx_%4i.mat"%(it), 'r') as f:
#         ep1 = sio.loadmat(f)
#         eps_xx = ep1.get('epsilon_xx')
#     with open("epsilon_zz_%4i.mat"%(it), 'r') as f:
#         ep2 = sio.loadmat(f)
#         eps_zz = ep2.get('epsilon_zz')
#     with open("epsilon_xz_%4i.mat"%(it), 'r') as f:
#         ep12 = sio.loadmat(f)
#         eps_xz = ep12.get('epsilon_xz')
        
    with open("sigma_xx_%4i.mat"%(it), 'r') as f:
        si1 = sio.loadmat(f)
        sig_xx = si1.get('sigma_xx')
    with open("sigma_zz_%4i.mat"%(it), 'r') as f:
        si2 = sio.loadmat(f)
        sig_zz = si2.get('sigma_zz')
    with open("sigma_xz_%4i.mat"%(it), 'r') as f:
        si12 = sio.loadmat(f)
        sig_xz = si12.get('sigma_xz')
        
    

# Making a small matrix (i.e surface of interest)

    srf_sig_xx1 = np.copy(np.flipud(sig_xx))
    srf_sig_zz1 = np.copy(np.flipud(sig_zz))
    srf_sig_xz1 = np.copy(np.flipud(sig_xz)) 
#     srf_eps_zz1 = np.copy(np.flipud(eps_zz)) 
#     srf_eps_xz1 = np.copy(np.flipud(eps_xz))
    
    srf_sig_xx = np.copy(srf_sig_xx1[ind_z_top:ind_z_bottom, ind_x_left:ind_x_right])
    srf_sig_zz = np.copy(srf_sig_zz1[ind_z_top:ind_z_bottom, ind_x_left:ind_x_right])
    srf_sig_xz = np.copy(srf_sig_xz1[ind_z_top:ind_z_bottom, ind_x_left:ind_x_right]) 
#     srf_eps_zz = np.copy(srf_eps_zz1[ind_z_top:ind_z_bottom, ind_x_left:ind_x_right]) 
#     srf_eps_xz = np.copy(srf_eps_xz1[ind_z_top:ind_z_bottom, ind_x_left:ind_x_right])

    
# Average of the whole surface        
#     ave_sig_xx[k] = (1/S) * np.trapz(np.trapz(sig_xx,dx=dz,axis=0),dx=dx)
#     ave_sig_zz[k] = (1/S) * np.trapz(np.trapz(sig_zz,dx=dz,axis=0),dx=dx)
#     ave_sig_xz[k] = (1/S) * np.trapz(np.trapz(sig_xz,dx=dz,axis=0),dx=dx)
    
#     ave_eps_zz[k] = (1/S) * np.trapz(np.trapz(eps_zz,dx=dz,axis=0),dx=dx)
#     ave_eps_xz[k] = (1/S) * np.trapz(np.trapz(eps_xz,dx=dz,axis=0),dx=dx)
    
    
# Average of the surface of interest        
#     ave_sig_xx[k] = (1/S) * np.trapz(np.trapz(srf_sig_xx,dx=dz,axis=0),dx=dx)
#     ave_sig_zz[k] = (1/S) * np.trapz(np.trapz(srf_sig_zz,dx=dz,axis=0),dx=dx)
#     ave_sig_xz[k] = (1/S) * np.trapz(np.trapz(srf_sig_xz,dx=dz,axis=0),dx=dx)
    
#     ave_eps_zz[k] = (1/S) * np.trapz(np.trapz(srf_eps_zz,dx=dz,axis=0),dx=dx)
#     ave_eps_xz[k] = (1/S) * np.trapz(np.trapz(srf_eps_xz,dx=dz,axis=0),dx=dx)

# # Summing the whole surface  
#     su_sig_xx[k] = np.sum(sig_xx)
#     su_sig_zz[k] = np.sum(sig_zz)
#     su_sig_xz[k] = np.sum(sig_xz)
            
#     su_eps_zz[k] = np.sum(eps_zz)
#     su_eps_xz[k] = np.sum(eps_xz)


# # Boundary for surface of interest    
    bndAve_sig_zz[k] = (1/S) * (((zcent) * np.trapz(srf_sig_zz[0,:], dx=dx)) + 
                                (-(-zcent) * np.trapz(srf_sig_zz[-1,:], dx=dx)) +
                               (-np.trapz((srf_sig_xz[:,0]*kkkk_z), dx=dx)) + 
                               (np.trapz((srf_sig_xz[:,-1]*kkkk_z), dx=dx)))

    bndAve_sig_xz[k] = (1/(2*S)) * (((zcent) * np.trapz(srf_sig_xz[0,:], dx=dx)) + 
                                    (np.trapz((srf_sig_zz[0,:]*(-kkkk)), dx=dx)) + 
                                    (-np.trapz((srf_sig_xx[:,0]*(kkkk_z)), dx=dz)) +
                                   (-(-xi/2) * np.trapz(srf_sig_xz[:,0], dx=dz)) + 
                                   (-(-zcent) * np.trapz(srf_sig_xz[-1,:], dx=dx)) + 
                                    (-np.trapz((srf_sig_zz[-1,:]*(-kkkk)), dx=dx)) +
                                    (np.trapz((srf_sig_xx[:,-1]*(kkkk_z)), dx=dz)) + 
                                    ((xi/2) * np.trapz(srf_sig_xz[:,-1], dx=dz)))
    
    

    uz2[k] = np.trapz(U_z[it,0:len(U_z[0,:])/2],xcT) - np.trapz(U_z[it,len(U_z[0,:])/2:len(U_z[0,:])],xcB)
    ux2[k] = np.trapz(U_x[it,0:len(U_x[0,:])/2],xcT) - np.trapz(U_x[it,len(U_x[0,:])/2:len(U_x[0,:])],xcB)
    
    uz1[k] = np.trapz(U_z_h[it,0:len(U_z_h[0,:])/2],xcT) - np.trapz(U_z_h[it,len(U_z_h[0,:])/2:len(U_z_h[0,:])],xcB)
    ux1[k] = np.trapz(U_x_h[it,0:len(U_x_h[0,:])/2],xcT) - np.trapz(U_x_h[it,len(U_x_h[0,:])/2:len(U_x_h[0,:])],xcB) 
    
    uz[k] = uz2[k] - uz1[k] #
    ux[k] = ux2[k] - ux1[k] #

    disp_mat[k,0,0] = (1/S) * ux[k]
    disp_mat[k,1,0] = (1/S) * uz[k]
    
    sigma_mat[k,0,0] = bndAve_sig_xz[k]
    sigma_mat[k,1,1] = bndAve_sig_zz[k]
    
#     Z_mat[k,:,:] = np.matmul(np.linalg.inv(sigma_mat[k,:,:]),disp_mat[k,:,:])
    
    Z11[k] = ((1/S) * ux[k]) / bndAve_sig_xz[k]
    Z22[k] = ((1/S) * uz[k]) / bndAve_sig_zz[k]
#     Z11[k] = ((1/S) * ux[k]) / ave_sig_xz[k]
#     Z22[k] = ((1/S) * uz[k]) / ave_sig_zz[k]    

    S_ij[k,0,0] = S_b_ij[0,0]
    S_ij[k,0,1] = S_ij[k,1,0] = S_b_ij[0,1]
    S_ij[k,1,1] = S_b_ij[1,1] + Z22[k]
    S_ij[k,2,2] = S_b_ij[2,2] + Z11[k]
#     S_ij[k,1,1] = S_b_ij[1,1] + Z_mat[k,1,0]
#     S_ij[k,2,2] = S_b_ij[2,2] + Z_mat[k,0,0]
    
    
#     S_ij[k,1,1] = (bndAve_eps_zz[k] - S_b_ij[0,1] * bndAve_sig_xx[k]) / bndAve_sig_zz[k]
#     S_ij[k,2,2] = bndAve_eps_xz[k] / bndAve_sig_xz[k]
    
#     S_ij[k,1,1] = (ave_eps_zz[k] - S_b_ij[0,1] * ave_sig_xx[k]) / ave_sig_zz[k]
#     S_ij[k,2,2] = ave_eps_xz[k] / ave_sig_xz[k]   
    

    
#     S_ij[k,1,1] = (strain_x_z_xz_sum[1] - S_b_ij[0,1] * stress_x_z_xz_sum[0]) / stress_x_z_xz_sum[1]
#     S_ij[k,2,2] = strain_x_z_xz_sum[2] / stress_x_z_xz_sum[2]
#     if math.isnan(S_ij[k,1,1]):
#         S_ij[k,1,1] = S_b_ij[1,1]
#     if math.isnan(S_ij[k,2,2]):
#         S_ij[k,2,2] = S_b_ij[2,2]
        
    
    C_ij[k,:,:] = np.linalg.inv(S_ij[k,:,:])
    
    
    S_vel[k] = np.sqrt(abs(C_ij[k,2,2]) / rho)
    P_vel[k] = np.sqrt(abs(C_ij[k,1,1]) / rho)

    
    print P_vel[k] , timeStep[k] , S_ij[k,1,1] , C_ij[k,1,1] #, Z11, Z22

print "done"



sio.savemat('00000eff_vp.mat', {'eff_vp': P_vel, 'timestep': timeStep})



