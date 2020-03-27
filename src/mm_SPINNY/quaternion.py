import numpy as np

def quaternion(phi,theta,psi):

    ### calculations from "Euler Angles Quaternions Transformation Matrices"
    ### https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf
    ### Using ZXZ axis rotation sequence (pg.26)
    ### NOTE: There is an error in the formula given on page 26 for the rotation matrix, M[0,1]
    ### The corrected formula is given here 

    q0 = np.cos(0.5*theta)*np.cos(0.5*(phi+psi))
    q1 = np.sin(0.5*theta)*np.cos(0.5*(phi-psi))
    q2 = np.sin(0.5*theta)*np.sin(0.5*(phi-psi))
    q3 = np.cos(0.5*theta)*np.sin(0.5*(phi+psi))

    quat = np.array([q0,q1,q2,q3])

    rot = np.empty((3,3))

    rot[0,0]= -np.sin(phi)*np.cos(theta)*np.sin(psi)+np.cos(phi)*np.cos(psi)
    rot[0,1]= -np.sin(phi)*np.cos(theta)*np.cos(psi)-np.cos(phi)*np.sin(psi)
    rot[0,2]=  np.sin(phi)*np.sin(theta)
    rot[1,0]=  np.cos(phi)*np.cos(theta)*np.sin(psi)+np.sin(phi)*np.cos(psi)
    rot[1,1]=  np.cos(phi)*np.cos(theta)*np.cos(psi)-np.sin(phi)*np.sin(psi)
    rot[1,2]= -np.cos(phi)*np.sin(theta)
    rot[2,0]=  np.sin(theta)*np.sin(psi)
    rot[2,1]=  np.sin(theta)*np.cos(psi)
    rot[2,2]=  np.cos(theta)

    return(quat,rot)

def mat2quat(mat):
    
    q0 = np.sqrt(1+mat[0,0]+mat[1,1]+mat[2,2])/2
    q1 = (mat[2,1]-mat[1,2])/(4.0*q0)
    q2 = (mat[0,2]-mat[2,0])/(4.0*q0)
    q3 = (mat[1,0]-mat[0,1])/(4.0*q0)

    quat = np.array([q0,q1,q2,q3])
    return(quat)

def quat2mat(quat):
    rot = np.empty([3,3])
    qr = quat[0]
    qi = quat[1]
    qj = quat[2]
    qk = quat[3]

    rot[0,0]= 1.-2.*(qj*qj+qk*qk)
    rot[0,1]=   2.*(qi*qj-qk*qr)
    rot[0,2]=   2.*(qi*qk+qj*qr)
    rot[1,0]=   2.*(qi*qj+qk*qr)
    rot[1,1]= 1.-2.*(qi*qi+qk*qk)
    rot[1,2]=   2.*(qj*qk-qi*qr)
    rot[2,0]=   2.*(qi*qk-qj*qr)
    rot[2,1]=   2.*(qj*qk+qi*qr)
    rot[2,2]= 1.-2.*(qi*qi+qj*qj)
    
    return(rot)
        
def quat2euler(quat): ## using the ZXZ/313 rotation sequence for Euler Angles, again from the NASA publication
    
    rot = quat2mat(quat)
    
    phi = np.arctan2(rot[0,2],rot[1,2])
    theta = np.arctan2( np.sqrt(1-rot[2,2]**2),rot[2,2])
    psi = np.arctan2(rot[2,0],rot[2,1])
    
    return(phi, theta, psi)
    
    
    
    
    
    