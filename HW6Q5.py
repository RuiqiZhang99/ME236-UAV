import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, tan, sqrt, pi
from scipy.spatial.transform import Rotation as R

def rotation_vector_to_matrix(rotation_vector):
    rotation = R.from_rotvec(rotation_vector)
    return rotation.as_matrix()

w_BE = np.array([1, 10, 0])
omega_BE = np.array([[0, 0, 10], [0, 0, -1], [-10, 1, 0]])
T_BE0 = np.eye(3)
RPY_BE0 = np.zeros_like(w_BE)
delta_t = 0.01
omega_no_ul = np.array([
    [0, 1, 10, 0],
    [-1, 0, 0, 10],
    [-10, 0, 0, -1],
    [0, -10, 1, 0]])
q_BE0 = np.array([1, 0, 0, 0]).T
lambda_ground_truth = -w_BE
modula = np.sqrt(np.sum(np.square(lambda_ground_truth)))
# q_ground_truth = np.array([cos(0.5*modula*pi/180), sin(0.5*modula*pi/180)*lambda_ground_truth[0]/modula, 
                           # sin(0.5*modula*pi/180)*lambda_ground_truth[1]/modula, sin(0.5*modula*pi/180)*lambda_ground_truth[2]/modula])
# print(q_ground_truth)
TM_ground_truth = rotation_vector_to_matrix(lambda_ground_truth)

def iter_TransformMatrix(initial_T=T_BE0, delta_t=0.01, tEnd=1):

    iter_times = int(tEnd / delta_t)
    T_current = initial_T
    for i in range(iter_times):
        delta_T = np.dot(T_current, omega_BE) * delta_t
        T_current = T_current + delta_T

    return T_current

def iter_RPY(initial_RPY=RPY_BE0.T, delta_t=0.01, tEnd=1):
    iter_times = int(tEnd / delta_t)
    rpy_current = initial_RPY
    for i in range(iter_times):
        r, p, y = rpy_current[0], rpy_current[1], rpy_current[2]
        delta_rpy = np.dot(np.array([
            [1, sin(r)*tan(p), cos(r)*tan(p)],
            [0, cos(r), -sin(r)],
            [0, sin(r)/cos(p), cos(r)/cos(p)]
            ]), w_BE.T) * delta_t
        rpy_current = rpy_current + delta_rpy
    
    return rpy_current


def iter_RotationVector(initial_rotation=q_BE0, delta_t=0.01, tEnd=1):
    iter_times = int(tEnd / delta_t)
    rv_current = initial_rotation
    for i in range(iter_times):
        delta_rv = 0.5 * np.dot(omega_no_ul.T, rv_current) * delta_t
        rv_current = rv_current + delta_rv
    return rv_current

TM_result = iter_TransformMatrix()
RPY_result = iter_RPY()
RV_result = iter_RotationVector()
print("Trans_Matrix_iter:\n ", TM_result)
print("RPY_iter:\n ", RPY_result)
print("RV_iter:\n ", RV_result)

def RPY2TM(rpy):
    r, p, y = rpy[0], rpy[1], rpy[2]
    trans_matrix = np.array([
        [cos(y)*cos(p), sin(y)*cos(p), -sin(p)],
        [cos(y)*sin(p)*sin(r)-sin(y)*cos(r), sin(y)*sin(p)*sin(r)+cos(y)*cos(r), cos(p)*sin(r)],
        [cos(y)*sin(p)*cos(r)+sin(y)*sin(r), sin(y)*sin(p)*cos(r)-cos(y)*sin(r), cos(p)*cos(r)]
    ])
    return trans_matrix

def RV2TM(rv):
    q0, q1, q2, q3 = rv[0], rv[1], rv[2], rv[3]
    trans_matrix = np.array([
        [q0**2+q1**2-q2**2-q3**2, 2*(q1*q2+q0*q3), 2*(q1*q3-q0*q2)],
        [2*(q1*q2-q0*q3), q0**2-q1**2+q2**2-q3**3, 2*(q2*q3+q0*q1)],
        [2*(q1*q3+q0*q2), 2*(q2*q3-q0*q1), q0**2-q1**2-q2**2+q3**2]
    ])
    return trans_matrix

def Lamda2TM(Lambda):
    lambda0, lambda1, lambda2 = Lambda[0], Lambda[1], Lambda[2]
    K_lambda = np.array([
        [0, -lambda2, lambda1],
        [lambda2, 0, -lambda0],
        [-lambda1, lambda0, 0]
        ])
    trans_matrix = np.exp(-K_lambda)
    return trans_matrix
    

TM_from_RPY = RPY2TM(RPY_result)
TM_from_RV = RV2TM(RV_result)
print("Trans_Matrix from RPY:\n ", TM_from_RPY)
print("Trans_Matrix from RV:\n ",TM_from_RV)
print("Ground Truth Trans_Matrix:\n ",TM_ground_truth)
    
def epsilon(A, B):
    C = A - B # the difference
    D = C*C.T # squaring the matrix
    e = np.max(np.linalg.eigvals(D))
    return np.sqrt(e)

error_matrix = epsilon(TM_result, TM_ground_truth)
error_RPY = epsilon(TM_from_RPY, TM_ground_truth)
error_RV = epsilon(TM_from_RV, TM_ground_truth)

print("Transformation Matrix Error:", error_matrix)
print("Euler Angles Error:", error_RPY)
print("Euler Symmetric Parameter Error:", error_RV)