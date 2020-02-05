import numpy as np

def Create_Jacobian_Matrix(nodes_points, integral_point, universal_element):

    '''Calculating Jacobian matrix'''

    jacobian = np.zeros((4))

    elements_nodes_x = nodes_points[0]
    elements_nodes_y = nodes_points[1]
    
    jacobian[0] = np.dot(universal_element.dN_ksi[integral_point], elements_nodes_x)
    jacobian[1] = np.dot(universal_element.dN_ksi[integral_point], elements_nodes_y)  
    jacobian[2] = np.dot(universal_element.dN_eta[integral_point], elements_nodes_x) 
    jacobian[3] = np.dot(universal_element.dN_eta[integral_point], elements_nodes_y) 
    
    return jacobian

def Clalucate_Jacobian_Det(jacobian):

    '''Calculating Jacobian det'''

    jacobian_for_det = jacobian.reshape(2,2)

    return np.linalg.det(jacobian_for_det)

def Create_Jacobian_div_Det_Matrix(jacobian):

    jacobian_div_det_matrix = np.zeros((4))

    jacobian_det = Clalucate_Jacobian_Det(jacobian)

    jacobian_div_det_matrix[0] = jacobian[3]/jacobian_det
    jacobian_div_det_matrix[1] = jacobian[1]/jacobian_det
    jacobian_div_det_matrix[2] = jacobian[2]/jacobian_det
    jacobian_div_det_matrix[3] = jacobian[0]/jacobian_det

    return jacobian_div_det_matrix