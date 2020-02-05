import Functions_Data as fd 

import numpy as np

class UniversalElement():
    
    '''Universal Element of Grid'''

    def __init__(self):

        self.dN_ksi = np.zeros((4,4))
        self.dN_eta = np.zeros((4,4))
        self.N_matrix = np.zeros((4,4))

        self.dN_ksi_create()
        self.dN_eta_create()
        self.N_matrix_create()


    def dN_ksi_create(self): 

        integral_points1 = fd.integral_points1

        array = np.zeros((4))
        for i in range(4):

            if i < 2:
                eta = float(integral_points1[0][0])
            
            else: eta = float(integral_points1[1][0])

            for j in range(4):
                array[j] = fd.dN_dksi[j](eta)

            self.dN_ksi[i] = array
            
        return 0
    
    def dN_eta_create(self): 

        integral_points1 = fd.integral_points1

        array = np.zeros((4))
        for i in range(4):

            if i == 0 or i == 3:
                ksi = float(integral_points1[0][0])
            
            else: ksi = float(integral_points1[1][0])

            for j in range(4):
                array[j] = fd.dN_deta[j](ksi)

            self.dN_eta[i] = array
            
        
        return 0
    
    def N_matrix_create(self):

        integral_points1 = fd.integral_points1

        array = np.zeros((4))
        for i in range(4):

            if i == 0 or i == 3:
                ksi = float(integral_points1[0][0])

            else: ksi = float(integral_points1[1][0])

            if i < 2:
                eta = float(integral_points1[0][0])
            
            else: eta = float(integral_points1[1][0])

            for j in range(4):
                array[j] = fd.shape_function2D[j](ksi, eta)

            self.N_matrix[i] = array
        
        return 0  