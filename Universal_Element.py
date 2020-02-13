import numpy as np

class UniversalElement():
    
    '''Universal Element of Grid with 4 matrices'''

    def __init__(self):

        #dN1/dksi, dN2/dksi, dN3/dksi, dN4/dksi 2D
        self.dN_dksi = (lambda eta: -0.25*(1-eta), \
                        lambda eta: 0.25*(1-eta),  \
                        lambda eta: 0.25*(1+eta),  \
                        lambda eta: -0.25*(1+eta))
        #dN1/deta, dN2/deta, dN3/deta, dN4/deta 2D
        self.dN_deta = (lambda ksi: -0.25*(1-ksi), \
                        lambda ksi: -0.25*(1+ksi), \
                        lambda ksi: 0.25*(1+ksi),  \
                        lambda ksi: 0.25*(1-ksi))
        #N1, N2, N3, N4 funkcje kształtu 2D
        self.shape_function2D = (lambda ksi, eta: 0.25*(1-ksi)*(1-eta), \
                                 lambda ksi, eta: 0.25*(1+ksi)*(1-eta), \
                                 lambda ksi, eta: 0.25*(1+ksi)*(1+eta), \
                                 lambda ksi, eta: 0.25*(1-ksi)*(1+eta))
        #Punkty całkowania dla dwupunktowego schematu:
        self.integral_points = [[-1./(3.**0.5), 1.], [1./(3.**0.5), 1.]]

        self.dN_ksi_matrix = np.zeros((4,4))
        self.dN_eta_matrix = np.zeros((4,4))
        self.N_matrix = np.zeros((4,4))

        self.dN_ksi_create()
        self.dN_eta_create()
        self.N_matrix_create()

    def dN_ksi_create(self): 

        array = np.zeros((4))
        for i in range(4):

            if i < 2:
                eta = float(self.integral_points[0][0])
            
            else: eta = float(self.integral_points[1][0])

            for j in range(4):
                array[j] = self.dN_dksi[j](eta)

            self.dN_ksi_matrix[i] = array
    
    def dN_eta_create(self): 

        array = np.zeros((4))
        for i in range(4):

            if i == 0 or i == 3:
                ksi = float(self.integral_points[0][0])
            
            else: ksi = float(self.integral_points[1][0])

            for j in range(4):
                array[j] = self.dN_deta[j](ksi)

            self.dN_eta_matrix[i] = array
    
    def N_matrix_create(self):

        array = np.zeros((4))
        for i in range(4):

            if i == 0 or i == 3:
                ksi = float(self.integral_points[0][0])

            else: ksi = float(self.integral_points[1][0])

            if i < 2:
                eta = float(self.integral_points[0][0])
            
            else: eta = float(self.integral_points[1][0])

            for j in range(4):
                array[j] = self.shape_function2D[j](ksi, eta)

            self.N_matrix[i] = array