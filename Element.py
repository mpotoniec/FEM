from Functions_Data import shape_function1D, shape_function2D, integral_points1
from Jacobian import Create_Jacobian_Matrix, Clalucate_Jacobian_Det, Create_Jacobian_div_Det_Matrix

import numpy as np

class Element():

    def __init__(self, element_number :int):
        
        self.element_number = element_number
        self.nodes_array = []

        self.H_boundary_matrix = np.zeros((4,4))
        self.H_matrix = np.zeros((4,4))
        self.C_matrix = np.zeros((4,4))
        self.P_vector = np.zeros((4,1))

    def sort_array(self):

        tmp_node = self.nodes_array[1]
        self.nodes_array[1] = self.nodes_array[2]
        self.nodes_array[2] = tmp_node

        tmp_node = self.nodes_array[3]
        self.nodes_array[3] = self.nodes_array[2]
        self.nodes_array[2] = tmp_node

    def print_element(self):
        
        print('element number:',self.element_number)
        
        self.nodes_array[0].print_node()
        self.nodes_array[3].print_node()
        self.nodes_array[1].print_node()
        self.nodes_array[2].print_node()

    def get_nodes_points(self):

        x = []
        y = []

        for node in self.nodes_array:
            
            x.append(node.x)
            y.append(node.y)

        return x, y
    
    def calculate_edges_with_boundary_condition(self):

        edges_with_boundary_conditions = [None for i in range(4)]

        for i in range(4):
            if self.nodes_array[i].boundary_condition and self.nodes_array[(i + 1)%4].boundary_condition:
                edges_with_boundary_conditions[i] = [self.nodes_array[i], self.nodes_array[(i + 1)%4]]

        return edges_with_boundary_conditions

    def calculate_side_length(self, node1, node2): 
        return abs(((node2.x - node1.x)**2 + (node2.y - node1.y)**2)**0.5)

    def H_boundary_matrix_create(self, alpha):

        fields_matrices = []

        edges_with_boundary_condition = self.calculate_edges_with_boundary_condition()
        jacobian_det_matrix = np.zeros((4))

        side_length = np.zeros((4))

        for i in range(4):
            if edges_with_boundary_condition[i] != None:
                side_length[i] = self.calculate_side_length(edges_with_boundary_condition[i][0], edges_with_boundary_condition[i][1])

        for integral_point in range(4):
            jacobian_det_matrix[integral_point] = side_length[integral_point]/2

        if edges_with_boundary_condition[0] != None: #Powierzchnia 1
            field1_N = []
            for integral_point in range(2):
                ksi = integral_points1[integral_point][0]
                eta = integral_points1[integral_point][1]*(-1)

                N1 = shape_function2D[0](ksi, eta)
                N2 = shape_function2D[1](ksi, eta)
        
                tmp_array = np.zeros((4,4))

                tmp_array[0][0] = N1*N1*alpha
                tmp_array[0][1] = N1*N2*alpha
                tmp_array[1][0] = N2*N1*alpha
                tmp_array[1][1] = N2*N2*alpha

                field1_N.append(tmp_array)

            field1 = field1_N[0]+field1_N[1]
            field1 = field1*jacobian_det_matrix[0]
            fields_matrices.append(field1)

        if edges_with_boundary_condition[1] != None: #Powierzchnia 2
            field2_N = []
            for integral_point in range(2):
                ksi = integral_points1[integral_point][1]
                eta = integral_points1[integral_point][0]

                N2 = shape_function2D[1](ksi, eta)
                N3 = shape_function2D[2](ksi, eta)
        
                tmp_array = np.zeros((4,4))

                tmp_array[1][1] = N2*N2*alpha
                tmp_array[1][2] = N2*N3*alpha
                tmp_array[2][1] = N3*N2*alpha
                tmp_array[2][2] = N3*N3*alpha

                field2_N.append(tmp_array)

            field2 = field2_N[0]+field2_N[1]
            field2 = field2*jacobian_det_matrix[1]
            fields_matrices.append(field2)

        if edges_with_boundary_condition[2] != None: #Powierzchnia 3
            field3_N = []
            for integral_point in range(2):
                ksi = integral_points1[1 - integral_point][0]
                eta = integral_points1[integral_point][1]

                N3 = shape_function2D[2](ksi, eta)
                N4 = shape_function2D[3](ksi, eta)
        
                tmp_array = np.zeros((4,4))

                tmp_array[2][2] = N3*N3*alpha
                tmp_array[2][3] = N3*N4*alpha
                tmp_array[3][2] = N4*N3*alpha
                tmp_array[3][3] = N4*N4*alpha

                field3_N.append(tmp_array)

            field3 = field3_N[0]+field3_N[1]
            field3 = field3*jacobian_det_matrix[2]
            fields_matrices.append(field3)

        if edges_with_boundary_condition[3] != None: #Powierzchnia 4
            field4_N = []
            for integral_point in range(2):
                ksi = integral_points1[integral_point][1]*(-1)
                eta = integral_points1[1 - integral_point][0] 

                N1 = shape_function2D[0](ksi, eta)
                N4 = shape_function2D[3](ksi, eta)
        
                tmp_array = np.zeros((4,4))

                tmp_array[0][0] = N1*N1*alpha
                tmp_array[0][3] = N1*N4*alpha
                tmp_array[3][0] = N4*N1*alpha
                tmp_array[3][3] = N4*N4*alpha

                field4_N.append(tmp_array)

            field4 = field4_N[0]+field4_N[1]
            field4 = field4*jacobian_det_matrix[3]
            fields_matrices.append(field4)

        for miatrix in fields_matrices:
            self.H_boundary_matrix +=  miatrix

    def H_matrix_create(self, k, universal_element):

        dN_dx = np.zeros((4,4))
        dN_dy = np.zeros((4,4))

        dN_dxT = []
        dN_dyT = []
        dN_dxdy = []

        jacobian_det_matrix = np.zeros((4))
        jacobian_div_det_matrix = np.zeros((4,4))
    
        for integral_point in range(4):
            jacobian_det_matrix[integral_point] = Clalucate_Jacobian_Det(Create_Jacobian_Matrix(self.get_nodes_points(), integral_point, universal_element))
            jacobian_div_det_matrix[integral_point] = Create_Jacobian_div_Det_Matrix(Create_Jacobian_Matrix(self.get_nodes_points(), integral_point, universal_element))

        #Tworzenie dN/dx i dN/dy
        for integral_point in range(4):
        
            for i in range(4):

                dN_dx[integral_point][i] = universal_element.dN_ksi[integral_point][i]*jacobian_div_det_matrix[integral_point][0] + universal_element.dN_eta[integral_point][i]*jacobian_div_det_matrix[integral_point][1] 
                dN_dy[integral_point][i] = universal_element.dN_ksi[integral_point][i]*jacobian_div_det_matrix[integral_point][2] + universal_element.dN_eta[integral_point][i]*jacobian_div_det_matrix[integral_point][3]

        #Tworzenie {dN/x}{dN/dx}T {dN/y}{dN/dy}T
        for integral_point in range(4):
        
            el_arr_x = np.zeros((4,4))
            el_arr_y = np.zeros((4,4))

            for i in range(4):
                for j in range(4):
                
                    el_arr_x[i][j] = dN_dx[integral_point][i]*dN_dx[integral_point][j]
                    el_arr_y[i][j] = dN_dy[integral_point][i]*dN_dy[integral_point][j]

            dN_dxT.append(el_arr_x)
            dN_dyT.append(el_arr_y)

        #Mnożenie dN/dx i dN/dy przez wyznaczkin macierzy jakobianu
        for integral_point in range(4):

            dN_dxT[integral_point] = dN_dxT[integral_point]*jacobian_det_matrix[integral_point]
            dN_dyT[integral_point] = dN_dyT[integral_point]*jacobian_det_matrix[integral_point]

        for integral_point in range(4):
            dN_dxdy.append(k*(dN_dxT[integral_point] + dN_dyT[integral_point]))

        for integral_point in range(4):
            self.H_matrix += dN_dxdy[integral_point]

    def C_matrix_create(self, c, ro, universal_element):
    
        jacobian_det_matrix = np.zeros((4))

        for integral_point in range(4):
            jacobian_det_matrix[integral_point] = Clalucate_Jacobian_Det(Create_Jacobian_Matrix(self.get_nodes_points(), integral_point, universal_element))

        def Create_C_Matrix_For_One_Integral_Point(jacobian_det, c, ro, integral_point, universal_element):
        
            C_matrix_for_one_integral_point = np.zeros((4,4))

            for i in range(4):
                for j in range(4):
                    C_matrix_for_one_integral_point[i][j] = c*ro* jacobian_det* universal_element.N_matrix[integral_point][i]*universal_element.N_matrix[integral_point][j]

            return C_matrix_for_one_integral_point
    
        for integral_point in range(4):
            self.C_matrix += Create_C_Matrix_For_One_Integral_Point(jacobian_det_matrix[integral_point], c, ro, integral_point, universal_element)

    def P_vector_create(self, alpha, t_ambient):

        edges_with_boundary_condition = self.calculate_edges_with_boundary_condition()
        jacobian_det_matrix = np.zeros((4))

        side_length = np.zeros((4))

        for i in range(4):
            if edges_with_boundary_condition[i] != None:
                side_length[i] = self.calculate_side_length(edges_with_boundary_condition[i][0], edges_with_boundary_condition[i][1])

        for integral_point in range(4):
            jacobian_det_matrix[integral_point] = side_length[integral_point]/2

        N = np.zeros(shape = (2, 1))
        for i in range(2):
            N += np.array([[shape_function1D[0](integral_points1[i][0]), shape_function1D[1](integral_points1[i][0])]]).reshape(2, 1)


        for i in range(4):
            if edges_with_boundary_condition[i] != None:
                if i != 3:
                    self.P_vector[i] += N[0]* jacobian_det_matrix[i]* alpha* t_ambient
                    self.P_vector[i + 1] += N[1]* jacobian_det_matrix[i]* alpha* t_ambient
                else:
                    self.P_vector[i] += N[1]* jacobian_det_matrix[i]* alpha* t_ambient 
                    self.P_vector[i - 3] += N[0]* jacobian_det_matrix[i]* alpha* t_ambient

        self.P_vector*=(-1)   