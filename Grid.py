from Element import Element
from Node import Node
from Universal_Element import UniversalElement

import time

import numpy as np

class Grid():

    def __init__(self, filename = 'variables.txt'):

        self.filename = filename

        variables = self.read_variables(self.filename)

        self.high = float(variables[0][1]) #Wysokość siatki (H)
        self.width = float(variables[1][1]) #Szerokość siatki (W)
        self.nodes_in_high = int(variables[2][1])  #liczba węzłow na wysokość
        self.nodes_in_width = int(variables[3][1]) #Liczba węzłów na długość
        
        self.nodes_number = int(self.nodes_in_high*self.nodes_in_width)
        self.elements_number = int((self.nodes_in_high - 1)*(self.nodes_in_width - 1))

        self.nodes = []
        self.elements = []
        
        self.H_global_matrix = np.zeros((self.nodes_number, self.nodes_number))
        self.C_global_matrix = np.zeros((self.nodes_number, self.nodes_number))
        self.H_boundary_global_matrix = np.zeros((self.nodes_number, self.nodes_number))
        self.P_global_vector = np.zeros((self.nodes_number, 1))

        self.universal_element = UniversalElement()

        self.create_nodes()
        self.create_elements()
        self.sort_elements_node_arrays()
        self.add_values((variables[4][0]).split(' ')[1], variables[4][1])

    def read_variables(self, filename):

        file = open(filename, 'r')

        from_file = []

        for line in file:

            line = line.split('\n')
            line = line[0].split('=')
            line[1] = float(line[1])

            from_file.append(line)

        file.close()

        return from_file

    def create_nodes(self):

        node_number = int(1)

        y = 0
        x = 0

        dy = self.high/ (self.nodes_in_high - 1) #todo
        dx = self.width/ (self.nodes_in_width - 1) #todo

        for width in range(0, self.nodes_in_width, 1):

            y=0

            for high in range(0, self.nodes_in_high, 1):
            
                #Blok if-else dotyczy warunków brzegowych w siatce
                if((high == 0) or (high == self.nodes_in_high - 1)): k = 1
                elif((width == 0) or width == self.nodes_in_width -1): k = 1
                else: k = 0

                tmp_node = Node(node_number, x, y, k)
                self.nodes.append(tmp_node)

                node_number+=1

                y+=dy
            
            x+=dx

    def create_elements(self): #todo

        element_number = int(1)
        high = 1

        a = -1
        b = self.nodes_in_high

        for element in range(self.elements_number): #Pętla w Po liczbie elementów mE max liczba elementów. Pierwszy element = 0.
            
            self.elements.append(Element(element_number))

            if high == self.nodes_in_high: #Warunek czy doszliśmy do końca wysokości.
                a+=1
                b+=1
                high = 1

            for j in range(4): #Pętla po 4 nodach w elemencie

                if((j == 0) or (j == 1)): self.elements[element].nodes_array.append(self.nodes[j + element_number + a])
                else: self.elements[element].nodes_array.append(self.nodes[j + element_number + b - 3])

            element_number+=1
            high+=1

    def sort_elements_node_arrays(self):
        
        for element in self.elements:
            element.sort_array()

    def add_values(self, value_name, value):

        for node in self.nodes:
            node.value_name = value_name
            node.value = value

    def print_grid(self):

        for element in self.elements:
            element.print_element()
            print('######################################################################')

    def H_global_matrix_create(self, k):
    
        for element in self.elements:
        
            element.H_matrix_create(k, self.universal_element)

            nodes_ids = [node.node_number - 1 for node in element.nodes_array]
            for row in range(4):
                for column in range(4):
                    self.H_global_matrix[nodes_ids[row], nodes_ids[column]] += element.H_matrix[row, column]

    def C_global_matrix_create(self, c, ro):

        for element in self.elements:

            element.C_matrix_create(c, ro, self.universal_element)

            nodes_ids = [node.node_number - 1 for node in element.nodes_array]

            for row in range(4):
                for column in range(4):
                    self.C_global_matrix[nodes_ids[row], nodes_ids[column]] += element.C_matrix[row, column]

    def P_global_vector_create(self, alpha, t_ambient):
        
        '''This grid function is resposible of calculating global P vector of grid.
            It takes a alpha coefficient and ambient temperature.'''
        
        for element in self.elements:

            element.P_vector_create(alpha, t_ambient)

            nodes_ids = [node.node_number - 1 for node in element.nodes_array]

            for row in range(4):
                self.P_global_vector[nodes_ids[row]] += element.P_vector[row]

    def H_boundary_global_matrix_create(self, alpha):

        for element in self.elements:

            element.H_boundary_matrix_create(alpha)

            nodes_ids = [node.node_number - 1 for node in element.nodes_array]

            for row in range(4):
                for column in range(4):
                    self.H_boundary_global_matrix[nodes_ids[row], nodes_ids[column]] += element.H_boundary_matrix[row, column] 

    def calculate_transform(self, simulation_time=None,simulation_step=None,k=None,c=None,ro=None,alpha=None,t_ambient=None,inital_temperature=None): #todo

        variables = self.read_variables(self.filename)

        if simulation_time == None: simulation_time = float(variables[5][1])
        if simulation_step == None: simulation_step = float(variables[6][1])
        if alpha == None: alpha = float(variables[7][1])
        if k == None: k = float(variables[8][1])
        if c == None: c = float(variables[9][1])
        if ro == None: ro = float(variables[10][1])
        if t_ambient == None: t_ambient = float(variables[11][1])
        if inital_temperature == None: inital_temperature = float(variables[4][1])

        steps = int(simulation_time/simulation_step)

        self.H_global_matrix_create(k)
        self.C_global_matrix_create(c, ro)
        self.H_boundary_global_matrix_create(alpha)
        self.P_global_vector_create(alpha, t_ambient)

        self.P_global_vector*=(-1)

        file = open('results.txt', 'w')

        dT = simulation_step
        for step in range(steps):

            temperatures = np.zeros((self.nodes_number))
            for i in range(self.nodes_number):
                temperatures[i] = self.nodes[i].value

            #Calculating {[H]+[H_BC]+[C]/dT}:
            H = self.H_global_matrix + self.H_boundary_global_matrix + self.C_global_matrix/dT

            #Wektor {[P]+{[C]/dT}*T0}:
            P = np.zeros((self.nodes_number, 1))

            for i in range(self.nodes_number):

                C_columns_tmp = 0
                for j in range(self.nodes_number):
                   
                    C_columns_tmp += self.C_global_matrix[i][j] /dT * temperatures[j]
                    
                P[i] = self.P_global_vector[i] + C_columns_tmp

            #Tworzenie macierzy H + P (nodes x nodes + 1)
            addedMatrix = np.zeros((self.nodes_number, self.nodes_number + 1))

            for i in range(self.nodes_number):
                addedMatrix[i][self.nodes_number] = P[i]

                for j in range(self.nodes_number):
                    addedMatrix[i][j] = H[i][j]

            #Metoda Eliminacji Gausa
            temp = np.zeros((self.nodes_number + 1))
            k = int(0)
            for i in range(self.nodes_number - 1):
                
                k+=1
                
                for l in range(k):
                    
                    for j in range(self.nodes_number + 1):
                        temp[j] = addedMatrix[l][j]*addedMatrix[k][l]/addedMatrix[l][l]

                    for j in range(self.nodes_number + 1):
                        addedMatrix[k][j] -= temp[j] 

            for i in range(self.nodes_number, 0, -1):

                sub = addedMatrix[i - 1, self.nodes_number]
                for j in range(self.nodes_number, i, -1):
                    sub -= addedMatrix[i - 1][j - 1]*temperatures[j - 1]
                
                temperatures[i - 1] = sub/addedMatrix[i - 1, i - 1]

            for i in range(self.nodes_number):
                self.nodes[i].value = temperatures[i]

            max_temperature = "%.2f" %temperatures.max()
            min_temperature = "%.2f" %temperatures.min()

            iteration_info = str('Iteration: ' + str(step + 1) + ' Time: ' + str(int(dT*(step + 1)))+'/'+str(int(simulation_time)))
            iteration_result = str('Minimum temperature = ' + str(min_temperature) + ' Maximum temperature = ' + str(max_temperature))

            print(iteration_info)
            print(iteration_result)

            file.write(iteration_info)
            file.write('\n')
            file.write(iteration_result)
            file.write('\n')

        file.close()            