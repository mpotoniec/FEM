

class Node():
    
    '''Node of grid'''

    def __init__(self, node_number :int, x :float, y :float, boundary_condition :bool):
        
        self.node_number = node_number
        self.x = x
        self.y = y
        self.boundary_condition = boundary_condition #Wartość tylko 0 lub 1
        self.value = None
        self.value_name = 'None'

    
    def print_node(self):
        
        print('node number:',self.node_number)
        print('x =',self.x)
        print('y =',self.y)
        print('Boundary Condition =',self.boundary_condition)
        print(self.value_name + ':', self.value)


