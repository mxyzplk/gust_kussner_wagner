import os
#
#
class Gust:
    def __init__(self, mac, a0):
        
        # Gust data & initial condition
        self.ude = 0
        self.tas = 0
        self.eas = 0
        self.mach = 0
        self.h = 0
        self.nchords = 0
        self.mac = mac
        self.G = 0
        self.a0 = a0

        # Gust arrays        
        self.s = None
        self.t = None
        self.u = None
        self.np = 0
        self.dt = 0
        
    def get_dirs(self):
        
        main_dir = os.path.dirname(os.path.abspath(__file__))
        resources_dir = os.path.join(main_dir, '../resources')
        self.path = os.path.join(resources_dir, 'condition.txt')


    def read_ac_data(self):
             
        with open(self.path, 'r') as file:
            
            line = file.readline()
            temp = line.split() 
            self.eas = float(temp[0])
            
            line = file.readline()
            temp = line.split()             
            self.tas = float(temp[0])

            line = file.readline()
            temp = line.split() 
            self.mach = float(temp[0])
            
            line = file.readline()
            temp = line.split()             
            self.h = float(temp[0])

            line = file.readline()
            temp = line.split() 
            self.ude = float(temp[0])

            line = file.readline()
            temp = line.split() 
            self.nchords = float(temp[0])
            
        self.G = self.nchords * self.mac
            
        file.close()
