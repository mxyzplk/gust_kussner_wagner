import os
#
#
class Mass:
    def __init__(self):
        self.weight = 0
        self.cg = 0
        self.x = 0
        self.y = 0
        self.z = 0
        self.Ixx = 0
        self.Iyy = 0
        self.Izz = 0
        self.Ixy = 0
        self.Iyz = 0
        self.Ixz = 0
        self.x = 0
        self.y = 0
        self.z = 0
        self.A = 0
        self.C = 0
        self.D = 0
        self.E = 0
        self.F = 0
        self.path = None
        
        self.get_dirs()
        self.read_inertia()
        
        #print("[mass] m : {:.1f} | ixx : {:.1f} | iyy : {:.1f} | izz : {:.1f} | ixz : {:.1f}".format(self.weight, self.Ixx, self.Iyy, self.Izz, self.Ixz))


    def get_dirs(self):
        
        main_dir = os.path.dirname(os.path.abspath(__file__))
        resources_dir = os.path.join(main_dir, '../resources')
        self.path = os.path.join(resources_dir, 'mass.txt')


    def read_inertia(self):
             
        with open(self.path, 'r') as file:
            
            line = file.readline()
            temp = line.split() 
            self.weight = float(temp[0])
            
            line = file.readline()
            temp = line.split()             
            self.cg = float(temp[0])

            line = file.readline()
            temp = line.split()             
            self.x = float(temp[0])
            
            line = file.readline()
            temp = line.split()             
            self.y = float(temp[0])

            line = file.readline()
            temp = line.split()             
            self.z = float(temp[0])            

            line = file.readline()
            temp = line.split() 
            self.Ixx = float(temp[0])
            
            line = file.readline()
            temp = line.split()             
            self.Iyy = float(temp[0])

            line = file.readline()
            temp = line.split() 
            self.Izz = float(temp[0])

            line = file.readline()
            temp = line.split() 
            self.Ixy = float(temp[0])

            line = file.readline()
            temp = line.split() 
            self.Iyz = float(temp[0])

            line = file.readline()
            temp = line.split() 
            self.Ixz = float(temp[0])
            
        file.close()    

        if ( self.Ixx != 0 and self.Izz != 0 and self.Ixz != 0 ):
            self.A = ((self.Iyy - self.Izz) / self.Ixx) - ((self.Ixz * self.Ixz) / (self.Ixx * self.Izz))
            self.C = ((self.Ixz * self.Ixz) / (self.Ixx * self.Izz)) + (self.Ixx - self.Iyy) / self.Izz
            self.D = (self.Ixz * (self.Iyy - self.Izz)) / (self.Ixx * self.Izz) - (self.Ixz / self.Izz)
            self.F = (self.Ixz * (self.Ixx - self.Iyy)) / (self.Ixx * self.Izz) + (self.Ixz / self.Ixx)
            self.E = 1.0 / (1.0 - ((self.Ixz * self.Ixz) / (self.Ixx * self.Izz)))

