import os
#

class Condition:
    def __init__(self):
        self.eas = 0
        self.tas = 0
        self.mach = 0
        self.ude = 0
        self.chords = 0
        self.altitude = 0
        self.path = None
        
        self.get_dirs()
        self.read_condition()

        self.pdyn = 0.5 * 1.225 * self.eas * self.eas        

    def get_dirs(self):
        
        main_dir = os.path.dirname(os.path.abspath(__file__))
        resources_dir = os.path.join(main_dir, '../resources')
        self.path = os.path.join(resources_dir, 'condition.txt')
        
    def read_condition(self):
             
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
            self.altitude = float(temp[0])

            line = file.readline()
            temp = line.split() 
            self.ude = float(temp[0])

            line = file.readline()
            temp = line.split() 
            self.chords = float(temp[0])
            
        file.close()
        