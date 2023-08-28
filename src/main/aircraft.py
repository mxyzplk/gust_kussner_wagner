import os
from mass import Mass
#
#
class Aircraft:
    def __init__(self):
        self.m = Mass()
        
        self.sw = 0
        self.sh = 0
        self.cw = 0
        self.ch = 0
        self.xlew = 0
        self.xleh = 0
        self.zlew = 0
        self.zleh = 0
        self.xlaw = 0
        self.xlah = 0
        self.zlaw = 0
        self.zlah = 0        
        
        self.get_dirs()
        self.read_ac_data()
        
        
    def get_dirs(self):
        
        main_dir = os.path.dirname(os.path.abspath(__file__))
        resources_dir = os.path.join(main_dir, '../resources')
        self.path = os.path.join(resources_dir, 'aircraft.txt')


    def read_ac_data(self):
             
        with open(self.path, 'r') as file:
            
            line = file.readline()
            temp = line.split() 
            self.sw = float(temp[0])
            
            line = file.readline()
            temp = line.split()             
            self.sh = float(temp[0])

            line = file.readline()
            temp = line.split() 
            self.cw = float(temp[0])
            
            line = file.readline()
            temp = line.split()             
            self.ch = float(temp[0])

            line = file.readline()
            temp = line.split() 
            self.xlew = float(temp[0])
            self.zlew = float(temp[1])

            line = file.readline()
            temp = line.split() 
            self.xlaw = float(temp[0])
            self.zlaw = float(temp[1])

            line = file.readline()
            temp = line.split() 
            self.xleh = float(temp[0])
            self.zleh = float(temp[1])

            line = file.readline()
            temp = line.split() 
            self.xlah = float(temp[0])
            self.zlah = float(temp[1])
            
        file.close()
