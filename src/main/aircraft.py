class Aircraft:
    def __init__(self):
        self.dclda     = None
        self.Sw        = None
        self.mac       = None


    def set_ac_parameters(self, Sw, mac):
        self.Sw        = Sw
        self.mac       = mac


    def set_aero(self, dclda):
        self.dclda     = dclda