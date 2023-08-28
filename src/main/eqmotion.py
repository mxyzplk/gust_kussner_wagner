from aircraft import Aircraft
import math

class ACModel:
    def __init__(self):
        # Aircraft
        self.ac = Aircraft()
        
        # Wind Axis = f(aoa)
        self.clwb_wa = 0
        self.cdwb_wa = 0
        self.cmwb_wa = 0
        self.clht_wa = 0
        self.cdht_wa = 0
        self.cmht_wa = 0   
        self.cltab_wa = 0
        self.cl_wa = 0
        self.cd_wa = 0
        self.cm_wa = 0

        # Wind Axis = f(aoaht)  
        self.clht_wah = 0
        self.cle_wah = 0 
        self.cdht_wah = 0

        # Arms - Wind Axis
        self.dx_wa = 0
        self.dz_wa = 0
        
        self.eps = 0
        self.aoaht = 0
        self.clq_s = 0
        self.cmq_s = 0
        
        self.ftx = 21997.14
        self.myt = -9826.88
        
    def get_clwb_wind_axis(self, aoa):
        self.clwb_wa = 0.090353983 * aoa + 0.209164559

    def get_cdwb_wind_axis(self, aoa):
        self.cdwb_wa = 0.000445905 * aoa * aoa + -0.000069685 * aoa + 0.025997036

    def get_cmwb_wind_axis(self, aoa):
        self.cmwb_wa = 0.006194872 * aoa + 5.81195E-05

    def get_clht_aoaht_wind_axis_ht(self, aoaht):
        self.clht_wa = 0.073129296 * aoaht * self.ac.sh / self.ac.sw
    
    def get_clht_elev_wind_axis_ht(self, elev):
        self.cle_wa = 0.040317106 * elev * self.ac.sh / self.ac.sw
    
    def get_cdht_aoaht_wind_axis_ht(self, aoaht):
        self.cdht_wa = (0.000084918 * 0.000084918 * aoaht - 0.000050100 * aoaht - 0.000009588) * self.ac.sh / self.ac.sw
    
    def get_aoa_ht(self, aoa, q):
        self.eps = 0.3122 * aoa
        self.aoaht = aoa - self.eps

    def get_clq(self, q, tas):
        qc2v = 0.5 * self.ac.cw * q / tas
        self.clq_s = 6.189331083 * qc2v
    
    def get_cmq(self, q, tas):
        qc2v = 0.5 * self.ac.cw * q / tas
        self.cmq_s = -32.50291511 * qc2v
        
    def get_cltab(self, tab):
        self.cltab_wa = 0.62 * tab * math.pi / 180
        
    def get_clht_wa(self, aoa):
        cl = self.clht_wah + self.cle_wah
        cd = self.cdht_wah
        
        self.clht_wa = cl * math.cos(self.eps * math.pi / 180) - cd * math.sin(self.eps * math.pi / 180) + self.cltab_wa
        self.cdht_wa = cl * math.sin(self.eps * math.pi / 180) + cd * math.cos(self.eps * math.pi / 180)
        
        # Body Axis
        dx = (self.ac.xlah - self.ac.xlaw)
        dz = (self.ac.zlah - self.ac.zlaw)
        
        self.dx_wa = dx * 
        self.dz_wa = 
        
        self.cmht_wa = self.cdht_wa * (self.ac.zlah - self.ac.zlaw) - self.clht_wa
        
    def get_coefs_wa(self):        
        self.cl_wa = self.clht_wa + self.clwb_wa
        self.cd_wa = self.cdht_wa + self.cdwb_wa
        self.cm_wa = self.cmwb_wa + self.cmht_wa
        


        








