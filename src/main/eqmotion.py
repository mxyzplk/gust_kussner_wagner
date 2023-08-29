from aircraft import Aircraft
from condition import Condition
import math
from mass import Mass

class Trimming1g:
    def __init__(self):
        # Aircraft
        self.ac = Aircraft()
        
        # Mass
        self.mass = Mass()
        
        # Controls
        self.elev = 0
        
        # Angles
        self.theta = 0
        self.aoa = 0
        self.aoa_ht = 0
        self.aoa_g = 0
        self.aoa_ht_g = 0
        
        # Condition
        self.cond = Condition()
        

class ACModel:
    def __init__(self):
        # Aircraft
        self.ac = Aircraft()
        
        # Mass
        self.mass = Mass()
        
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

        # Distance between HT and Wing aerodynamic centers - Wind Axis
        self.dx_wa = 0
        self.dz_wa = 0
        
        # Distance between HT and Wing aerodynamic centers - Body Axis 
        self.dx_b = 0
        self.dz_b = 0
        
        # Distance between HT aerodynamic center and c.g. - Body Axis
        self.dxcght_b = 0

        # Distance between Wing aerodynamic center and c.g. - Body Axis        
        self.dxcgw_b = 0
        self.dzcgw_b = 0
        
        # Body Axis
        self.cl_b = 0
        self.cd_b = 0
        self.cm_b = 0
        self.fx_b = 0
        self.fz_b = 0
        self.my_b = 0
        
        # Angles
        self.eps = 0
        self.aoaht = 0
        
        # Engine
        self.ftx_b = 21997.14
        self.myt_b = -9826.88
        
    def get_clwb_wind_axis(self, aoa, aoa_a):
        self.clwb_wa = 0.090353983 * (aoa + aoa_a) + 0.209164559

    def get_cdwb_wind_axis(self, aoa, aoa_a):
        self.cdwb_wa = 0.000445905 * (aoa + aoa_a) * (aoa + aoa_a) - 0.000069685 * (aoa + aoa_a) + 0.025997036

    def get_cmwb_wind_axis(self, aoa, aoa_a):
        self.cmwb_wa = 0.006194872 * (aoa + aoa_a) + 5.81195E-05

    def get_clht_aoaht_wind_axis_ht(self, aoaht):
        self.clht_wa = 0.073129296 * aoaht * self.ac.sh / self.ac.sw
    
    def get_clht_elev_wind_axis_ht(self, elev):
        self.cle_wa = 0.040317106 * elev * self.ac.sh / self.ac.sw
    
    def get_cdht_aoaht_wind_axis_ht(self, aoaht):
        self.cdht_wa = (0.000084918 * 0.000084918 * aoaht - 0.000050100 * aoaht - 0.000009588) * self.ac.sh / self.ac.sw
    
    def get_aoa_ht(self, aoa, q, aoa_ht_g, tas):
        self.eps = 0.3122 * aoa
        self.dxcght_b = self.ac.xlah - self.mass.x
        self.aoaht = aoa - self.eps + aoa_ht_g + self.dxcg_b * q * (180 / math.pi) / tas

    def get_clq(self, q, tas):
        qc2v = 0.5 * self.ac.cw * q / tas
        self.clq_s = 6.189331083 * qc2v
    
    def get_cmq(self, q, tas):
        qc2v = 0.5 * self.ac.cw * q / tas
        self.cmq_s = -32.50291511 * qc2v
        
    def get_cltab(self, tab):
        self.cltab_wa = 0.62 * tab * math.pi / 180
        
    def get_wing_coefficients_wa(self, aoa, aoa_g, q, tas):
        self.get_clwb_wind_axis(aoa, aoa_g)
        self.get_cdwb_wind_axis(aoa, aoa_g)
        self.get_cmwb_wind_axis(aoa, aoa_g)
        self.get_clq(q, tas)
        self.get_cmq(q, tas)
        
    def get_ht_coefficients_wa(self, aoa, q, aoa_ht_g, tas, elev, tab):
        self.get_aoa_ht(aoa, q, aoa_ht_g)
        tab = self.trimtab(elev, self.aoaht)
        self.get_cltab(tab)
        self.get_clht_aoaht_wind_axis_ht(self.aoaht)
        self.get_cdht_aoaht_wind_axis_ht(self.aoaht)
        self.get_clht_elev_wind_axis_ht(elev)
        
        cl = self.clht_wah + self.cle_wah
        cd = self.cdht_wah
        
        # HT aerodynamic coefficients - Wind Axis
        self.clht_wa = cl * math.cos(self.eps * math.pi / 180) - cd * math.sin(self.eps * math.pi / 180) + self.cltab_wa
        self.cdht_wa = cl * math.sin(self.eps * math.pi / 180) + cd * math.cos(self.eps * math.pi / 180)
        
        # Distance between HT and Wing aerodynamic centers - Body Axis
        self.dx_b = (self.ac.xlah - self.ac.xlaw)
        self.dz_b = (self.ac.zlah - self.ac.zlaw)
        
        # Distance between HT and Wing aerodynamic centers - Wind Axis
        self.dx_wa = self.dx_b * math.cos(aoa * math.pi / 180) + self.dz_b * math.sin(aoa * math.pi / 180)
        self.dz_wa = -1.0* self.dx_b * math.sin(aoa * math.pi / 180) + self.dz_b * math.sin(aoa * math.pi / 180)
        
        # HT Pitch moment
        self.cmht_wa = self.cdht_wa * self.dz_w - self.clht_wa * self.dx_wa
    
    def eval_db(self, aoa, q, aoa_g, aoa_ht_g, elev, tas, dynp, tab):
        self.get_wing_coefficients_wa(aoa, aoa_g, q, tas)
        self.get_ht_coefficients_wa(aoa, q, aoa_ht_g, tas, elev, tab)
        self.sum_coefficients_wa()
        self.sum_coefficients_b()
        
    def sum_coefficients_wa(self):        
        self.cl_wa = self.clht_wa + self.clwb_wa
        self.cd_wa = self.cdht_wa + self.cdwb_wa
        self.cm_wa = self.cmwb_wa + self.cmht_wa + self.cmq_s
        
    def sum_coefficients_b(self, aoa, pdyn):
        self.cl_b = -1.0 * self.cl_wa * math.sin(aoa * math.pi / 180) + self.cd_wa * math.cos(aoa * math.pi / 180) + self.clq_s
        self.cd_b = 1.0 * self.cd_wa * math.sin(aoa * math.pi / 180) + self.cl_wa * math.cos(aoa * math.pi / 180)
        # Cm at c.g.
        self.dxcgw_b = self.ac.xlaw - self.mass.x
        self.dzcgw_b = self.ac.zlaw - self.mass.z       
        self.cm_b = self.cmht_wa - self.cl_b * self.dxcgw_b + self.cd_b * self.dzcgw_b
        # Forces and moments
        self.fx_b = self.cd_b * pdyn * self.ac.sw + self.ftx_b
        self.fz_b = self.cl_b * pdyn * self.ac.sw
        self.my_b = self.cm_b * pdyn * self.ac.sw + self.myt_b
        
        print("[db] fxb : %10.1 | fzb : %10.1 | myb : %10.1")
        
    def trimtab(self, elev, aoa_ht):
        b0 = 0.000
        b1 = -0.004 # aoa ht
        b2 = -0.390 # elev
        b3 = -0.353 # tab
        deg2rad = 180 / math.pi
        rad2deg = math.pi / 180
        tab = (-1.0 * (b0 + b1 * aoa_ht * deg2rad + b2 * elev * deg2rad) / b3) * rad2deg
        if (b0 + b1 * aoa_ht * deg2rad + b2 * elev * deg2rad >= 0):
            tab = min(self.tab, 20)
        else:
            tab = max(self.tab, -20)
            
        return tab
        
    
        
        
    
        
        
        


        








