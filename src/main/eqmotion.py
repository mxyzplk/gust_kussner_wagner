from aircraft import Aircraft
from condition import Condition
import math
from mass import Mass
import scipy.optimize
import os

class Trimming1g:
    def __init__(self):
        # Aircraft
        self.ac = Aircraft()
        
        # Model
        self.model = ACModel()
        
        # Controls
        self.elev = 0.0
        
        # Angles
        self.theta = 0.0
        self.aoa = 0.0
        self.aoa_ht = 0.0
        self.aoa_g = 0.0
        self.aoa_ht_g = 0.0
        self.u = 0.0
        self.w = 0.0
        self.q = 0.0
        self.up = 0.0
        self.wp = 0.0
        self.qp = 0.0
        
        # Condition
        self.cond = Condition()
        
        self.trimming()
        self.get_dirs()
        
        self.u = self.cond.tas * math.cos(self.aoa * math.pi / 180)
        self.w = self.cond.tas * math.sin(self.aoa * math.pi / 180)
        self.nx = ((self.up + self.w * self.q) / 9.80665) + math.sin(self.theta)
        self.nz = ((self.wp - self.u * self.q) / 9.80665) - math.cos(self.theta)            
        
        with open(self.path_trim, 'w') as file:
            file.write("[eq] aoa : {:.4f} | elev : {:.3f} | theta : {:.4f} | up : {:.3f} | wp : {:.3f} | qp : {:.3f} | nx : {:.3f} | nz : {:.3f}".format(self.aoa, self.elev, self.theta, self.up, self.wp, self.qp, self.nx, self.nz))

        file.close()
        

    def get_dirs(self):
        
        main_dir = os.path.dirname(os.path.abspath(__file__))
        res_dir = os.path.join(main_dir, '../results')
        self.path_trim = os.path.join(res_dir, 'trimming_results.txt')
        
    
    def trimming_func(self, x):
        self.model.eval_db(x[0], 0, 0, 0, x[1], self.cond.tas, self.cond.pdyn)
        wp = -9.80665 * math.cos(x[2] * math.pi / 180) + self.model.fz_b / self.ac.m.weight
        up = -9.80665 * math.sin(x[2] * math.pi / 180) - self.model.fx_b / self.ac.m.weight
        qp = self.model.my_b / self.ac.m.Iyy

        err = abs(wp) + abs(up) + abs(qp)

        return err    
    
    def trimming(self):

        ibounds = scipy.optimize.Bounds([-15, -25, -5],[25, 15, 25])      

        res = scipy.optimize.differential_evolution(self.trimming_func, ibounds, tol=1e-7, disp=False)
        
        res.x
        
        wp = -9.80665 * math.cos(res.x[2] * math.pi / 180) + self.model.fz_b / self.ac.m.weight
        up = -9.80665 * math.sin(res.x[2] * math.pi / 180) - self.model.fx_b / self.ac.m.weight
        qp = self.model.my_b / self.ac.m.Iyy   
        
        err = abs(wp) + abs(up) + abs(qp)
            
        #print("[eq] aoa : {:.4f} | elev : {:.3f} | theta : {:.4f} | up : {:.3f} | wp : {:.3f} | qp : {:.3f} | err : {:.5f}".format(res.x[0], res.x[1], res.x[2], up, wp, qp, err))
        
        self.aoa = res.x[0]
        self.elev = res.x[1]
        self.theta = res.x[2]
        self.up = up
        self.wp = wp
        self.qp = qp

class ACModel:
    def __init__(self):
        # Aircraft
        self.ac = Aircraft()
        
        # Mass
        self.mass = Mass()
        
        # Wind Axis = f(aoa)
        self.clwb_wa = 0.0
        self.cdwb_wa = 0.0
        self.cmwb_wa = 0.0
        self.clht_wa = 0.0
        self.cdht_wa = 0.0
        self.cmht_wa = 0.0   
        self.cltab_wa = 0.0
        self.cl_wa = 0.0
        self.cd_wa = 0.0
        self.cm_wa = 0.0
        self.mywb = 0.0
        self.myht = 0.0

        # Wind Axis = f(aoaht)  
        self.clht_wah = 0.0
        self.cle_wah = 0.0
        self.cdht_wah = 0.0

        # Distance between HT and Wing aerodynamic centers - Wind Axis
        self.dx_wa = 0.0
        self.dz_wa = 0.0
        
        # Distance between HT and Wing aerodynamic centers - Body Axis 
        self.dx_b = 0.0
        self.dz_b = 0.0
        
        # Distance between HT aerodynamic center and c.g. - Body Axis
        self.dxcght_b = 0.0

        # Distance between Wing aerodynamic center and c.g. - Body Axis        
        self.dxcgw_b = 0.0
        self.dzcgw_b = 0.0
        
        # Body Axis
        self.cl_b = 0.0
        self.cd_b = 0.0
        self.cm_b = 0.0
        self.fx_b = 0.0
        self.fz_b = 0.0
        self.my_b = 0.0
        
        # Angles
        self.eps = 0.0
        self.aoaht = 0.0
        self.etab = 0.0
        self.deda = 0.3122
        
        # Engine
        self.ftx_b = 21997.14
        self.myt_b = -9826.88
        
        # Paths
        self.path_db = None
        
        self.get_dirs()
        
        self.np = 0.25 - 0.006194872 / 0.090353983
        
        with open(self.path_db, 'w') as file:
            file.write("        AOA          AOAHT        ELEV         ETAB         DYNP        FXHT         FZHT         MYHT         FXE          FZE          MYE          FXWB         FZWB         MYWB         FXTAB        FZTAB        MYTAB        \n")
        
        file.close()        

    def get_dirs(self):
        
        main_dir = os.path.dirname(os.path.abspath(__file__))
        logs_dir = os.path.join(main_dir, '../logs')
        self.path_db = os.path.join(logs_dir, 'log_coefs.txt')
        
    def get_clwb_wind_axis(self, aoa, aoa_a):
        self.clwb_wa = 0.090353983 * (aoa + aoa_a) + 0.209164559

    def get_cdwb_wind_axis(self, aoa, aoa_a):
        self.cdwb_wa = 0.000445905 * (aoa + aoa_a) * (aoa + aoa_a) - 0.000069685 * (aoa + aoa_a) + 0.025997036

    def get_cmwb_wind_axis(self, aoa, aoa_a):
        self.cmwb_wa = 0.006194872 * (aoa + aoa_a) + 5.81195E-05

    def get_clht_aoaht_wind_axis_ht(self, aoaht):
        self.clht_wah = 0.073129296 * aoaht * self.ac.sh / self.ac.sw
    
    def get_clht_elev_wind_axis_ht(self, elev):
        self.cle_wah = 0.040317106 * elev * self.ac.sh / self.ac.sw
    
    def get_cdht_aoaht_wind_axis_ht(self, aoaht):
        self.cdht_wah = (0.000084918 * 0.000084918 * aoaht - 0.000050100 * aoaht - 0.000009588) * self.ac.sh / self.ac.sw
    
    def get_aoa_ht(self, aoa, q, aoa_ht_g, tas):
        self.eps = 0.3122 * aoa
        self.dxcght_b = self.ac.xlah - self.mass.x
        self.aoaht = aoa - self.eps + aoa_ht_g + self.dxcght_b * q * (180 / math.pi) / tas
        #print("[aoaht] aoa: {:10.4f} aoaht: {:10.4f} q: {:10.4f} eps: {:10.4f}".format(aoa, self.aoaht, q, self.eps))

    def get_clq(self, q, tas):
        qc2v = 0.5 * self.ac.cw * q / tas
        self.clq_s = 6.189331083 * qc2v
    
    def get_cmq(self, q, tas):
        qc2v = 0.5 * self.ac.cw * q / tas
        self.cmq_s = -32.50291511 * qc2v
        
    def get_cltab(self, tab):
        self.cltab_wa = 0.62 * tab * (self.ac.sh / self.ac.sw) * math.pi / 180
        self.etab = tab
        
    def get_wing_coefficients_wa(self, aoa, aoa_g, q, tas):
        self.get_clwb_wind_axis(aoa, aoa_g)
        self.get_cdwb_wind_axis(aoa, aoa_g)
        self.get_cmwb_wind_axis(aoa, aoa_g)
        self.get_clq(q, tas)
        self.get_cmq(q, tas)
               
    def get_ht_coefficients_wa(self, aoa, q, aoa_ht_g, tas, elev):
        self.get_aoa_ht(aoa, q, aoa_ht_g, tas)
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

        self.clht_wa1 = self.clht_wah * math.cos(self.eps * math.pi / 180) - self.cdht_wah * math.sin(self.eps * math.pi / 180)
        self.cdht_wa1 = self.clht_wah * math.sin(self.eps * math.pi / 180) + self.cdht_wah * math.cos(self.eps * math.pi / 180)

        self.cle_wa = self.cle_wah * math.cos(self.eps * math.pi / 180)
        self.cde_wa = self.cle_wah * math.sin(self.eps * math.pi / 180)

        self.cmht_wa1 = self.cdht_wa1 * self.dz_wa - self.clht_wa1 * self.dx_wa        
        
        # Distance between HT and Wing aerodynamic centers - Body Axis
        self.dx_b = (self.ac.xlah - self.ac.xlaw)
        self.dz_b = (self.ac.zlah - self.ac.zlaw)
        
        # Distance between HT and Wing aerodynamic centers - Wind Axis
        self.dx_wa = self.dx_b * math.cos(aoa * math.pi / 180) + self.dz_b * math.sin(aoa * math.pi / 180)
        self.dz_wa = -self.dx_b * math.sin(aoa * math.pi / 180) + self.dz_b * math.cos(aoa * math.pi / 180)
        
        # HT Pitch moment
        self.cmht_wa = self.cdht_wa * self.dz_wa - self.clht_wa * self.dx_wa
            
    def eval_db(self, aoa, q, aoa_g, aoa_ht_g, elev, tas, dynp):
        self.get_wing_coefficients_wa(aoa, aoa_g, q, tas)
        self.get_ht_coefficients_wa(aoa, q, aoa_ht_g, tas, elev)
        self.sum_coefficients_wa()
        self.sum_coefficients_b(aoa, elev, dynp)
        
    def sum_coefficients_wa(self):        
        self.cl_wa = self.clht_wa + self.clwb_wa
        self.cd_wa = self.cdht_wa + self.cdwb_wa
        self.cm_wa = self.cmwb_wa * self.ac.cw + self.cmht_wa + self.cmq_s * self.ac.cw
        
    def sum_coefficients_b(self, aoa, elev, pdyn):
        self.cd_b = -1.0 * self.cl_wa * math.sin(aoa * math.pi / 180) + self.cd_wa * math.cos(aoa * math.pi / 180)
        self.cl_b = 1.0 * self.cd_wa * math.sin(aoa * math.pi / 180) + self.cl_wa * math.cos(aoa * math.pi / 180) + self.clq_s
        # Cm at c.g.
        self.dxcgw_b = self.ac.xlaw - self.mass.x
        self.dzcgw_b = self.ac.zlaw - self.mass.z       
        self.cm_b = self.cm_wa - self.cl_b * self.dxcgw_b + self.cd_b * self.dzcgw_b
        # Forces and moments
        self.fx_b = self.cd_b * pdyn * self.ac.sw - self.ftx_b
        self.fz_b = self.cl_b * pdyn * self.ac.sw
        self.my_b = self.cm_b * pdyn * self.ac.sw + self.myt_b
        # Logs
        self.mywb = self.cmwb_wa * self.ac.cw * pdyn * self.ac.sw
        self.cdwb_g = -1.0 * self.clwb_wa * math.sin(aoa * math.pi / 180) + self.cdwb_wa * math.cos(aoa * math.pi / 180)
        self.clwb_g = 1.0 * self.cdwb_wa * math.sin(aoa * math.pi / 180) + self.clwb_wa * math.cos(aoa * math.pi / 180) + self.clq_s
        self.fzwb_g = self.clwb_g * pdyn * self.ac.sw
        self.fxwb_g = self.cdwb_g * pdyn * self.ac.sw
        self.cdht_g = -1.0 * self.clht_wa * math.sin(aoa * math.pi / 180) + self.cdht_wa * math.cos(aoa * math.pi / 180)
        self.clht_g = 1.0 * self.cdht_wa * math.sin(aoa * math.pi / 180) + self.clht_wa * math.cos(aoa * math.pi / 180)
        self.cdht_g1 = -1.0 * self.clht_wa1 * math.sin(aoa * math.pi / 180) + self.cdht_wa1 * math.cos(aoa * math.pi / 180)
        self.clht_g1 = 1.0 * self.cdht_wa1 * math.sin(aoa * math.pi / 180) + self.clht_wa1 * math.cos(aoa * math.pi / 180)        
        self.fzht_g = self.clht_g * pdyn * self.ac.sw
        self.fxht_g = self.cdht_g * pdyn * self.ac.sw
        self.cde_g = -1.0 * self.cle_wa * math.sin(aoa * math.pi / 180) + self.cde_wa * math.cos(aoa * math.pi / 180)
        self.cle_g = 1.0 * self.cde_wa * math.sin(aoa * math.pi / 180) + self.cle_wa * math.cos(aoa * math.pi / 180)
        self.etabcdb_g = -1.0 * self.cltab_wa * math.sin(aoa * math.pi / 180)
        self.etabclb_g = self.cltab_wa * math.cos(aoa * math.pi / 180)        
        self.fze_g = self.cle_g * pdyn * self.ac.sw
        self.fxe_g = self.cde_g * pdyn * self.ac.sw
        self.fztab_g = self.etabclb_g * pdyn * self.ac.sh
        self.fxtab_g = self.etabcdb_g * pdyn * self.ac.sh         
        self.fzht_g1 = self.clht_g1 * pdyn * self.ac.sw
        self.fxht_g1 = self.cdht_g1 * pdyn * self.ac.sw         
        self.dxcgh_b = self.ac.xlah - self.mass.x
        self.dzcgh_b = self.ac.zlah - self.mass.z          
        self.mywbcg = self.mywb + (- self.clwb_g * self.dxcgw_b + self.cdwb_g * self.dzcgw_b) * pdyn * self.ac.sw
        self.myecg = - self.fze_g * self.dxcgh_b + self.fxe_g * self.dzcgh_b
        self.mytabcg = - self.fztab_g * self.dxcgh_b + self.fxtab_g * self.dzcgh_b
        self.myhtcg = - self.fzht_g1 * self.dxcgh_b + self.fxht_g1  * self.dzcgh_b
        self.myht = self.cmht_wa  * pdyn * self.ac.sw
        
        with open(self.path_db, 'a') as file:
            file.write("{:12.3f} {:12.3f} {:12.3f} {:12.3f} {:12.0f} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e}\n".format(
                float(aoa),
                float(self.aoaht), 
                float(elev), 
                float(self.etab), 
                float(pdyn), 
                float(self.fxht_g1), 
                float(self.fzht_g1), 
                float(self.myhtcg),                 
                float(self.fxe_g), 
                float(self.fze_g), 
                float(self.myecg), 
                float(self.fxwb_g), 
                float(self.fzwb_g), 
                float(self.mywbcg), 
                float(self.fxtab_g), 
                float(self.fztab_g), 
                float(self.mytabcg)))                       
        
        file.close()
               
        
    def trimtab(self, elev, aoa_ht):
        b0 = 0.000
        b1 = -0.004 # aoa ht
        b2 = -0.390 # elev
        b3 = -0.353 # tab
        deg2rad = 180 / math.pi
        rad2deg = math.pi / 180
        tab = (-1.0 * (b0 + b1 * aoa_ht * deg2rad + b2 * elev * deg2rad) / b3) * rad2deg
        if (b0 + b1 * aoa_ht * deg2rad + b2 * elev * deg2rad >= 0):
            tab = min(tab, 20)
        else:
            tab = max(tab, -20)
            
        return tab
        
    
                

    
        
        
        


        








