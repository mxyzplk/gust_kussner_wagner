from eqmotion import Trimming1g, ACModel
import math
import numpy as np
from rk import Rk4
from kussner import Kussner
from wagner import Wagner
import os
import matplotlib.pyplot as plt
#
#
class Gust:
    def __init__(self):

        # 1g Trimming        
        self.ini = Trimming1g()
        self.mac = self.ini.ac.cw
        
        # Gust arrays        
        self.ude = self.ini.cond.ude
        self.tas = self.ini.cond.tas
        self.chords = self.ini.cond.chords
        self.g = self.chords * self.mac
        self.np = int()
        self.dt = 0.002

        fac = 10 ** (int(math.log10(1.2 * (self.chords * self.ini.ac.cw + self.ini.ac.xlew) / self.tas)) - 1)
        self.tlimit = int((1.2 * (self.chords * self.ini.ac.cw + self.ini.ac.xlew) / self.tas) + 1) * fac
        self.nps = int(self.tlimit / self.dt) + 2

        self.s = np.empty(self.nps)
        self.t = np.empty(self.nps)
        self.ug = np.empty(self.nps)
        
        # Motion parameters
        self.w = np.empty(self.nps)
        self.u = np.empty(self.nps)
        self.q = np.empty(self.nps)
        self.theta = np.empty(self.nps)
        self.wp = np.empty(self.nps)
        self.up = np.empty(self.nps)
        self.qp = np.empty(self.nps)
        self.thetap = np.empty(self.nps)
        self.aoawb = np.empty(self.nps)
        self.aoawbtot = np.empty(self.nps)
        self.aoaht = np.empty(self.nps)
        self.aoagwb = np.empty(self.nps)
        self.aoaght = np.empty(self.nps)
        self.nx = np.empty(self.nps)
        self.nz = np.empty(self.nps)        
        
        self.t[0] = 0.0
        self.ug[0] = 0.0
        self.s[0] = 0.0
        
        self.rkgwb = Kussner(self.tas, self.ini.ac.cw)
        self.rkght = Kussner(self.tas, self.ini.ac.cw)
        self.rkwwb = Wagner(self.tas, self.ini.ac.cw)
        self.rkwht = Wagner(self.tas, self.ini.ac.cw)
        y0 = [self.ini.w, self.ini.u, self.ini.q, self.ini.theta]
        self.rkm = Rk4(4, y0) # up wp qp thetap

        self.amodel = ACModel()
        
        for i in range(int(self.tlimit) - 1):
            self.t[i + 1] = self.t[i] + self.dt
            self.s[i + 1] = self.tas * self.t[i + 1] / self.mac
            self.ug[i + 1] = 0.5 * self.ude * (1 - math.cos((2 * math.pi * self.s[i + 1]) / self.g))
        
        self.delaywb = self.ini.ac.xlew / self.tas
        self.delayht = self.ini.ac.xleh / self.tas
        
        self.get_dirs()
        self.run_gust()
        self.plot_gust(self.aoawbtot, "AOA WB TOTAL", self.path_fig_gust1)
        self.plot_gust(self.aoawbtot, "AOA WB GUST", self.path_fig_gust2)
        self.plot_gust(self.aoawbtot, "AOA HT GUST", self.path_fig_gust3)
        self.plot_gust(self.aoawbtot, "NZ", self.path_fig_gust4)
        self.plot_gust(self.aoawbtot, "Q", self.path_fig_gust5)
        self.plot_gust(self.aoawbtot, "QP", self.path_fig_gust6)        

    def get_dirs(self):
        
        main_dir = os.path.dirname(os.path.abspath(__file__))
        res_dir = os.path.join(main_dir, '../results')
        self.path_gust = os.path.join(res_dir, 'gust_results.txt')
        self.path_fig_gust1 = os.path.join(res_dir, 'graph_aoa_wb_tot.png')
        self.path_fig_gust2 = os.path.join(res_dir, 'graph_aoa_wb_g.png')
        self.path_fig_gust3 = os.path.join(res_dir, 'graph_aoa_ht_g.png')
        self.path_fig_gust4 = os.path.join(res_dir, 'graph_nz.png')
        self.path_fig_gust5 = os.path.join(res_dir, 'graph_q.png')
        self.path_fig_gust6 = os.path.join(res_dir, 'graph_qp.png')
        
        with open(self.path_gust, 'w') as file:
            file.write("TIME         U            AOAWBTOT     AOAWB        AOAWB_G      AOAHT        AOAHT_G      NX           NZ           ")
            
        file.close()
        
        
    def run_gust(self):
        
        glwb = (self.amodel.np - self.ini.ac.m.cg) * self.mac
        glht = self.ini.ac.xlah - (self.inic.ac.xlew + self.amodel.np * self.mac)
        
        dydt = np.empty(4)
        
        self.w[0] = self.ini.w
        self.u[0] = self.ini.u
        self.q[0] = self.ini.q
        self.theta[0] = self.ini.theta
        self.wp[0] = self.ini.wp
        self.up[0] = self.ini.up
        self.qp[0] = self.ini.qp
        self.thetap[0] = self.ini.q
        self.aoawb[0] = self.ini.aoa
        self.aoawbtot[0] = self.ini.aoa
        self.aoaht[0] = self.ini.aoa_ht
        self.aoagwb[0] = 0.0
        self.aoaght[0] = 0.0         

        self.nx[0] = ((self.up[0] + self.w[0] * self.q[0]) / 9.80665) + math.sin(self.theta[0])
        self.nz[0] = ((self.wp[0] - self.u[0] * self.q[0]) / 9.80665) - math.cos(self.theta[0])

        for i in range(1, self.nps, 1):
            
            for j in range(4): # Runge Kutta
            
                if (j == 0):
                    trk = self.t[i]
                elif (j == 1):
                    trk = self.t[i] + 0.5 * self.dt
                elif (j == 2):
                    trk = self.t[i] + 0.5 * self.dt
                elif (j == 3):
                    trk = self.t[i] + self.dt
            
                if (trk >= self.delaywb):
                    timewb = trk - self.delaywb
                    ugwb = np.interp(timewb, self.t, self.ug)
                    awbg = math.atan(ugwb / self.tas)
                else:
                    timewb = 0.0
                    ugwb = 0.0
                    awbg = 0.0
                        
                if (trk >= self.delayht):
                    timeht = trk - self.delaywb
                    ught = np.interp(timeht, self.t, self.u)
                    ahtg = math.atan(ught / self.tas)    
                else:
                    timeht = 0.0
                    ught = 0.0
                    ahtg = 0.0
                    
                awbg1 = self.rkgwb.eval_rkstep(timewb, self.dt, j+1, awbg, self.tas)
                ahtg1 = self.rkght.eval_rkstep(timeht, self.dt, j+1, ahtg, self.tas)
                
                if (i == 0):
                    aqwb = 0
                    aqht = 0
                else:
                    aqwb = math.atan((self.w[i] - self.w[i-1]) / self.tas + (self.q[i] - self.q[i-1]) * glwb / self.tas)
                    aqht = math.atan((self.w[i] - self.w[i-1]) / self.tas + (self.q[i] - self.q[i-1]) * (glwb + glht) / self.tas)
                    
                awbw1 = self.rkwwb.eval_rkstep(timewb, self.dt, j+1, aqwb, self.tas)
                ahtw1 = self.rkwht.eval_rkstep(timeht, self.dt, j+1, aqht, self.tas)
                
                awbwag = awbw1 - aqwb
                ahtwag = ahtw1 - aqht
                
                self.aoagwb[i] = (awbg1 + awbwag) * 180 / math.pi
                self.aoaght[i] = (ahtg1 + ahtwag) * 180 / math.pi
                
                if (trk == 0):
                    pass
                else:
                    self.amodel.eval_db(self.aoawb[i-1], self.q[i-1], self.aoagwb[i], self.aoaght[i], self.ini.elev, self.tas, self.ini.cond.pdyn)
                    self.wp[i] = -9.80665 * math.cos(self.theta * math.pi / 180) + self.amodel.model.fz_b / self.ini.ac.m.weight
                    self.up[i] = -9.80665 * math.sin(self.theta * math.pi / 180) - self.amodel.model.fx_b / self.ini.ac.m.weight 
                    self.qp[i] = self.amodel.my_b / self.ini.ac.m.Iyy
                    self.thetap[i] = self.q[i]
                    
                    dydt[0] = self.wp[i]
                    dydt[1] = self.up[i]
                    dydt[2] = self.qp[i]
                    dydt[3] = self.thetap[i]
                    
                    self.rkm.rk(j+1, self.dt, self.t[i], dydt)
                    
                    self.w[i] = self.rkm.y[j + 1, 0]
                    self.u[i] = self.rkm.y[j + 1, 1]
                    self.w[i] = self.rkm.y[j + 1, 2]
                    self.theta[i] = self.rkm.y[j + 1, 3]
                    
                    if (j == 3):
                        self.aoawbtot[i] = self.aoawb[i] + self.aoagwb[i]
                        with open(self.path_gust, 'a') as file:
                            file.write("{:12.5f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f}".format(self.time[i], self.ug[i], self.aoawbtot[i], self.aoawb[i], self.aoagwb[i], self.aoaht[i], self.aoaght[i],  self.nx[i], self.nz[i]))
                            file.close()
                    
                    self.aoawb[i] = math.atan(self.w[i] / self.tas) * 180 / math.pi
                    self.aoaht[i] = self.amodel.aoaht
                    self.nx[i] = ((self.up[i] + self.w[i] * self.q[i]) / 9.80665) + math.sin(self.theta[i])
                    self.nz[i] = ((self.wp[i] - self.u[i] * self.q[i]) / 9.80665) - math.cos(self.theta[i])
                
                
    def plot_gust(self, yarray, label, file_path):
        plt.style.use('_mpl-gallery')
        
        # plot
        plt.plot(self.t[:], yarray)
        
        plt.xlabel('TIME')
        plt.ylabel(label)

        plt.savefig(file_path)                 

                    
                
                
                    