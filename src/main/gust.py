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
        self.tlimit = ((1.2 * (self.chords * self.ini.ac.cw + self.ini.ac.xlew) / self.tas / fac) + 1) * fac
        self.nps = int(self.tlimit / self.dt) + 2

        self.s = np.empty(self.nps)
        self.t = np.empty(self.nps)
        self.ug = np.empty(self.nps)
        
        # Equations of Motion
        self.w = np.empty(self.nps)
        self.u = np.empty(self.nps)
        self.q = np.empty(self.nps)
        self.theta = np.empty(self.nps)
        self.wp = np.empty(self.nps)
        self.up = np.empty(self.nps)
        self.qp = np.empty(self.nps)
        self.thetap = np.empty(self.nps)
        self.nx = np.empty(self.nps)
        self.nz = np.empty(self.nps)    
        
        # AoA Results
        self.aoawb = np.empty(self.nps)
        self.aoawbtot = np.empty(self.nps)
        self.aoaht = np.empty(self.nps)
        self.aoagwb = np.empty(self.nps)   # AoA Gust - Wing
        self.aoaght = np.empty(self.nps)   # AoA Gust - HT

        # Intermediate parameters
        self.aoa_g = np.empty(self.nps)
        self.aoa_g_ht = np.empty(self.nps)
        self.aoa_kussner_wb = np.empty(self.nps)
        self.aoa_wagner_wb = np.empty(self.nps)        
        self.aoa_kussner_ht = np.empty(self.nps)
        self.aoa_wagner_ht = np.empty(self.nps)
        self.aoa_wagner_wb2 = np.empty(self.nps)         
        self.aoa_wagner_ht2 = np.empty(self.nps)        
        self.time_wb = np.empty(self.nps)
        self.time_ht = np.empty(self.nps)  
        self.aoa_w = np.empty(self.nps)
        self.aoa_w_ht = np.empty(self.nps)    
        
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
        
        for i in range(int(self.nps) - 1):
            self.t[i + 1] = self.t[i] + self.dt
            self.s[i + 1] = self.tas * self.t[i + 1] / self.mac
            self.ug[i + 1] = 0.5 * self.ude * (1 - math.cos((2 * math.pi * self.s[i + 1]) / self.chords))
            
        self.delaywb = self.ini.ac.xlew / self.tas
        self.delayht = self.ini.ac.xleh / self.tas
        
        self.get_dirs()
        self.run_gust()
        self.plot_gust(self.aoawbtot[0:self.nps-1], "AOA WB TOTAL", self.path_fig_gust1)
        self.plot_gust(self.aoawb, "AOA WB", self.path_fig_gust2)
        self.plot_gust(self.aoaht, "AOA HT", self.path_fig_gust3)
        self.plot_gust(self.nz, "NZ", self.path_fig_gust4)
        self.plot_gust(self.q, "Q", self.path_fig_gust5)
        self.plot_gust(self.theta, "THETA", self.path_fig_gust9)
        self.plot_gust(self.qp[0:self.nps-1], "QP", self.path_fig_gust6)
        self.plot_gust(self.wp[0:self.nps-1], "WP", self.path_fig_gust12)
        self.plot_gust(self.aoa_g[0:self.nps-1], "AG", self.path_fig_gust7)
        self.plot_gust(self.ug[0:self.nps-1], "UG", self.path_fig_gust8)
        self.plot_gust(self.aoa_kussner_wb[0:self.nps-1], "AGWBK", self.path_fig_gust10)
        self.plot_gust(self.aoa_wagner_wb[0:self.nps-1], "AGWBK", self.path_fig_gust11)
        self.plot_gust(self.w, "W", self.path_fig_gust13)            

    def get_dirs(self):
        
        main_dir = os.path.dirname(os.path.abspath(__file__))
        res_dir = os.path.join(main_dir, '../results')
        logs_dir = os.path.join(main_dir, '../logs')
        self.path_gust = os.path.join(res_dir, 'gust_results.txt')
        self.path_logs = os.path.join(logs_dir, 'log_gust.txt')
        self.path_fig_gust1 = os.path.join(res_dir, 'graph_aoa_wb_tot.png')
        self.path_fig_gust2 = os.path.join(res_dir, 'graph_aoa_wb_g.png')
        self.path_fig_gust3 = os.path.join(res_dir, 'graph_aoa_ht_g.png')
        self.path_fig_gust4 = os.path.join(res_dir, 'graph_nz.png')
        self.path_fig_gust5 = os.path.join(res_dir, 'graph_q.png')
        self.path_fig_gust6 = os.path.join(res_dir, 'graph_qp.png')
        self.path_fig_gust7 = os.path.join(res_dir, 'graph_ag.png')
        self.path_fig_gust8 = os.path.join(res_dir, 'graph_ug.png')
        self.path_fig_gust9 = os.path.join(res_dir, 'graph_theta.png')
        self.path_fig_gust10 = os.path.join(res_dir, 'graph_aoa_wb_kussner.png')
        self.path_fig_gust11 = os.path.join(res_dir, 'graph_aoa_wb_wagner.png')
        self.path_fig_gust12 = os.path.join(res_dir, 'graph_wp.png')
        self.path_fig_gust13 = os.path.join(res_dir, 'graph_w.png')

        
        with open(self.path_gust, 'w') as file:
            file.write("PERIOD: {:12.4f} | LENGTH: {:12.4f} | SIMULATED TIME: {:12.4f} | NPS: {:8d}\n". format(self.g/self.tas, self.g, self.tlimit, self.nps))
            file.write("      TIME         U            AOAWBTOT     AOAWB        AOAWB_G      AOAHT        AOAHT_G      U            W            Q            UP           WP           FX           FZ           MY           NX           NZ           \n")
            
        file.close()
        
        with open(self.path_logs, 'w') as file:
            file.close()
        
    def run_gust(self):
        
        glwb = (self.amodel.np - self.ini.ac.m.cg) * self.mac
        glht = self.ini.ac.xlah - (self.ini.ac.xlew + self.amodel.np * self.mac)
        
        dydt = np.empty(4)
        
        # Equations of Motion
        self.w[0] = self.ini.w
        self.u[0] = self.ini.u
        self.q[0] = self.ini.q
        self.theta[0] = self.ini.theta
        self.wp[0] = self.ini.wp
        self.up[0] = self.ini.up
        self.qp[0] = self.ini.qp
        self.thetap[0] = self.ini.q
        
        # Results - AoA
        self.aoawb[0] = self.ini.aoa
        self.aoawbtot[0] = self.ini.aoa
        self.aoaht[0] = self.ini.aoa_ht
        self.aoagwb[0] = 0.0
        self.aoaght[0] = 0.0        
        
        # Auxiliary - AoA
        self.aoa_kussner_wb[:] = 0.0
        self.aoa_wagner_wb[:] = 0.0        
        self.aoa_kussner_ht[:] = 0.0
        self.aoa_wagner_ht[:] = 0.0
        self.aoa_wagner_ht2[:] = 0.0 
        self.aoa_wagner_wb2[:] = 0.0              
        self.time_wb[:] = 0.0
        self.time_ht[0] = 0.0
        self.aoa_g[0] = 0.0
        self.aoa_g_ht[0] = 0.0
        self.aoa_w[0] = 0.0
        self.aoa_w_ht[0] = 0.0               

        self.nx[0] = -(((self.up[0] + self.w[0] * self.q[0]) / 9.80665) + math.sin(self.theta[0] * math.pi / 180))
        self.nz[0] = -(((self.wp[0] - self.u[0] * self.q[0]) / 9.80665) - math.cos(self.theta[0] * math.pi / 180))
        
        with open(self.path_gust, 'a') as file:
            file.write("{:12.5f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.1f} {:12.1f} {:12.1f} {:12.4f} {:12.4f}\n".format(self.t[0], self.ug[0], self.aoawbtot[0], self.aoawb[0], self.aoagwb[0], self.aoaht[0], self.aoaght[0], self.u[0], self.w[0], self.q[0], self.up[0], self.wp[0], self.ini.model.fx_b, self.ini.model.fz_b, self.ini.model.my_b, self.nx[0], self.nz[0]))
            file.close()

        for i in range(self.nps-1):
            
            for j in range(4): # Runge Kutta
            
                # Time Step
                if (j == 0):
                    trk = self.t[i]
                elif (j == 1):
                    trk = self.t[i] + 0.5 * self.dt
                elif (j == 2):
                    trk = self.t[i] + 0.5 * self.dt
                elif (j == 3):
                    trk = self.t[i] + self.dt
            
                # Time WB - Kussner
                if (trk >= self.delaywb):
                    self.time_wb[i] = trk - self.delaywb
                    ugwb = np.interp(self.time_wb[i], self.t, self.ug)
                    self.aoa_g[i] = math.atan(ugwb / self.tas)   # Gust AOA - Free stream
                else:
                    self.time_wb[i] = 0.0
                    ugwb = 0.0
                    self.aoa_g[i]  = 0.0

                # Kussner - Wing
                self.aoa_kussner_wb[i] = self.rkgwb.eval_rkstep(self.time_wb[i], self.dt, j+1, self.aoa_g[i], self.tas)
                
                # Time HT - Kussner
                if (trk >= self.delayht):
                    self.time_ht[i] = trk - self.delayht
                    self.aoa_g_ht[i] = np.interp(self.time_ht[i], self.t, self.aoa_g) - self.amodel.deda * np.interp(self.time_ht[i], self.t, self.aoa_kussner_wb) # Gust AOA - Kussner Effect
                else:
                    self.time_ht[i] = 0.0
                    self.aoa_g_ht[i] = 0.0
                
                # Kussner - HT
                self.aoa_kussner_ht[i] = self.rkght.eval_rkstep(self.time_ht[i], self.dt, j+1, self.aoa_g_ht[i], self.tas)
                
                # Time WB - Wagner
                if (i == 0):
                    self.aoa_w[i] = 0
                else:
                    self.aoa_w[i] = math.atan((self.w[i] - self.w[i-1]) / self.tas + (self.q[i] - self.q[i-1]) * glwb / self.tas)

                # Wagner - Wing
                self.aoa_wagner_wb[i] = self.rkwwb.eval_rkstep(self.time_wb[i], self.dt, j+1, self.aoa_w[i], self.tas)

                # Time HT - Wagner 
                if (i == 0):
                    self.aoa_w_ht[i] = 0
                else:
                    self.aoa_w_ht[i] = math.atan((self.w[i] - self.w[i-1]) / self.tas + (self.q[i] - self.q[i-1]) * (glwb + glht) / self.tas) - self.amodel.deda * np.interp(self.time_ht[i], self.t, self.aoa_wagner_wb)
                
                # Wagner - HT
                self.aoa_wagner_ht[i] = self.rkwht.eval_rkstep(self.time_ht[i], self.dt, j+1, self.aoa_w_ht[i], self.tas)
                
                self.aoa_wagner_wb2[i] = self.aoa_wagner_wb[i]  - self.aoa_w[i]
                self.aoa_wagner_ht2[i] = self.aoa_wagner_ht[i]  - self.aoa_w_ht[i]
                
                self.aoagwb[i] = (self.aoa_kussner_wb[i] + self.aoa_wagner_wb2[i]) * 180 / math.pi
                self.aoaght[i] = (self.aoa_kussner_ht[i] + self.aoa_wagner_ht2[i]) * 180 / math.pi             
                
                #print("[amodel] aoa: {:10.4f} aoaht: {:10.4f} q: {:10.4f}".format(self.aoagwb[i], self.aoaht[i-1], self.q[i-1]))
                self.amodel.eval_db(self.aoawb[i], self.q[i], self.aoagwb[i], self.aoaght[i], self.ini.elev, self.tas, self.ini.cond.pdyn)
                self.wp[i] = - self.q[i] * self.u[i] - 9.80665 * math.cos(self.theta[i] * math.pi / 180) + self.amodel.fz_b / self.ini.ac.m.weight
                self.up[i] = - self.w[i] * self.q[i] - 9.80665 * math.sin(self.theta[i] * math.pi / 180) + self.amodel.fx_b / self.ini.ac.m.weight 
                self.qp[i] = self.amodel.my_b / self.ini.ac.m.Iyy
                self.thetap[i] = self.q[i]
                
                dydt[0] = self.wp[i]
                dydt[1] = self.up[i]
                dydt[2] = self.qp[i]
                dydt[3] = self.thetap[i]
                
                self.rkm.rk(j+1, self.dt, self.t[i], dydt)
                
                self.w[i + 1] = self.rkm.y[j - 1, 0]
                self.u[i + 1] = self.rkm.y[j - 1, 1]
                self.q[i + 1] = self.rkm.y[j - 1, 2]
                self.theta[i + 1] = self.rkm.y[j - 1, 3]
                                    
                self.aoawb[i + 1] = math.atan(self.w[i + 1] / self.tas) * 180 / math.pi
                self.aoaht[i + 1] = self.amodel.aoaht
                self.nx[i + 1] = (((self.up[i] + self.w[i + 1] * self.q[i]) / 9.80665) + math.sin(self.theta[i + 1] * math.pi / 180))
                self.nz[i + 1] = ((self.wp[i] + self.u[i + 1] * self.q[i]) / 9.80665) + math.cos(self.theta[i + 1] * math.pi / 180)  
 

                if (j == 3):
                    self.aoawbtot[i] = self.aoawb[i] + self.aoagwb[i]
                    with open(self.path_gust, 'a') as file:
                        file.write("{:12.5f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.1f} {:12.1f} {:12.1f} {:12.4f} {:12.4f}\n".format(self.t[i], self.ug[i], self.aoawbtot[i], self.aoawb[i], self.aoagwb[i], self.aoaht[i], self.aoaght[i], self.u[i], self.w[i], self.q[i], self.up[i], self.wp[i], self.amodel.fx_b, self.amodel.fz_b, self.amodel.my_b, self.nx[i], self.nz[i]))
                        file.close()
                    
                with open(self.path_logs, 'a') as file:
                    file.write("TIME: {:12.4f} | TIME_WB: {:12.4f} | TIME_HT: {:12.4f} |  UW: {:8.4f} | RKSTEP: {:3d} | AGW: {:12.6f} |  AGH: {:12.6f} | AGWK: {:12.6f} |  AGHK: {:12.6f} | AGWW1: {:12.6f} |  AGHW1: {:12.6f} | AGWW2: {:12.6f} |  AGHW2: {:12.6f}\n". format(trk, self.time_wb[i], self.time_ht[i], ugwb, j+1, self.aoa_g[i], self.aoa_g_ht[i], self.aoa_kussner_wb[i], self.aoa_kussner_ht[i], self.aoa_wagner_wb[i], self.aoa_wagner_ht[i], self.aoa_wagner_wb2[i], self.aoa_wagner_ht2[i]))
                    
                file.close()   
                
                
    def plot_gust(self, yarray, label, file_path):
        plt.style.use('ggplot')
        
        ilen = len(yarray)
        
        # plot
        plt.plot(self.t[0:ilen], yarray)
        
        plt.xlabel('TIME')
        plt.ylabel(label)

        plt.savefig(file_path)    

        plt.clf()

        plt.close()             

                    
                
                
                    