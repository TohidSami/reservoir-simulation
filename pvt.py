# pvt.py
import numpy as np

class FluidModel:
    def __init__(self, input):
                 # 1. Saturation Endpoints
                 self.swc=input.get('swc',0.1)    # Connate Water Saturation
                 self.sor=input.get('sor',0.3)    # Residual Oil Saturation
                 
                 # 2. Capillary Pressure Parameters: Pc = C * (1-Snw)^n
                 self.pc_coeff=input.get('pc_coeff',5.0) # C
                 self.pc_exp=input.get('pc_exp',2.0)   # n
                 
                 # 3. Relative Permeability Parameters: Kr = K_max * S^n
                 self.krw_max=input.get('krw_max',0.4)  # C for water
                 self.krw_exp=input.get('krw_exp',1.2)  # n for water
                 self.kro_max=input.get('kro_max',0.35) # C for oil
                 self.kro_exp=input.get('kro_exp',2.5)  # n for oil
                 
                 # 4. Oil FVF (Linear): Bo = slope * P + intercept
                 self.bo_slope=input.get('bo_slope',-6e-6)
                 self.bo_intercept=input.get('bo_intercept',1.275)
                 
                 # 5. Water FVF (Quadratic): Bw = Bwi * (1 + Cw*dP + 0.5*Cw^2*dP^2)
                 self.bw_base=input.get('bw_base',1.03) # Bwi
                 self.cw=input.get('cw',3e-6)      # Compressibility
                 self.p_ref=input.get('p_ref',4000.0) # Reference Pressure
                 
                 # 6. Water Viscosity (Constant)
                 self.muw=input.get('muw',0.4)
                 
                 # 7. Oil Viscosity (3 Regions)
                 # Region 1: P <= p_low -> mu = mu_oil_low
                 # Region 3: P >= p_high -> mu = mu_oil_high
                 # Region 2: Linear interpolation in between
                 self.p_low=input.get('p_low',2500.0)
                 self.p_high=input.get('p_high',3500.0)
                 self.mu_oil_low=input.get('mu_oil_low',0.5) 
                 self.mu_oil_high=input.get('mu_oil_high',0.6) 
                 
        
        
        # self.swc = swc
        # self.sor = sor
        
        # self.pc_coeff = pc_coeff
        # self.pc_exp = pc_exp
        
        # self.krw_max = krw_max
        # self.krw_exp = krw_exp
        # self.kro_max = kro_max
        # self.kro_exp = kro_exp
        
        # self.bo_slope = bo_slope
        # self.bo_intercept = bo_intercept
        
        # self.bw_base = bw_base
        # self.cw = cw
        # self.p_ref = p_ref
        
        # self.muw = muw
        
        # # Oil Viscosity Parameters
        # self.p_low = p_low
        # self.p_high = p_high
        # self.mu_oil_low = mu_oil_low
        # self.mu_oil_high = mu_oil_high
        
     
        # Slope m = (y2 - y1) / (x2 - x1)
                 if self.p_high != self.p_low:
                    self.mu_slope = (self.mu_oil_high - self.mu_oil_low) / (self.p_high - self.p_low)
                 else:
                    self.mu_slope = 0

    # ==========================
    # Main Calculation Functions
    # ==========================

    def snw_cal(self, swi):
        """Normalized Water Saturation"""
       
        denominator = 1.0 - self.swc - self.sor
        if denominator == 0: return 0
        
        snw = (swi - self.swc) / denominator
        
        if snw < 0: snw = 0
        if snw > 1: snw = 1
        return snw

    def pc_cal(self, snwi):
        """Capillary Pressure"""
        return self.pc_coeff * (1 - snwi) ** self.pc_exp

    def krw_cal(self, snwi):
        """Relative Permeability Water"""
        return self.krw_max * (snwi) ** self.krw_exp

    def kro_cal(self, snwi):
        """Relative Permeability Oil"""
        return self.kro_max * (1 - snwi) ** self.kro_exp

    def bo_cal(self, pi):
        """Oil FVF (Linear)"""
        return self.bo_slope * pi + self.bo_intercept

    def bw_cal(self, pi):
        """Water FVF (Quadratic)"""
        dp = pi - self.p_ref
        term1 = self.cw * dp
        term2 = 0.5 * (self.cw ** 2) * (dp ** 2)
        return self.bw_base * (1 + term1 + term2)

    def mio_cal(self, pi):
        """Oil Viscosity (3 Regions)"""
        if pi <= self.p_low:
            return self.mu_oil_low
        elif pi >= self.p_high:
            return self.mu_oil_high
        else:
            
            # y = y1 + m * (x - x1)
            return self.mu_oil_low + self.mu_slope * (pi - self.p_low)

    def miw_cal(self, pi):
        """Water Viscosity (Constant)"""
        return self.muw

    # ==========================
    # Derivative Functions (Div)
    # ==========================

    def div_snw(self, swi):
        """d(Snw) / d(Sw)"""
      
        return 1.0 / (1.0 - self.swc - self.sor)

    def div_pc(self, snwi):
        """d(Pc) / d(Snw)"""
        # d/dx [ C * (1-x)^n ] = C * n * (1-x)^(n-1) * (-1)
        term = (1 - snwi)
        if term < 0: term = 0 
        
        if self.pc_exp == 2:
             return self.pc_coeff * self.pc_exp * term * (-1)
        else:
             
             if term == 0: return 0
             return -1 * self.pc_coeff * self.pc_exp * (term ** (self.pc_exp - 1))

    def div_krw(self, snwi):
        """d(Krw) / d(Snw)"""
        # d/dx [ C * x^n ] = C * n * x^(n-1)
        if snwi <= 0: return 0
        return self.krw_max * self.krw_exp * (snwi ** (self.krw_exp - 1))

    def div_kro(self, snwi):
        """d(Kro) / d(Snw)"""
        # d/dx [ C * (1-x)^n ] = -C * n * (1-x)^(n-1)
        term = 1 - snwi
        if term <= 0: return 0
        return -1 * self.kro_max * self.kro_exp * (term ** (self.kro_exp - 1))

    def div_bo(self, pi):
        """d(Bo) / d(P)"""
      
        return self.bo_slope

    def div_bw(self, pi):
        """d(Bw) / d(P)"""
       
        # Bw = Bwi * (1 + Cw*dP + 0.5*Cw^2*dP^2)
        # dBw/dP = Bwi * (Cw + Cw^2 * dP)
        dp = pi - self.p_ref
        return self.bw_base * (self.cw + (self.cw ** 2) * dp)

    def div_mio(self, pi):
        """d(Muo) / d(P)"""
        if pi <= self.p_low:
            return 0
        elif pi >= self.p_high:
            return 0
        else:
            return self.mu_slope 

    def div_miw(self, pi):
        """d(Muw) / d(P)"""
        return 0 