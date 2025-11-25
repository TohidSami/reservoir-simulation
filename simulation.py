import numpy as np
import math

class ReservoirSolver:
    def __init__(self, grid, fluid_model, well_controls):     # Max Production Rate (Optional)

        self.grid = grid
        self.fm = fluid_model
        
        # Store Well Constraints
        self.inj_rate_target = well_controls.get('inj_rate_target',1000)
        self.inj_bhp_max = well_controls.get('inj_bhp_max',5000)
        self.prod_bhp_min = well_controls.get('prod_bhp_min',3000)
        self.prod_rate_target = well_controls.get('prod_rate_target',500)
        LA = self.grid.LA
        
        # Properties
        self.Bo = np.zeros((LA, 1))
        self.Bw = np.zeros((LA, 1))
        self.mio = np.zeros((LA, 1))
        self.miw = np.zeros((LA, 1))
        self.kro = np.zeros((LA, 1))
        self.krw = np.zeros((LA, 1))
        self.Pc  = np.zeros((LA, 1))
        self.PV  = np.zeros((LA, 1)) # Pore Volume
        
        # Derivatives (Used for Jacobian)
        self.dBo  = np.zeros((LA, 1))
        self.dBw  = np.zeros((LA, 1))
        self.dmio = np.zeros((LA, 1))
        self.dmiw = np.zeros((LA, 1)) # Likely zero if constant
        self.dKro = np.zeros((LA, 1))
        self.dKrw = np.zeros((LA, 1))
        self.dPc  = np.zeros((LA, 1))
        self.dPV  = np.zeros((LA, 1))

    def update_properties(self, x):
        """
        Calculates PVT and SCAL properties for all active blocks based on state vector x.
        Corresponds to the first ~110 lines of the old EQU.py and X2CALC logic.
        
        Args:
            x: State vector of size (2*LA, 1). Even indices=Pressure, Odd=Sw.
        """
        LA = self.grid.LA
        
        for i in range(LA):
            # Extract State Variables
            # x is a flattened vector: [P0, Sw0, P1, Sw1, ...]
            Pi = x[2*i, 0]
            Swi = x[2*i+1, 0]
            
            # Safety clamp for Saturation (to avoid negative values in formulas)
            if Swi < 0.0: Swi = 0.0
            if Swi > 1.0: Swi = 1.0

            # 1. Calculate PVT Properties (using FluidModel)
            self.Bo[i]  = self.fm.bo_cal(Pi)
            self.Bw[i]  = self.fm.bw_cal(Pi)
            self.mio[i] = self.fm.mio_cal(Pi)
            self.miw[i] = self.fm.miw_cal(Pi)
            
            # Calculate Pore Volume (assuming constant for now, but structure allows compressibility)
            # Note: In old code PV was calculated per block in GRID, but modified in X2CALC.
            # Here we take base PV from grid and can apply rock compressibility if needed.
            self.PV[i] = self.grid.PV[i][0] # / 5.615 done in FluidModel or Grid logic if needed
            
            # 2. Calculate SCAL (Saturation Dependent Properties)
            # First, normalize saturation
            Snw = self.fm.snw_cal(Swi)
            dSnw_dSw = self.fm.div_snw(Swi)
            self.krw[i] = self.fm.krw_cal(Snw) * dSnw_dSw
            self.kro[i] = self.fm.kro_cal(Snw) * dSnw_dSw
            self.Pc[i]  = self.fm.pc_cal(Snw) * dSnw_dSw
            
            # 3. Calculate Derivatives (Div Functions)
            self.dBo[i]  = self.fm.div_bo(Pi)
            self.dBw[i]  = self.fm.div_bw(Pi)
            self.dmio[i] = self.fm.div_mio(Pi)
            self.dmiw[i] = self.fm.div_miw(Pi)
            
            self.dKro[i] = self.fm.div_kro(Snw)
            self.dKrw[i] = self.fm.div_krw(Snw)
            self.dPc[i]  = self.fm.div_pc(Snw)
            
            # Rock compressibility derivative (currently zero in your model)
            self.dPV[i]  = 0

    def _calculate_density_oil(self, Bo_1, Bo_2):
        """
        Calculates average oil density between two blocks.
        Formula: rho = (rho_oil_std + Rs * Crs * rho_gas_std) / Bo
        """
        # Standard densities (lb/ft3) from original code
        ros = 37.457      # Oil density at std conditions
        rogs = 0.062428   # Gas density at std conditions
        Rs = 1000         # GOR (scf/stb)
        Crs = 1 / 5.615   # Conversion factor
        
        # Calculate density for block 1 and block 2
        rho1 = (ros + Rs * Crs * rogs) / Bo_1
        rho2 = (ros + Rs * Crs * rogs) / Bo_2
        
        # Return average
        return (rho1 + rho2) / 2

    def _calculate_density_water(self, Bw_1, Bw_2):
        """
        Calculates average water density between two blocks.
        Formula: rho = rho_water_std / Bw
        """
        rows = 62.366     # Water density at std conditions
        
        rho1 = rows / Bw_1
        rho2 = rows / Bw_2
        
        return (rho1 + rho2) / 2

    def calculate_potentials(self, x):
        """
        Calculates Phase Potentials (Pressure - Gravity Head).
        Logic iterates from top to bottom to accumulate hydrostatic pressure.
        """
        LA = self.grid.LA
        
        # Initialize Potentials with current Pressure
        # POTo = P - gamma_o * h
        # POTw = P - Pc - gamma_w * h
        self.POTo = np.zeros((LA, 1))
        self.POTw = np.zeros((LA, 1))
        
        g_factor = 1 / 144.0  # Gravity conversion factor (psi/ft approx)
        dz = self.grid.dz     # Grid block height (ft)

        for j in range(LA):
            # Current block pressure and Pc
            current_P = x[2*j, 0]
            current_Pc = self.Pc[j, 0]
            
            # Initialize base potential
            self.POTo[j, 0] = current_P
            self.POTw[j, 0] = current_P - current_Pc
            
            # Retrieve geometric layer index (K index)
            # IJ_IA maps Active Index -> [I, J, K]
            # Index 2 corresponds to 'z' or 'k'
            z_index = self.grid.IJ_IA[j][2] 
            
            # Accumulate Gravity Head from top to current block
            # This recreates the 'while z > 0' logic from original code
            
            temp_z = z_index
            current_neighbor_idx = j # Start from current block
            
            while temp_z > 0:
                # Find the neighbor ABOVE (Index 4 in NL corresponds to Top/Z-1)
                # Note: NL stores 1-based indexes, so we subtract 1.
                # If NL value is -1 (no neighbor), it handles strictly via logic check,
                # but here we assume valid grid for gravity column.
                
                neighbor_above_global = self.grid.NL[current_neighbor_idx][4]
                
                if neighbor_above_global > 0:
                    neighbor_above_idx = neighbor_above_global - 1 # 0-based
                    
                    # Fetch Bo/Bw for average density calculation
                    Bo_curr = self.Bo[current_neighbor_idx, 0]
                    Bo_above = self.Bo[neighbor_above_idx, 0]
                    
                    Bw_curr = self.Bw[current_neighbor_idx, 0]
                    Bw_above = self.Bw[neighbor_above_idx, 0]
                    
                    # Calculate Densities
                    rho_o_avg = self._calculate_density_oil(Bo_curr, Bo_above)
                    rho_w_avg = self._calculate_density_water(Bw_curr, Bw_above)
                    
                    # Add Hydrostatic Pressure (rho * g * h)
                    # Note: We ADD because we are going deeper? 
                    # Original code logic: POTo += rho * g * dz
                    # Meaning potential increases with depth relative to the surface datum.
                    self.POTo[j, 0] += rho_o_avg * g_factor * dz
                    self.POTw[j, 0] += rho_w_avg * g_factor * dz
                    
                    # Move pointer up for next iteration of while loop
                    current_neighbor_idx = neighbor_above_idx
                    temp_z -= 1
                else:
                    # Should not happen in a regular grid if z > 0, but safe break
                    break


    def _calculate_well_index(self, block_idx):
        """
        Helper to calculate Peaceman Well Index (WI) or Cf.
        Formula: WI = (c * theta * K * h) / (ln(ro/rw) + S)
        """
        # Grid properties for the block
        ix, iy, iz = self.grid.IJ_IA[block_idx]
        kx = self.grid.PERM_X[iz][iy][ix]
        ky = self.grid.PERM_Y[iz][iy][ix]
        
        dx, dy, dz = self.grid.dx, self.grid.dy, self.grid.dz
        
        # Peaceman equivalent radius (ro)
        # ro = 0.28 * [ (dy^2 * (kx/ky)^0.5) + (dx^2 * (ky/kx)^0.5) ] / [ (kx/ky)^0.5 + (ky/kx)^0.5 ]
        # Note: Your code simplified K ratio terms. I'll reproduce your exact formula logic:
        
        term1 = (kx * ky) ** 0.5
        ratio_yx = (ky / kx) ** 0.5
        ratio_xy = (kx / ky) ** 0.5
        
        numerator = (dy**2 * ratio_xy) + (dx**2 * ratio_yx)
        denominator = ratio_xy + ratio_yx
        
        ro = 0.28 * numerator / denominator
        
        # Well constants
        rw = 0.5          # Well radius (ft)
        s = 0             # Skin factor
        c = 0.008527      # Unit conversion constant
        theta = 6.2832    # 2*pi (Full circle)
        h = dz            # Net pay thickness (assuming full block height)
        
        WI = (c * theta * term1 * h) / (math.log(ro/rw) + s)
        return WI



    def calculate_residuals(self, x, dt, aco_old, acw_old):
        LA = self.grid.LA
        RO = np.zeros((LA, 1))
        RW = np.zeros((LA, 1))
        
        well_results = {
            'q_inj_w': 0, 'pbh_inj': 0,
            'q_prod_o': 0, 'q_prod_w': 0, 'pbh_prod': 0
        }

        # --- MAIN LOOP ---
        for i in range(LA):
            
            # 1. Flux Calculation (Transmissibility)
            flux_o = 0; flux_w = 0
            for k in range(6):
                neighbor_global = self.grid.NL[i][k]
                if neighbor_global > 0:
                    n = neighbor_global - 1
                    trans = self.grid.TRANS[i][k]
                    
                    # Upwind Oil
                    if self.POTo[i, 0] > self.POTo[n, 0]:
                        mob_o = self.kro[i, 0] / (self.mio[i, 0] * self.Bo[i, 0])
                    else:
                        mob_o = self.kro[n, 0] / (self.mio[n, 0] * self.Bo[n, 0])   
                    flux_o += trans * mob_o * (self.POTo[n, 0] - self.POTo[i, 0])

                    # Upwind Water
                    if self.POTw[i, 0] > self.POTw[n, 0]:
                        mob_w = self.krw[i, 0] / (self.miw[i, 0] * self.Bw[i, 0])
                    else:
                        mob_w = self.krw[n, 0] / (self.miw[n, 0] * self.Bw[n, 0])    
                    flux_w += trans * mob_w * (self.POTw[n, 0] - self.POTw[i, 0])

            # 2. Accumulation Term
            const_acc = self.grid.Vb / (5.615 * dt)
            poro = self.grid.POR
            
            aco_curr = (poro * (1 - x[2*i+1, 0])) / self.Bo[i, 0]
            acw_curr = (poro * x[2*i+1, 0]) / self.Bw[i, 0]
            
            acc_term_o = const_acc * (aco_curr - aco_old[i, 0])
            acc_term_w = const_acc * (acw_curr - acw_old[i, 0])
            
            RO[i, 0] = flux_o - acc_term_o
            RW[i, 0] = flux_w - acc_term_w

            # -----------------------------
            # 4. WELLS (Safe Logic)
            # -----------------------------
            
            # --- INJECTOR ---
            if i == 0:
                if self.inj_rate_target == 0:
                    real_rate = 0
                    real_pbh = x[2*i, 0]
                else:
                    WI_inj = self._calculate_well_index(i)
                    mob_inj = 1.0 / (self.miw[i, 0] * self.Bw[i, 0])
                    
                    calc_pbh = x[2*i, 0] + (self.inj_rate_target / (WI_inj * mob_inj))
                    
                    if calc_pbh > self.inj_bhp_max:
                        real_pbh = self.inj_bhp_max
                        real_rate = WI_inj * mob_inj * (real_pbh - x[2*i, 0])
                    else:
                        real_pbh = calc_pbh
                        real_rate = self.inj_rate_target
                
                RW[i, 0] += real_rate
                well_results['q_inj_w'] = real_rate
                well_results['pbh_inj'] = real_pbh

            # --- PRODUCER ---

            # --- PRODUCER (Well 2) ---
            elif i == LA - 1:
                target_oil = self.prod_rate_target
                
                # --- FIX: Stronger check for Zero Rate ---
                if target_oil is None or abs(target_oil) < 1e-10:
                    q_oil = 0
                    q_wat = 0
                    real_pbh = x[2*i, 0]
                else:
                    WI_prod = self._calculate_well_index(i)
                    mob_o_well = self.kro[i, 0] / (self.mio[i, 0] * self.Bo[i, 0])
                    mob_w_well = self.krw[i, 0] / (self.miw[i, 0] * self.Bw[i, 0])
                    
                    if mob_o_well <= 1e-9:

                        q_oil = 0
                        q_wat = 0
                        real_pbh = x[2*i, 0]
                    else:
                        calc_pbh = x[2*i, 0] - (target_oil / (WI_prod * mob_o_well))
                        
                        if calc_pbh < self.prod_bhp_min:
                            real_pbh = self.prod_bhp_min
                        else:
                            real_pbh = calc_pbh
                        
                        delta_p = real_pbh - x[2*i, 0]
                        q_oil = WI_prod * mob_o_well * delta_p
                        q_wat = WI_prod * mob_w_well * delta_p
                
                RO[i, 0] += q_oil
                RW[i, 0] += q_wat
                
                well_results['q_prod_o'] = q_oil
                well_results['q_prod_w'] = q_wat
                well_results['pbh_prod'] = real_pbh

        return RO, RW, well_results

    def build_jacobian(self, x, dt, well_results):
        """
        Constructs the Jacobian Matrix (Derivatives of Residuals).
        Size: (2*LA, 2*LA)
        Row 2i   : Oil Equation
        Row 2i+1 : Water Equation
        Col 2j   : Deriv wrt Pressure
        Col 2j+1 : Deriv wrt Saturation
        """
        LA = self.grid.LA
        Jac = np.zeros((2 * LA, 2 * LA))
        
        # Pre-calculate constant accumulation factor
        const_acc = self.grid.Vb / (5.615 * dt)
        poro = self.grid.POR

        # --- MAIN LOOP ---
        for i in range(LA):
            
            # --- 1. ACCUMULATION TERM DERIVATIVES (Diagonal Only) ---
            # Acc_o = const * phi * (1-Sw) / Bo
            # d(Acc_o)/dP = const * phi * (1-Sw) * (-1/Bo^2) * dBo/dP
            # d(Acc_o)/dSw = const * phi * (1/Bo) * (-1)
            
            # Oil Accumulation Derivatives
            term_o = (poro * (1 - x[2*i+1, 0])) # phi * So
            dAccO_dP = const_acc * term_o * (-1/(self.Bo[i, 0]**2)) * self.dBo[i, 0]
            dAccO_dSw = const_acc * (poro / self.Bo[i, 0]) * (-1)
            
            # Water Accumulation Derivatives
            term_w = (poro * x[2*i+1, 0]) # phi * Sw
            dAccW_dP = const_acc * term_w * (-1/(self.Bw[i, 0]**2)) * self.dBw[i, 0]
            dAccW_dSw = const_acc * (poro / self.Bw[i, 0]) * (1)
            
            # Add to Diagonal (Subtracting Accumulation in Residual eq: Flux - Acc)
            # So deriv becomes -dAcc
            Jac[2*i, 2*i]     += -dAccO_dP
            Jac[2*i, 2*i+1]   += -dAccO_dSw
            Jac[2*i+1, 2*i]   += -dAccW_dP
            Jac[2*i+1, 2*i+1] += -dAccW_dSw

            # --- 2. TRANSMISSIBILITY FLUX DERIVATIVES ---
            for k in range(6):
                neighbor_global = self.grid.NL[i][k]
                
                if neighbor_global > 0:
                    n = neighbor_global - 1 # Neighbor index (j)
                    trans = self.grid.TRANS[i][k]
                    
                    # ---------------- OIL PHASE ----------------
                    # Check Upwind Direction for Oil
                    if self.POTo[i, 0] > self.POTo[n, 0]:
                        # Flow i -> n (Upstream is i)
                        upstream = i
                        sign = 1 # Flux leaving i is positive in calculation but...
                        # Wait, Standard: Residual = Flux_in - Flux_out ...
                        # My residual calc was: Flux_net_in - Acc.
                        # Let's stick to: Flux term = Trans * Mob * (Pn - Pi)
                        # Deriv wrt Pi: Trans * [ dMob/dPi * (Pn - Pi) + Mob * (-1) ]
                        # Deriv wrt Pn: Trans * [ Mob * (1) ]
                        
                        dPot = self.POTo[n, 0] - self.POTo[i, 0] # This is negative
                        
                        # Calculate Mobility Derivatives for Upstream (i)
                        # dLam/dP
                        d_mu_B = (self.dmio[i,0]*self.Bo[i,0]) + (self.mio[i,0]*self.dBo[i,0])
                        dLamO_dP = self.kro[i,0] * (-1/(self.mio[i,0]*self.Bo[i,0])**2) * d_mu_B
                        # dLam/dS
                        dLamO_dS = (1/(self.mio[i,0]*self.Bo[i,0])) * self.dKro[i,0]
                        
                        mob_o = self.kro[i,0] / (self.mio[i,0]*self.Bo[i,0])
                        
                        # Diagonal (wrt i variables)
                        Jac[2*i, 2*i]   += trans * (dLamO_dP * dPot + mob_o * (-1))
                        Jac[2*i, 2*i+1] += trans * (dLamO_dS * dPot)
                        
                        # Off-Diagonal (wrt n variables)
                        # Mobility is from i, so dMob/dPn = 0
                        Jac[2*i, 2*n]   += trans * (mob_o * 1)
                        # Jac[2*i, 2*n+1] is 0
                        
                    else:
                        # Flow n -> i (Upstream is n)
                        dPot = self.POTo[n, 0] - self.POTo[i, 0] # This is positive
                        
                        # Mobility Derivatives for Upstream (n)
                        d_mu_B = (self.dmio[n,0]*self.Bo[n,0]) + (self.mio[n,0]*self.dBo[n,0])
                        dLamO_dP = self.kro[n,0] * (-1/(self.mio[n,0]*self.Bo[n,0])**2) * d_mu_B
                        dLamO_dS = (1/(self.mio[n,0]*self.Bo[n,0])) * self.dKro[n,0]
                        
                        mob_o = self.kro[n,0] / (self.mio[n,0]*self.Bo[n,0])
                        
                        # Diagonal (wrt i variables)
                        Jac[2*i, 2*i]   += trans * (mob_o * (-1))
                        
                        # Off-Diagonal (wrt n variables)
                        Jac[2*i, 2*n]   += trans * (dLamO_dP * dPot + mob_o * 1)
                        Jac[2*i, 2*n+1] += trans * (dLamO_dS * dPot)

                    # ---------------- WATER PHASE ----------------
                    # Same logic, just careful with Pc in potential if needed, 
                    # but here derivatives are wrt P and Sw directly.
                    # Note: PotW = P - Pc. 
                    # d(PotW_i)/dP_i = 1.  d(PotW_i)/dS_i = -dPc/dSw.
                    
                    # Check Upwind Direction for Water
                    if self.POTw[i, 0] > self.POTw[n, 0]:
                        # Flow i -> n
                        upstream = i
                        dPot = self.POTw[n, 0] - self.POTw[i, 0]
                        
                        # Mobility Derivs (i)
                        d_mu_B = (self.dmiw[i,0]*self.Bw[i,0]) + (self.miw[i,0]*self.dBw[i,0])
                        dLamW_dP = self.krw[i,0] * (-1/(self.miw[i,0]*self.Bw[i,0])**2) * d_mu_B
                        dLamW_dS = (1/(self.miw[i,0]*self.Bw[i,0])) * self.dKrw[i,0]
                        mob_w = self.krw[i,0] / (self.miw[i,0]*self.Bw[i,0])
                        
                        # d(dPot)/dP_i = -1
                        # d(dPot)/dS_i = -(-dPc_i) = +dPc_i (Since PotW = P - Pc)
                        dPot_dPi = -1
                        dPot_dSi = self.dPc[i, 0] 
                        
                        # Diagonal (wrt i)
                        # Flux = Trans * Mob_i * (Pot_n - Pot_i)
                        # dFlux/dPi = Trans * [ dMob/dPi * dPot + Mob * (-1) ]
                        Jac[2*i+1, 2*i]   += trans * (dLamW_dP * dPot + mob_w * dPot_dPi)
                        
                        # dFlux/dSi = Trans * [ dMob/dSi * dPot + Mob * (dPot/dSi) ]
                        Jac[2*i+1, 2*i+1] += trans * (dLamW_dS * dPot + mob_w * dPot_dSi)
                        
                        # Off-Diagonal (wrt n)
                        # dFlux/dPn = Trans * Mob * (1)
                        Jac[2*i+1, 2*n]   += trans * (mob_w * 1)
                        # dFlux/dSn = Trans * Mob * (-dPot_n/dSn) = Trans * Mob * (- (-dPc_n)) = Trans*Mob*dPc_n
                        # Wait, Pot_n = Pn - Pc_n. dPot_n/dSn = -dPc/dSn
                        dPot_dSn = -self.dPc[n, 0]
                        Jac[2*i+1, 2*n+1] += trans * (mob_w * dPot_dSn)

                    else:
                        # Flow n -> i
                        dPot = self.POTw[n, 0] - self.POTw[i, 0]
                        
                        # Mobility Derivs (n)
                        d_mu_B = (self.dmiw[n,0]*self.Bw[n,0]) + (self.miw[n,0]*self.dBw[n,0])
                        dLamW_dP = self.krw[n,0] * (-1/(self.miw[n,0]*self.Bw[n,0])**2) * d_mu_B
                        dLamW_dS = (1/(self.miw[n,0]*self.Bw[n,0])) * self.dKrw[n,0]
                        mob_w = self.krw[n,0] / (self.miw[n,0]*self.Bw[n,0])
                        
                        # d(dPot)/dP_i = -1
                        # d(dPot)/dS_i = -(-dPc_i) = dPc_i
                        dPot_dPi = -1
                        dPot_dSi = self.dPc[i, 0]
                        
                        # Diagonal (wrt i) - Only affects Potential difference, not mobility
                        Jac[2*i+1, 2*i]   += trans * (mob_w * dPot_dPi)
                        Jac[2*i+1, 2*i+1] += trans * (mob_w * dPot_dSi)
                        
                        # Off-Diagonal (wrt n) - Affects both Mob and Potential
                        # d(dPot)/dPn = 1
                        # d(dPot)/dSn = -dPc_n
                        dPot_dPn = 1
                        dPot_dSn = -self.dPc[n, 0]
                        
                        Jac[2*i+1, 2*n]   += trans * (dLamW_dP * dPot + mob_w * dPot_dPn)
                        Jac[2*i+1, 2*n+1] += trans * (dLamW_dS * dPot + mob_w * dPot_dSn)

            # --- 3. WELL DERIVATIVES ---
            # Injector at i=0
            if i == 0:
                # If Rate Constrained (Fixed Rate): Q is constant, Derivs are 0.
                # If BHP Constrained: Q = WI * Mob * (BHP - P)
                # We need to know which mode is active.
                # We can check Pbh from well_results to infer mode or recalculate logic.
                
                target_rate = self.inj_rate_target
                WI_inj = self._calculate_well_index(i)
                mob_inj = 1.0 / (self.miw[i, 0] * self.Bw[i, 0])
                calc_pbh = x[2*i, 0] + (target_rate / (WI_inj * mob_inj))
                
                if calc_pbh > self.inj_bhp_max:
                    # BHP Control Active (Pressure Constrained)
                    # Q = WI * Mob * (BHP_limit - P)
                    # This adds to RW (Water Eq).
                    # dQ/dP = WI * [ dMob/dP * (BHP-P) + Mob * (-1) ]
                    
                    real_pbh = self.inj_bhp_max
                    d_mu_B = (self.dmiw[i,0]*self.Bw[i,0]) + (self.miw[i,0]*self.dBw[i,0])
                    dMob_dP = -1/(self.miw[i,0]*self.Bw[i,0])**2 * d_mu_B
                    
                    dQ_dP = WI_inj * (dMob_dP * (real_pbh - x[2*i, 0]) + mob_inj * (-1))
                    
                    # Add to Jacobian (Row 2i+1, Col 2i)
                    Jac[2*i+1, 2*i] += dQ_dP
                    
                    # dQ/dS = 0 (Mobility of water inj doesn't depend on S here, pure water)

            # Producer at i=LA-1
# --- PRODUCER (Well 2) ---
            elif i == LA - 1:
                target_oil = self.prod_rate_target
                
                if target_oil is None or abs(target_oil) < 1e-10:
                    pass 
                else:
                    WI_prod = self._calculate_well_index(i)
                    mob_o_well = self.kro[i, 0] / (self.mio[i, 0] * self.Bo[i, 0])
                    mob_w_well = self.krw[i, 0] / (self.miw[i, 0] * self.Bw[i, 0])
                    
                    if mob_o_well <= 1e-9:
                        pass 
                    else:
                      
                        calc_pbh = x[2*i, 0] - (target_oil / (WI_prod * mob_o_well))
                        
                        is_bhp_control = False
                        if calc_pbh < self.prod_bhp_min:
                            is_bhp_control = True
                            real_pbh = self.prod_bhp_min
                        else:
                        
                            is_bhp_control = False
                            real_pbh = calc_pbh
                        
                        if is_bhp_control:
                            # اگر روی فشار ثابت قفل شده باشیم، دبی با فشار تغییر می‌کند
                            # پس باید مشتقات را حساب کنیم.
                            
                            # Oil Derivatives (مشتقات دبی نفت نسبت به P و S)
                            # dQo/dP = WI * [ dMobO/dP * (Pbh - P) + MobO * (-1) ]
                            
                            d_mu_B_o = (self.dmio[i,0]*self.Bo[i,0]) + (self.mio[i,0]*self.dBo[i,0])
                            dMobO_dP = self.kro[i,0] * (-1/(self.mio[i,0]*self.Bo[i,0])**2) * d_mu_B_o
                            dMobO_dS = (1/(self.mio[i,0]*self.Bo[i,0])) * self.dKro[i,0]
                            
                            dQo_dP = WI_prod * (dMobO_dP * (real_pbh - x[2*i, 0]) + mob_o_well * (-1))
                            dQo_dS = WI_prod * (dMobO_dS * (real_pbh - x[2*i, 0]))
                            
                            Jac[2*i, 2*i]   += dQo_dP
                            Jac[2*i, 2*i+1] += dQo_dS
                            
                            # Water Derivatives (مشتقات دبی آب نسبت به P و S)
                            # dQw/dP = WI * [ dMobW/dP * (Pbh - P) + MobW * (-1) ]
                            
                            d_mu_B_w = (self.dmiw[i,0]*self.Bw[i,0]) + (self.miw[i,0]*self.dBw[i,0])
                            dMobW_dP = self.krw[i,0] * (-1/(self.miw[i,0]*self.Bw[i,0])**2) * d_mu_B_w
                            dMobW_dS = (1/(self.miw[i,0]*self.Bw[i,0])) * self.dKrw[i,0]
                            
                            dQw_dP = WI_prod * (dMobW_dP * (real_pbh - x[2*i, 0]) + mob_w_well * (-1))
                            dQw_dS = WI_prod * (dMobW_dS * (real_pbh - x[2*i, 0]))
                            
                            Jac[2*i+1, 2*i]   += dQw_dP
                            Jac[2*i+1, 2*i+1] += dQw_dS

        return Jac