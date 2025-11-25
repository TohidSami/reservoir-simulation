# accumulation.py
import numpy as np

def calculate_accumulation_at_t_n(grid, sim, x_old):
    """
    Calculates Accumulation terms at time step n (Old Time Level).
    Logic ported from old BLOCK_SAVE.py.
    """
    # 1. Update simulator properties using the OLD state variables
    # This ensures Bo and Bw are calculated at Pressure_old
    sim.update_properties(x_old)
    
    poro = grid.POR
    LA = grid.LA
    
    aco_old = np.zeros((LA, 1))
    acw_old = np.zeros((LA, 1))
    
    for i in range(LA):
        # x indices: Even=Pressure, Odd=Saturation (Sw)
        Sw_old = x_old[2*i+1, 0]
        So_old = 1 - Sw_old
        
        # Calculate (phi * S) / B
        # Note: sim.Bo is now updated to reflect Bo(P_old)
        aco_old[i, 0] = (poro * So_old) / sim.Bo[i, 0]
        acw_old[i, 0] = (poro * Sw_old) / sim.Bw[i, 0]
        
    return aco_old, acw_old