# run_model.py
import numpy as np
# import config as cfg
from mygrid import Grid
from pvt import FluidModel
from simulation import ReservoirSolver
from accumulation import calculate_accumulation_at_t_n
from exporter import write_results_to_txt
# from utils import run_animation

def get_initial_conditions(mode, grid,well_p):
    """
    Set initial P and Sw based on Scenario Mode.
    Uses Geometric Coordinates (ix, iy, iz) for robust splitting.
    """
    x = np.zeros((2 * grid.LA, 1))
    
    if mode == 1: # Segregation
        print(f"--- Mode 1: Segregation Test ---")
        print(f"Splitting Grid along Y-Axis (Midpoint: {grid.NY // 2})")
        well_p['prod_rate_target']=0
        well_p['inj_rate_target']=0
        mid_y = grid.NY // 2
        split_index = int(grid.LA * 0.98)
        for idx in range(grid.LA):
            ix, iy, iz = grid.IJ_IA[idx]
            
            x[2*idx, 0] = 4000.0 
            
            if iy < mid_y:
                x[2*idx+1, 0] = 0.12 # Oil Zone (Low Sw)
            else:
                x[2*idx+1, 0] = 0.68 # Water Zone (High Sw)
                
    elif mode == 2: # Stability
        print(f"--- Mode 2: Stability Test ---")
        well_p['prod_rate_target']=0
        well_p['inj_rate_target']=0
        for i in range(grid.LA):
            x[2*i, 0] = 4000.0
            x[2*i+1, 0] = 0.15 # Uniform
            
    elif mode == 3: # Main Model
        print(f"--- Mode 3: Main Reservoir Model ---")
        for i in range(grid.LA):
            x[2*i, 0] = 4000.0
            x[2*i+1, 0] = 0.12 # Initial Sw
            
    return x, well_p
def map_db_to_engine(db_case):
    case = {
        'grid': {
            'nx': getattr(db_case, 'nx', 15),
            'ny': getattr(db_case, 'ny', 15),
            'nz': getattr(db_case, 'nz', 3)
        },
        'pvt': {
            'swc': getattr(db_case, 'swc', 0.1),
            'sor': getattr(db_case, 'sor', 0.3),
            
            'pc_coeff': getattr(db_case, 'pc_coeff', 5.0),
            'pc_exp': getattr(db_case, 'pc_exp', 2.0),
            
            'krw_max': getattr(db_case, 'krw_max', 0.4),
            'krw_exp': getattr(db_case, 'krw_exp', 1.2),
            'kro_max': getattr(db_case, 'kro_max', 0.35),
            'kro_exp': getattr(db_case, 'kro_exp', 2.5),
            
            'bo_slope': getattr(db_case, 'bo_slope', -6e-6),
            'bo_intercept': getattr(db_case, 'bo_intercept', 1.275),
            'mu_oil_low': getattr(db_case, 'mu_oil_low', 0.5),
            'mu_oil_high': getattr(db_case, 'mu_oil_high', 0.6),
            'p_low': getattr(db_case, 'p_low', 2500.0),
            'p_high': getattr(db_case, 'p_high', 3500.0),
            
            'bw_base': getattr(db_case, 'bw_base', 1.03),
            'cw': getattr(db_case, 'cw', 3e-6),
            'muw': getattr(db_case, 'muw', 0.4),
            
            'p_ref': getattr(db_case, 'p_ref', 4000.0)
        },
        'wells': {
            'inj_rate_target': getattr(db_case, 'inj_rate_target', 1500.0),
            'inj_bhp_max': getattr(db_case, 'inj_bhp_max', 5000.0),
            'prod_rate_target': getattr(db_case, 'prod_rate_target', 300.0),
            'prod_bhp_min': getattr(db_case, 'prod_bhp_min', 3000.0)
        },
        'time': {
            'dt': getattr(db_case, 'dt', 10.0),
            'n_steps': getattr(db_case, 'n_steps', 50),
            'SCENARIO_MODE': getattr(db_case, 'scenario_mode', 3)
        }
    }
    return case
def run_simulation(case):
    grid_p=case['grid']
    pvt_p=case['pvt']
    well_p=case['wells']
    time_p=case['time']
    my_grid = Grid(nx=grid_p['nx'], ny=grid_p['ny'], nz=grid_p['nz'])

    my_fm = FluidModel(pvt_p)
    
    # Configure Wells
    mode=time_p['SCENARIO_MODE']
    dt=time_p['dt']
    n_steps=time_p['n_steps']
  
    
    x , well_p= get_initial_conditions(mode, my_grid,well_p)
    history_sw = []
    history_p=[]
    times=[]
    current_time=0.0
    well_properties={
              'q_inj_w': [], 'pbh_inj': [],
            'q_prod_o': [], 'q_prod_w': [], 'pbh_prod': []
    }
    sim = ReservoirSolver(grid=my_grid, fluid_model=my_fm,well_controls=well_p )  
    print(f"Running Simulation... (Steps: {n_steps}, dt: {dt})")
    
    for step in range(n_steps):
        aco_old, acw_old = calculate_accumulation_at_t_n(my_grid, sim, x)
        x_new = x.copy()
        converged = False
        current_time +=dt
        times.append(current_time)
        for iter_idx in range(20):
            sim.update_properties(x_new)
            sim.calculate_potentials(x_new)
            RO, RW, well_res = sim.calculate_residuals(x_new, dt, aco_old, acw_old)
                        
            F = np.zeros((2 * my_grid.LA, 1))
            for i in range(my_grid.LA):
                F[2*i, 0] = RO[i, 0]
                F[2*i+1, 0] = RW[i, 0]
            
            if np.sum(np.abs(F)) < 1e-4:
                converged = True
                break
            


            Jac = sim.build_jacobian(x_new, dt, well_res)
            try:
                delta_x = np.linalg.solve(Jac, -F)
                # x_new = x_new + delta_x
            except np.linalg.LinAlgError:
                print("Solver Failed: Singular Matrix")
                break

            dp = delta_x[0::2]  
            ds = delta_x[1::2]  
            dp_clamped = np.clip(dp, -500.0, 500.0)
            ds_clamped = np.clip(ds, -0.2, 0.2)
            delta_x[0::2] = dp_clamped
            delta_x[1::2] = ds_clamped   
            x_new = x_new + delta_x         
        x = x_new.copy()
        
        sw_grid3 = np.zeros((my_grid.NZ,my_grid.NY, my_grid.NX))
        sw_grid3.fill(np.nan)
        p_grid3 = np.zeros((my_grid.NZ,my_grid.NY, my_grid.NX))
        p_grid3.fill(np.nan)
        
        for idx in range(my_grid.LA):
            ix, iy, iz = my_grid.IJ_IA[idx]
            sw_grid3[iz,iy, ix] = x[2*idx+1, 0]
            p_grid3[iz,iy,ix]=x[2*idx,0]
                
        history_sw.append(sw_grid3.copy())
        history_p.append(p_grid3.copy())
        well_properties['pbh_inj'].append(well_res['pbh_inj'])
        well_properties['pbh_prod'].append(well_res['pbh_prod'])
        well_properties['q_inj_w'].append(well_res['q_inj_w'])
        well_properties['q_prod_o'].append(well_res['q_prod_o'])
        well_properties['q_prod_w'].append(well_res['q_prod_w'])
        # if step % 10 == 0:
        print(f"Step {step} | Inj: {well_res['q_inj_w']:.1f} | ProdO: {well_res['q_prod_o']:.1f} | ProdW: {well_res['q_prod_w']:.1f} | bhpP: {well_res['pbh_prod']:.1f} | phbI: {well_res['pbh_inj']:.1f}")
        ###### for WC #########
        if np.abs(well_res['q_prod_w'])> np.abs(well_res['q_prod_o']) or well_res['q_prod_o']==0:
            print(" wc is more than 50 percent or qo =0")
            write_results_to_txt(history_sw, history_p, my_grid, times)
            return history_sw, history_p, times,well_properties
    write_results_to_txt(history_sw, history_p, my_grid, times)    
    print("Simulation Finished.")
    return history_sw, history_p, times ,well_properties


def main():
    mock_db_data = {
        'grid': {'nx': 15, 'ny': 15, 'nz': 3},
        'pvt': {
            'swc': 0.1, 'sor': 0.3,
            'pc_coeff': 5.0, 'pc_exp': 2.0,
            'krw_max': 0.4, 'krw_exp': 1.2,
            'kro_max': 0.35, 'kro_exp': 2.5,
            'bo_slope': -6e-6, 'bo_intercept': 1.275,
            'bw_base': 1.03, 'cw': 3e-6,
            'mu_oil_low': 0.5, 'mu_oil_high': 0.6,
            'p_low': 2500, 'p_high': 3500
        },
        'wells': {
            'inj_rate_target': 1500.0,
            'inj_bhp_max': 5000.0,
            'prod_rate_target': 500.0,
            'prod_bhp_min': 3000.0
        },
        'time': {'dt': 10, 'n_steps': 50, 'SCENARIO_MODE':3}
    }

    # Run without config file!
    from utils import run_animation
    
    results = run_simulation(mock_db_data)
    
    temp_grid = Grid(15, 15, 3) 
    run_animation(results, temp_grid, "Refactored Run")

if __name__ == "__main__":
    print("--- Running in Test Mode with Mock Data ---")

    class MockCase:
        # --- Grid ---
        nx = 15
        ny = 15
        nz = 3
        
        # --- Saturation & Rock ---
        swc = 0.1
        sor = 0.3
        
        # --- PC & Rel Perm ---
        pc_coeff = 5.0
        pc_exp = 2.0
        krw_max = 0.4
        krw_exp = 1.2
        kro_max = 0.35
        kro_exp = 2.5
        
        # --- Oil PVT ---
        bo_slope = -6e-6
        bo_intercept = 1.275
        mu_oil_low = 0.5
        mu_oil_high = 0.6
        p_low = 2500.0
        p_high = 3500.0
        
        # --- Water PVT ---
        bw_base = 1.03
        cw = 3e-6
        muw = 0.4
        p_ref = 4000.0
        p_init = 4000.0  
        
        # --- Wells ---
        inj_rate_target = 1500.0
        inj_bhp_max = 6000.0
        prod_rate_target = 500.0
        prod_bhp_min = 2500.0
        
        # --- Time ---
        dt = 10.0
        n_steps = 5
        scenario_mode = 3
    mock_case = MockCase()
    
    case = map_db_to_engine(mock_case)
    
    run_simulation(case)