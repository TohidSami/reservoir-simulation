# config.py
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
# ==========================================
#   USER DASHBOARD (INPUT DATA)
# ==========================================

# --- 1. SIMULATION SCENARIO ---
# Mode 1: Segregation Test (No Wells, Gravity check)
# Mode 2: Stability Test (No Wells, Uniform)
# Mode 3: Main Model (Water Injection & Production)
SCENARIO_MODE = 3

# --- 2. GRID SETTINGS ---
GRID_NX = 15
GRID_NY = 15
GRID_NZ = 3 

# --- 3. FLUID & PVT PROPERTIES ---
# Saturation Endpoints
SWC = 0.1    # Connate Water Saturation
SOR = 0.3    # Residual Oil Saturation

# Capillary Pressure: Pc = C * (1-Snw)^n
PC_COEFF = 5.0
PC_EXP   = 2.0

# Relative Permeability: Kr = Max * S^n
KRW_MAX = 0.4
KRW_EXP = 1.2
KRO_MAX = 0.35
KRO_EXP = 2.5

# Oil PVT (Linear Model)
BO_SLOPE = -6e-6
BO_INTERCEPT = 1.275
MU_OIL_LOW = 0.5   
MU_OIL_HIGH = 0.6  
P_LOW = 2500.0
P_HIGH = 3500.0

# Water PVT
BW_BASE = 1.03
CW = 3e-6
MU_WATER = 0.4

# --- 4. WELL CONTROLS (Only for Mode 3) ---
# Injector
INJ_RATE_TARGET = 1500.0  # STB/Day
INJ_BHP_MAX     = 5000.0  # psi

# Producer
PROD_RATE_TARGET = 500.0  # STB/Day (Target Oil Rate)
PROD_BHP_MIN     = 3000.0 # psi

dt=10
time_step=500

## run with this code dont change that 
if __name__ == "__main__":
    print("Loading Simulation Engine...")
    try:
        import run_model
        print(f"Scenario Mode: {SCENARIO_MODE}")
        run_model.main()
    except ImportError as e:
        print("Error: Could not find 'run_model.py'. Make sure it is in the same folder.")
        print(e)
    except Exception as e:
        print("An unexpected error occurred:")
        print(e)
        # Keep window open if run from double-click
        input("Press Enter to exit...")