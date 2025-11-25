import pytest
from pvt import FluidModel
def test_oil_viscosity_low_pressure():
    input={
        'mu_oil_low': 0.4,
        'p_low': 2000
    }
    fm=FluidModel(input)
    muo=fm.mio_cal(1500)
    assert muo==0.4

def test_water_fvf_calculation():
    # input={
    #     'cw':3e-6
    # }
    # fm=FluidModel(input)
    # bw=fm.bw_cal(5000)
    # return print(bw)

    inputs = {'bw_base': 1.0, 'cw': 0.0, 'p_ref': 4000}
    fm = FluidModel(inputs)
    
    assert fm.bw_cal(5000) == 1.0