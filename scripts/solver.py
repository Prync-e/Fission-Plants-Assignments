from pyfluids import Fluid, FluidsList, Input
from scipy import constants as cs
import data.pipe_handler as pp
import data.assignment_data as dh
import scripts.correlations as rr
import math as m

def solver():
    # Problem setup
    water = Fluid(FluidsList.Water)
    
    # HX2 fixed parameter
    q2_HX2 = dh.Pth/dh.Atransfer_HX2
    Delta_Tsat = rr.inverse_McAdams(q2_HX2)
    h_out_HX2 = q2_HX2/Delta_Tsat
    Din_HX2 = dh.Dout_HX2 - 2*dh.thick_HX2
    
    # Iterative quantities
    err = 1
    toll = 1.e-5
    cc = 0
    COUNT = 1e4
    
    # Initial guesses
    mass_flow = 100     # kg/s
    T_hot = 160         # °C
    
    # Temperature independent parameters
    Delta_h = dh.Pth/mass_flow
    mass_flow_pipe_HX2 = mass_flow/dh.Ntubes_HX2

    while err > toll and cc < COUNT:
        print(f'{cc} iter: T hot = {T_hot:.4f}, relative err: {err:.3e}')
        
        # Fluid states definition
        w_hot = water.with_state(Input.temperature(T_hot), Input.pressure(dh.p_ISC))
        h_hot = w_hot.enthalpy
        h_cold = h_hot - Delta_h
        w_cold = water.with_state(Input.pressure(dh.p_ISC), Input.enthalpy(h_cold))
        T_cold = w_cold.temperature
        T_avg = (T_hot + T_cold)/2
        w_avg = water.with_state(Input.temperature(T_avg), Input.pressure(dh.p_ISC))
        
        # Fluid properties in HX2 (T_avg) 
        Nu = rr.Nusselt(w_avg, mass_flow_pipe_HX2, Din_HX2, dh.epsilon_rel_HX2)
        k_w = w_avg.conductivity
        cp_w = w_avg.specific_heat
        
        h_in_HX2 = Nu * k_w / Din_HX2
        
        # Solving HX2
        inverse_Uw_HX2 = dh.Dout_HX2/Din_HX2/h_in_HX2 + dh.Dout_HX2*m.log(dh.Dout_HX2/Din_HX2)/2/dh.k_HX2 + 1/h_out_HX2
        Uw_HX2 = pow(inverse_Uw_HX2 ,-1)
        NTUw_HX2 = Uw_HX2*dh.Atransfer_HX2/mass_flow/cp_w
           
        DeltaT_ML_HX2 = q2_HX2 / Uw_HX2
        
        # Fining the new T_hot
        T_hot_new = 100 + DeltaT_ML_HX2*NTUw_HX2/(1-m.exp(-NTUw_HX2))
        
        # Updating loop parameters
        cc += 1
        err = abs(T_hot-T_hot_new)/T_hot
        
        T_hot = T_hot_new
    
    print(f'{cc} iter: T hot = {T_hot:.4f}, relative err: {err:.3e}')
        