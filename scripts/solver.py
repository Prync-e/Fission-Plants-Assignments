from pyfluids import Fluid, FluidsList
from scipy import constants as cs
from data.pipe_handler import get_diameter as d_list
from data.pipe_handler import ROWS
import data.assignment_data as dh
from scripts.correlations import properties, Reynolds, friction_smooth

# Solution of point a
def solve_point_a(liquid_phase: Fluid, steam_phase: Fluid):  
    # Saturated water properties (density, enthalpy)
    rho_liquid, h_liquid, _ = properties(liquid_phase)
    
    # Saturated steam properties (density, enthalpy)
    rho_steam, h_steam, _ = properties(steam_phase)
    
    # Mass flow rate (knowing P and Delta h)
    mass_rate  = dh.Pth/(h_steam-h_liquid)
    
    # Compute H for each diameter of the pipe
    for i in range(ROWS):
        Pipe_name, D_pipe = d_list(i)
        
        # Speed of fluids based on mass flow rate
        v_liquid = mass_rate/rho_liquid*4/cs.pi/D_pipe**2
        v_steam = mass_rate/rho_steam*4/cs.pi/D_pipe**2
    
        # pressure drop due to localized losses, The velocities were assigned by prof's slide
        Dp_h = dh.k_h*rho_steam*v_steam**2/2
        Dp_c = dh.k_c*rho_liquid*v_liquid**2/2
    
        # Re
        Re_liquid = Reynolds(liquid_phase, v_liquid, D_pipe)
        Re_steam = Reynolds(steam_phase, v_liquid, D_pipe)
    
        # friction factor
        f_liquid = friction_smooth(Re_liquid)
        f_steam = friction_smooth(Re_steam)
        
        # distributed losses
        Dp_steam=dh.L/D_pipe/2*f_steam*rho_steam*v_steam**2
        Dp_liquid= dh.L/D_pipe/2*f_liquid*rho_liquid*v_liquid**2
                        
        h = (Dp_c+Dp_h+Dp_liquid+Dp_steam)/(cs.g*(rho_liquid-rho_steam))
        print(f"{Pipe_name}: h = {h} m")


def solver():
    water = Fluid(FluidsList.Water)

    # Liquid objects from pyFluids
    liquid_phase = water.bubble_point_at_pressure(dh.p)
    steam_phase = water.dew_point_at_pressure(dh.p)

    solve_point_a(liquid_phase, steam_phase)
  