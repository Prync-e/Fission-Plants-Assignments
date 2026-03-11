from pyfluids import Fluid, FluidsList
from scipy import constants as cs
import data.pipe_handler as pp
import data.assignment_data as dh
import scripts.correlations as rr

# Total localized losses
def total_loc(rho_liquid: float, rho_steam: float, v_liquid: float, v_steam: float) -> float:
    Dp_c = rr.localized_loss(dh.k_c, rho_liquid, v_liquid)              # Cooler (liquid)
    Dp_h = rr.localized_loss(dh.k_h, rho_steam, v_steam)                # Heater (steam)
    Dp_e_liquid = 2*rr.localized_loss(dh.k_e, rho_liquid, v_liquid)     # Pipe elbows (liquid)
    Dp_e_steam = 2*rr.localized_loss(dh.k_e, rho_steam, v_steam)        # Pipe elbows (steam)
    
    return Dp_c + Dp_h + Dp_e_liquid + Dp_e_steam

# Solution of point a
def solve_point_a(liquid_phase: Fluid, steam_phase: Fluid) -> dict:  
    # Saturated water properties (density, enthalpy)
    rho_liquid, h_liquid, _ = rr.properties(liquid_phase)
    
    # Saturated steam properties (density, enthalpy)
    rho_steam, h_steam, _ = rr.properties(steam_phase)
    
    # Mass flow rate (knowing P and Delta h)
    mass_rate  = dh.Pth/(h_steam-h_liquid)
    
    results = {}
    
    # Compute H for each diameter of the pipe
    for i in range(pp.ROWS):
        Pipe_name, D_pipe = pp.get_diameter(i)
        
        # Speed of fluids based on mass flow rate
        v_liquid = mass_rate/rho_liquid*4/cs.pi/D_pipe**2
        v_steam = mass_rate/rho_steam*4/cs.pi/D_pipe**2
    
        # pressure drop due to localized losses, the velocities were assigned by prof's slide    
        loc_losses = total_loc(rho_liquid, rho_steam, v_liquid, v_steam)
    
        # Re
        Re_liquid = rr.Reynolds(liquid_phase, v_liquid, D_pipe)
        Re_steam = rr.Reynolds(steam_phase, v_liquid, D_pipe)
    
        # friction factor
        f_liquid = rr.friction_smooth(Re_liquid)
        f_steam = rr.friction_smooth(Re_steam)
        
        # distributed losses
        Dp_liquid = rr.distributed_loss(dh.L, D_pipe, v_liquid, f_liquid, rho_liquid)
        Dp_steam = rr.distributed_loss(dh.L, D_pipe, v_steam, f_steam, rho_steam)
        
        dist_losses = Dp_liquid + Dp_steam
                        
        h = (loc_losses+dist_losses)/(cs.g*(rho_liquid-rho_steam))
        results[Pipe_name] = [D_pipe, h, v_steam]
        # print(f"{Pipe_name}: h = {h:5f} m, v_s = {v_steam:5f}")
        
    return results
        
# Solution of point b       
def solve_point_b(liquid_phase: Fluid, steam_phase: Fluid) -> dict:
    # Saturated water properties (density, enthalpy)
    rho_liquid, h_liquid, _ = rr.properties(liquid_phase)
    
    # Saturated steam properties (density, enthalpy)
    rho_steam, h_steam, _ = rr.properties(steam_phase)
    
    # Mass flow rate (knowing P and Delta h)
    mass_rate  = dh.Pth/(h_steam-h_liquid)
    
    results = {}
    
    # Compute H for each diameter of the pipe
    for i in range(pp.ROWS):
        Pipe_name, D_pipe = pp.get_diameter(i)
        
        # Speed of fluids based on mass flow rate
        v_liquid = mass_rate/rho_liquid*4/cs.pi/D_pipe**2
        v_steam = mass_rate/rho_steam*4/cs.pi/D_pipe**2
    
        # pressure drop due to localized losses, the velocities were assigned by prof's slide
        loc_losses = total_loc(rho_liquid, rho_steam, v_liquid, v_steam)
    
        # Re
        Re_liquid = rr.Reynolds(liquid_phase, v_liquid, D_pipe)
        Re_steam = rr.Reynolds(steam_phase, v_liquid, D_pipe)
    
        # friction factor
        f_liquid = rr.friction_smooth(Re_liquid)
        f_steam = rr.friction_smooth(Re_steam)
        
        err = 1
        COUNT = 10000
        cc = 0
        h = 0
        
        while err > 1.e-4 and cc < COUNT:
            # distributed losses
            Dp_liquid = rr.distributed_loss(h, D_pipe, v_liquid, f_liquid, rho_liquid)
            Dp_steam = rr.distributed_loss(h, D_pipe, v_steam, f_steam, rho_steam)
            
            dist_losses = Dp_steam + Dp_liquid
                            
            h_new = (loc_losses+dist_losses)/(cs.g*(rho_liquid-rho_steam))
                
            err = 1 - h/h_new
            cc += 1
            h = h_new
        
        results[Pipe_name] = [D_pipe, h, v_steam]
        # print(f"{Pipe_name}: h = {h:5f} m, iter = {cc}, err = {err:2e}")
    
    return results

# Solution of point c       
def solve_point_c(liquid_phase: Fluid, steam_phase: Fluid) -> dict:
    # Saturated water properties (density, enthalpy)
    rho_liquid, h_liquid, _ = rr.properties(liquid_phase)
    
    # Saturated steam properties (density, enthalpy)
    rho_steam, h_steam, _ = rr.properties(steam_phase)
    
    # Mass flow rate (knowing P and Delta h)
    mass_rate  = dh.Pth/(h_steam-h_liquid)
    
    results = {}
    
    # Compute H for each diameter of the pipe
    for i in range(pp.ROWS):
        Pipe_name, D_pipe = pp.get_diameter(i)
        
        # Speed of fluids based on mass flow rate
        v_liquid = mass_rate/rho_liquid*4/cs.pi/D_pipe**2
        v_steam = mass_rate/rho_steam*4/cs.pi/D_pipe**2
    
        # pressure drop due to localized losses, the velocities were assigned by prof's slide
        loc_losses = total_loc(rho_liquid, rho_steam, v_liquid, v_steam)
    
        # Re
        Re_liquid = rr.Reynolds(liquid_phase, v_liquid, D_pipe)
        Re_steam = rr.Reynolds(steam_phase, v_liquid, D_pipe)
    
        # friction factor
        f_liquid = rr.friction_rough(Re_liquid, dh.epsilon, D_pipe)
        f_steam = rr.friction_rough(Re_steam, dh.epsilon, D_pipe)
        
        err = 1
        COUNT = 10000
        cc = 0
        h = 0
        
        while err > 1.e-4 and cc < COUNT:
            # distributed losses
            Dp_liquid = rr.distributed_loss(h, D_pipe, v_liquid, f_liquid, rho_liquid)
            Dp_steam = rr.distributed_loss(h, D_pipe, v_steam, f_steam, rho_steam)
            
            dist_losses = Dp_steam + Dp_liquid
                            
            h_new = (loc_losses+dist_losses)/(cs.g*(rho_liquid-rho_steam))
                
            err = 1 - h/h_new
            cc += 1
            h = h_new
        
        results[Pipe_name] = [D_pipe, h, v_steam]
        # print(f"{Pipe_name}: h = {h:5f} m, iter = {cc}, err = {err:2e}")
    
    return results

def solver() -> list:
    water = Fluid(FluidsList.Water)

    # Liquid objects from pyFluids
    liquid_phase = water.bubble_point_at_pressure(dh.p)
    steam_phase = water.dew_point_at_pressure(dh.p)
    
    # Computing sound speed of steam phase (physical limit)
    rr.SOUND_SPEED = steam_phase.sound_speed

    # Computing the single requests of the problem
    result = []

    result.append(solve_point_a(liquid_phase, steam_phase))
    result.append(solve_point_b(liquid_phase, steam_phase))
    result.append(solve_point_c(liquid_phase, steam_phase))
    
    return result
  