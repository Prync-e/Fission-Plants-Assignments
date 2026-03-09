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



def solver() -> list:
    water = Fluid(FluidsList.Water)

    # Liquid objects from pyFluids
    liquid_phase = water.bubble_point_at_pressure(dh.p)
    steam_phase = water.dew_point_at_pressure(dh.p)

    # Computing the single requests of the problem
    result = []
    
    return result
  