from pyfluids import Fluid, FluidsList
import data.assignment_data as dh

def solver():
    water = Fluid(FluidsList.Water)

    # Liquid objects from pyFluids
    liquid_phase = water.bubble_point_at_pressure(dh.p)
    steam_phase = water.dew_point_at_pressure(dh.p)

    # Saturated water properties (density, enthalpy)
    rho_liquid = liquid_phase.density
    h_liquid = liquid_phase.enthalpy
    
    # Saturated steam properties (density, enthalpy)
    rho_steam = steam_phase.density
    h_steam = steam_phase.enthalpy
    
    # Mass flow rate (knowing P and Delta h)
    mass_rate  = dh.Pth/(h_steam-h_liquid)
    
    # Volume flow rates
    Q_liquid = mass_rate/rho_liquid
    Q_steam = mass_rate/rho_steam
    
    # Dato un D dobbiamo trovare h tale che V sia compatibile con Q
    
    # Considerato che non c'è più un tratto orizzontale L=h
    
    # Punto b ma con roughness

    print(rho_liquid)
    print(rho_steam)
    print(mass_rate)
    