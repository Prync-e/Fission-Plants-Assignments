from pyfluids import Fluid
import numpy as np

SOUND_SPEED = 0

def properties(fluid: Fluid) -> float|float|float:
    """Produces properties of fluids from the pyFluid library

    Args:
        fluid (Fluid): 

    Returns:
        float|float|float: fluid density, specific enthalpy, dynamic viscosity
    """
    return fluid.density, fluid.enthalpy, fluid.dynamic_viscosity

def Reynolds(fluid: Fluid, vel: float, diam: float) -> float:
    rho, _, mu = properties(fluid)
    Re = rho * vel * diam / mu
    return Re

def friction_smooth(Re: float) -> float:
    return 64/Re if Re < 3000 else 0.316/pow(Re, 0.25)

def friction_rough(Re: float, epsilon: float, diam: float) -> float:
    # Colebrook equation
    err = 1
    f = 1
    cc = 0
    COUNT = 1000
    
    while err > 1.e-3 and cc < COUNT:
        f_new = pow(-2*np.log10(2.51/Re/np.sqrt(f) + epsilon/diam/3.71),-2)

        err = 1 - f_new/f
        f = f_new
        cc += 1
        
    return 64/Re if Re < 3000 else f

def localized_loss(k: float, rho: float, v: float) -> float:
    return 0.5*k*rho*v**2

def distributed_loss(l: float, d: float, v: float, f: float, rho: float) -> float:
    return 0.5*f*l/d*rho*v**2
