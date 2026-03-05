from pyfluids import Fluid

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

