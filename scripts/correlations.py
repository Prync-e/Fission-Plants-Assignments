from pyfluids import Fluid
import numpy as np

# Fluid properties function
def properties(fluid: Fluid) -> float|float|float:
    """Produces properties of fluids from the pyFluid library

    Args:
        fluid (Fluid): 

    Returns:
        float|float|float: fluid density, specific enthalpy, dynamic viscosity
    """
    return fluid.density, fluid.enthalpy, fluid.dynamic_viscosity

# Reynolds number
def Reynolds(fluid: Fluid, vel: float, diam: float) -> float:
    rho, _, mu = properties(fluid)
    Re = rho * vel * diam / mu
    return Re

def Reynolds_massrate(fluid: Fluid, mass_rate: float, diam: float) -> float:
    _,_,mu = properties(fluid)
    Re = 4 * mass_rate / (mu * diam * np.pi)
    return Re

# Prandtl number
def Prandtl(fluid: Fluid):
    return fluid.prandtl

# Nusselt number
def Nusselt(fluid: Fluid, mass_rate: float, diam: float, epsilon_rel: float) -> float:
    Re = Reynolds_massrate(fluid, mass_rate, diam)
    Pr = Prandtl(fluid)
    f = friction_rough(Re, epsilon_rel)
    # Gnielinski correlation
    Nu = (f/8*(Re-1000)*Pr)/(1+12.7*pow(f/8,0.5)*(pow(Pr,2/3)-1)) 
    return Nu

# Friction coefficient for smooth pipes
def friction_smooth(Re: float) -> float:
    return 64/Re if Re < 3000 else 0.316/pow(Re, 0.25)

# Friction coeeficient for rough pipes
def friction_rough(Re: float, epsilon_rel: float) -> float:
    # Colebrook equation
    err = 1
    f = 1
    cc = 0
    COUNT = 1000
    
    while err > 1.e-4 and cc < COUNT:
        f_new = pow(-2*np.log10(2.51/Re/np.sqrt(f) + epsilon_rel/3.71),-2)

        err = 1 - f_new/f
        f = f_new
        cc += 1
        
    return 64/Re if Re < 3000 else f

# Localized losses
def localized_loss(k: float, rho: float, v: float) -> float:
    return 0.5*k*rho*v**2

# Distributed losses
def distributed_loss(l: float, d: float, v: float, f: float, rho: float) -> float:
    return 0.5*f*l/d*rho*v**2

# Equivalent shell-side diameter for triangular lattice (HX)
def shell_D_equiv(pitch: float, diam: float) -> float:
    return 2*np.sqrt(3)*pitch**2 / (np.pi*diam) - diam

# Heat transfer coefficent HX1
def h_HX1(Re: float, pitch: float, diam: float, Pr: float, kth_fluid: float) -> float:
    return 0.351*pow(Re, 0.55)*kth_fluid/shell_D_equiv(pitch, diam)*pow(Pr, 1/3)

# Shell side flow area
def shell_flow_area(Dshell: float, pitch: float, diam: float, lbaffles: float):
    return Dshell/pitch*(pitch-diam)*lbaffles

# Local loss coeefficient for HX1 
def local_loss_HX1(Re: float, dshell: float, dtube: float, pitch: float, Nb = 2.) -> float:
    return 8*0.227/pow(Re, 0.193)*dshell/shell_D_equiv(pitch, dtube)*(Nb+1)

# McAdams correlation (inverse to get DeltaTsat)
def inverse_McAdams(q2: float) -> float:
    return pow(q2/2.257, 1/3.86)    # [°C - K]
