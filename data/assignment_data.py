import scipy.constants as cs

# Problem data
H = 168*cs.inch                 # m, active height
lambda_tr = 0.29e-2             # m, transport length
Dc = lambda_tr/3                # -, diffusion coefficient in the core
Dr = 0.16                       # -, diffusion coefficient in the reflector
Lr = 2.85e-2                    # m, diffusion length in the reflector
delta = Dc/Dr*Lr                # m, reflector savings
heat_in_fuel = 0.974            # -, heat generated in fuel
Fq = 2.6                        # -, heat flux hot channel factor
P = 3400e6                      # W, reactor core heat output

# Fuel pellets
N_rods = 41448                  # -, number of fuel rods
D_fuel_pellet = 0.3225*cs.inch  # m, fuel pellet diameter
D_fuel_road = 0.374*cs.inch     # m, fuel road outer diameter
t_gap = 0.0065*cs.inch          # m, gap between fuel and cladding
t_cladding = 0.0225*cs.inch     # m, cladding thickness

# Coolant flow
mass_flow_coolant = 106.8*cs.lb/3600    # kg/s, effective mass flow rate in the core
A_flow_coolant = 41.8*cs.foot**2        # m^2, effective flow area
p_coolant = 2250*cs.psi         # Pa, coolant pressure
T_in_coolant = cs.convert_temperature(535.0, 'F', 'C')  # °C, inlet coolant temperature
pitch = 0.496*cs.inch           # m, core pitch

# Numerics
n_points = 11
