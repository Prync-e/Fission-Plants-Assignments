# Problem data
Pth = 6e6                   # Wth
DeltaP_Vessel = 1.2e5       # Pa
mass_rate_reference = 3200  # kg/s

# PSC
p_PSC = 75e5                # Pa
Dout_PSC = 16               # in, schedule 100
L_PSC = 16                  # m
epsilon_rel_PSC = 2e-4      # -, epsilon/D
k_loss_PSC = 0.45
k_valve_PSC = 0.12
H1_PSC = 7                  # m
H2_PSC = 3                  # m

# HX1
Dout_HX1 = 19.05e-3         # m
thick_HX1 = 1.24e-3         # m
Ntubes_HX1 = 897            # -, tubes in the heat exchanger
pitch_HX1 = 28.5e-3         # m
Din_shell_HX1 = 1.5         # m
Nbaffles_HX1 = 2
lbaffles_HX1 = 1.6          # m
Ltubes_HX1 = 9.314          # m
k_HX1 = 15                  # W m-1 K-1
Aheader_HX1 = 0.883         # m2
Atransfer_HX1 = 500         # m2
epsilon_rel_HX1 = 1e-4      # -, epsilon/D
FT_HX1 = 0.7                # -, correction factor for MLDT

# ISC
p_ISC = 70e5                # Pa
Dout_ISC = 16               # in, schedule 100
L_ISC = 40                  # m
epsilon_rel_PSC = 2e-4      # -, epsilon/D
k_loss_ISC = 0.45
H_ISC = 10                  # m

# HX2
Dout_HX2 = 25.4e-3          # m
thick_HX2 = 1.24e-3         # m
Ntubes_HX2 = 770            # -, tubes in the heat exchanger
Ltubes_HX2 = 7              # m
Dmanifold_HX2 = 16          # in
k_HX2 = 15                  # W m-1 K-1
Atransfer_HX2 = 430         # m2
epsilon_rel_HX2 = 1e-4      # -, epsilon/D
