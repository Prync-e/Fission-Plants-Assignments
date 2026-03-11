# Read csv and converts to IS units
import pandas as pd
import scipy.constants as cs

df = pd.read_csv("data/tubi_ansi.csv")
ROWS, _ = df.shape

# Returning the inner diameter of the i-th row of the j-th col and its name
def get_diameter(row: int, col = "STD") -> str|float:
    # Reading from csv
    thickness = df.loc[row][col]
    D_out = df.loc[row]["Outside_Diameter"]
    
    # Inner diameter calculation and convertion do IS units
    D_in = (D_out - 2 * thickness) * cs.inch
    
    # Saving the name of the pipe
    name = df.loc[row]["NPS"]

    return name,D_in
