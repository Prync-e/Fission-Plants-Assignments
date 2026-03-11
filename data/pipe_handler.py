# Read csv and converts to IS units
import pandas as pd
import scipy.constants as cs

df = pd.read_csv("data/tubi_ansi.csv")
ROWS, _ = df.shape

# Returning the inner diameter of the i-th row of the j-th col and its name
def get_diameter(row: int, col = "100") -> str|float:
    # Reading from csv
    thickness = df.loc[row][col]
    D_out = df.loc[row]["Outside_Diameter"]
    
    # Inner diameter calculation and convertion do IS units
    D_in = (D_out - 2 * thickness) * cs.inch
    
    # Saving the name of the pipe
    name = df.loc[row]["NPS"]

    return name,D_in

# Returning the row of the pipe, given its outer diameter in inches
def get_row(D_out: float) -> int:
    for r in range(ROWS):
        if float(df.loc[r]["Outside_Diameter"]) == D_out:
            break
        
    return r
