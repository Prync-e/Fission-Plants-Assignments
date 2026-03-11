# Main assignment handler
from scripts.solver import solver
from scripts.plotter import plotter


def main():
    results = solver()
    plotter(results)
    

if __name__ == "__main__":
    main()
    