import matplotlib.pyplot as plt
import data.assignment_data as dh

# Handles the points results (dictionaries)
def split_results(res: dict) -> list|list|list:
    names = res.keys()
    diam = [val[0] for val in res.values()]
    h = [val[1] for val in res.values()]
    
    return names, diam, h

# Specific plots for point a
def plots_point_a(point_a: dict):
    Pipe_names, D_pipes, computed_h = split_results(point_a)   
    
    # Look for the first accepteable h (h<L)
    index_h = 0
    for h in computed_h:
        if h < dh.L:
            index_h = computed_h.index(h) - 1
            break
    
    # Slicing the result lists
    Pipe_names = list(Pipe_names)[index_h:]
    D_pipes = D_pipes[index_h:]
    computed_h = computed_h[index_h:]
    
    plt.plot(D_pipes, computed_h, "ro-", label="Minimum h")
    plt.xticks(D_pipes, Pipe_names)
    plt.axhline(y=dh.L, color='k', linestyle='--', linewidth=2, label="Max pipe lenght")
    plt.xlabel("D [in]")
    plt.ylabel("h [m]")
    plt.legend()
    plt.grid()
    plt.show()


def plotter(results: list):
    point_a, point_b, point_c = results
    plots_point_a(point_a)
    
    