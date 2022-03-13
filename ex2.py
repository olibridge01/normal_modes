from numpy import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib as mpl
import os
import sys

# Physical constants and unit conversions
E_J = 4.3597e-18
m_u = 1.6605e-27
r_m = 1e-10
theta_rad = (2*pi)/360
c = 2.9979e10

def plot_energy_surface(geometry_list,molecule):
    """
    Function creates a 3D surface plot of the potential energy.
    Input: Geometry list containing 3-tuples (r, theta, energy); essentially a list of 3D coordinates for the PE surface
    Output: Plots a 3D surface from all the 3D points
    """
    # Extract data from 3-tuples and store them as x, y and z arrays
    x = [item[0] for item in geometry_list]
    y = [item[1] for item in geometry_list]
    z = [item[2] for item in geometry_list]

    # Plot the 3D surface
    ax = plt.axes(projection='3d')
    ax.set_xlabel('r / Angstroms')
    ax.set_ylabel('Theta / degrees')
    ax.set_zlabel('Energy / Hartrees')
    plt.title(molecule + " PE Surface")
    ax.view_init(35,300)
    ax.plot_trisurf(x, y, z, cmap=mpl.cm.coolwarm, alpha = 1.0)
    plt.savefig(molecule+'_PE_surface.pdf')
    plt.show()
    plt.close()
    return None


def load_geometry(file, geometry_list):
    """
    Function extracts the geometry and energy of the triatomic molecule from an input file, and appends these values as a 3-tuple onto an input list of all geometries.
    Input: File to be opened and list to append to
    Output: Appends 3-tuple extracted from the file onto the list
    """
    # Extract the bond length and angle from the file name.
    filename = os.path.basename(file.name)
    r = float(filename[5:9])
    if len(filename) == 22:
        theta = float(filename[14:18])
    else:
        theta = float(filename[14:19])

    # Read the file to find the computed energy for this particular geometry.
    for line in file:
        if "SCF Done" in line:
            line_list = line.split()
            energy = float(line_list[4])

    geometry_list.append((r, theta, energy))
    return None


def get_equilibrium(geometry_list):
    """
    Sorts the list of geometries by energy and stores the equilibrium r, theta and energy as variables to return.
    Input: List of 3-tuples containing geometry data for a molecule
    Output: Variables containing the equilibrium r, theta and energy
    """
    # Sort the geometry list by energy
    sorted_list = geometry_list
    sorted_list.sort(key=lambda x: x[2])

    # Store the equilibrium values (first element in sorted list)
    eq_r, eq_theta, eq_energy = sorted_list[0]

    return eq_r, eq_theta, eq_energy


def print_equilibrium(eq_r, eq_theta, eq_energy, molecule):
    """
    Produces a table from the obtained equilbirium data.
    Input: equilibrium geometry data for a named molecule
    Output: Prints a table showing the equilibrium geometry
    """
    # Print table with equilibrium geometry data
    print()
    print(f"Equilibrium geometry for {molecule}:")

    print("{:<20} {:<20} {:<20}".format("r / Angstroms","Theta / Degrees", "Energy / Hartrees"))
    print("-------------------------------------------------------------")

    print("    {:<20}  {:<15}  {:<15}".format(eq_r, eq_theta, eq_energy))
    print("-------------------------------------------------------------")
    print()


def plot_quadratic_fit(vals, xvar,molecule):
    """
    Input: List of values to plot alongside a quadratic fit
    Output: Scatter plot of values superimposed on plot of quadratic fit
    """
    xvals = [val[0] for val in vals]
    yvals = [val[1] for val in vals]
    p = polyfit(xvals, yvals, 2)

    # Format and display the figure
    fig, ax = plt.subplots()
    xfit = arange(min(xvals),max(xvals),(max(xvals) - min(xvals))/100)
    yfit = polyval(p, xfit)
    plt.title(f"Quadratic fit for small oscillations of {xvar}")
    if xvar == "r":
        ax.set_xlabel(f"{xvar} / Angstrom")
    elif xvar == "Theta":
        ax.set_xlabel(f"{xvar} / Degrees")
    ax.set_ylabel("Energy / Hartrees")
    ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%e'))
    plt.scatter(xvals, yvals, color='black', marker='x')
    plt.plot(xfit, yfit, color='red')
    fig = plt.gcf()
    fig.set_size_inches(12,8)
    plt.savefig(molecule+'_fit_'+xvar+'.pdf')
    plt.show()
    plt.close()


def get_quadratic_fit(vals):
    """
    Input: List of values to fit a quadratic to
    Output: Array containing coefficients a,b and c for quadratic fit ax^2 + bx + c
    """
    xvals = [val[0] for val in vals]
    yvals = [val[1] for val in vals]
    p = polyfit(xvals, yvals, 2)
    return p


def get_normal_mode(vals, xvar, eq_r):
    """
    Input: List of values from the quadratic plots i.e. the five lowest energy values for constant r or theta
    Output: Frequency of the corresponding normal mode
    """
    # Create copy of list of values in SI units using the conversion variables defined above
    SI_vals = []
    for index, item in enumerate(vals):
        SI_vals.append([item[0], item[1]])
        SI_vals[index][1] *= E_J
        if xvar == "r":
            SI_vals[index][0] *= r_m
        elif xvar == "Theta":
            SI_vals[index][0] *= theta_rad

    # Fit the SI_data to a quadratic
    polynom = get_quadratic_fit(SI_vals)

    # Compute the relevant normal mode using formula from handout
    if xvar == "r":
        freq = (1 / (2 * pi * c)) * sqrt((2 * polynom[0]) / (2 * m_u))
        return freq
    elif xvar == "Theta":
        freq = (1 / (2 * pi * c)) * sqrt((2 * polynom[0]) / (power(eq_r*r_m, 2) * 0.5 * m_u))
        return freq


def user_help():
    """
    Print a help guide for the user to know what to input into the command line.
    """
    print()
    print("User input: python3 ex2.py /[valid Gaussian directory name]/\n"
          "NB: Make sure the Gaussian output file directory exists and is placed in the same directory as the python file")
    print()


# Create empty list to add geometry data
geoms = []
filepath = os.path.dirname(os.path.abspath(__file__))

if len(sys.argv) != 2 or os.path.exists(filepath + sys.argv[1]) == False:
    user_help()
    quit()

# Create path for Gaussian files
path = filepath + sys.argv[1]

# Extract geometry data
molecule = os.listdir(path)[0][0:3]

for filename in os.listdir(path):
    f = open(f'{path}{filename}', 'r')
    load_geometry(f, geoms)

# Compute and print equilibrium geometry
eq_r, eq_theta, eq_energy = get_equilibrium(geoms)
print_equilibrium(eq_r, eq_theta, eq_energy, molecule)

# Plot PE surface
plot_energy_surface(geoms, molecule)

# Create list of 5 values around equilibrium at constant r and theta for the quadratic fits
theta_vals = []
for item in geoms:
    if item[0] == eq_r:
        theta_vals.append((item[1],item[2]))
    if len(theta_vals) == 5:
        break

r_vals = []
for item in geoms:
    if item[1] == eq_theta:
        r_vals.append((item[0], item[2]))
    if len(r_vals) == 5:
        break

# Plot the quadratic fits and compute and print the normal modes
plot_quadratic_fit(r_vals, "r", molecule)
plot_quadratic_fit(theta_vals, "Theta", molecule)

print(f"Symmetric stretch: v1 = {round(get_normal_mode(r_vals, 'r', eq_r), 1)} cm-1")
print(f"Symmetric bend: v2 = {round(get_normal_mode(theta_vals, 'Theta', eq_r), 1)} cm-1")
print()
