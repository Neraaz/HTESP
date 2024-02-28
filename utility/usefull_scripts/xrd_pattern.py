from pymatgen.core import structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import matplotlib.pyplot as plt
import sys
filename = sys.argv[1]

# Load crystal structure from CIF file
structure = structure.Structure.from_file(filename)

# Set up XRD calculator
xrd_calculator = XRDCalculator()

# Calculate XRD pattern
xrd_pattern = xrd_calculator.get_pattern(structure)

# Plot the XRD pattern
plt.figure(figsize=(8, 6))
xrd_calculator.show_plot(structure, two_theta_range=(0, 90), annotate_peaks=True, ax=plt.gca())
plt.title("XRD Pattern")
plt.xlabel(r"2$\theta$ (degrees)")
plt.ylabel("Intensity (a.u.)")
plt.savefig("xrd.pdf")

