"""
A tool to inflate/deflate bilayers in preparation for backmapping, to avoid stereoclashes in the result.

Usage:
$> python inflate.py input_file output_file scale_factor

"""

import sys
import mdtraj as md
import numpy as np

sys.argv = ["inflate.py", "eq4_lipids.gro", "scaled.gro", -2]

if len(sys.argv) != 4:
    # print("ERROR: Incorrect number of arguments given.\nUsage:\n$> python inflate.py input_file output_file scale_factor")
    raise ValueError("ERROR: Incorrect number of arguments given.\nUsage:\n$> python inflate.py input_file output_file scale_factor")

input_file = sys.argv[1]
output_file = sys.argv[2]

try:
    scale_factor = float(sys.argv[3])
except ValueError:
    raise ValueError("ERROR: Scale factor must be a number")




# Load in bilayer
try:
    bilayer = md.load(input_file)
except FileNotFoundError:
    raise FileNotFoundError("ERROR: Couldn't load file: %s" % input_file)
    
# Get COM of system
system_com = np.mean(bilayer.xyz[0, :, :2], axis=0)

# Get COM of each residue
res_com = np.asarray([np.mean(coords, axis=0) for coords in [bilayer.xyz[0, [i.index for i in res.atoms], :2] for res in bilayer.topology.residues]])

# Get inflated position of each residue COM
# Scale the magnitude of the XY vector from the system COM to the residue COM
# If the scale factor is negative, we want to shrink the system
if scale_factor > 0:
    scaled_com = scale_factor * np.asarray([(r - system_com) for r in res_com])
elif scale_factor < 0:
    scaled_com = (1/scale_factor) * np.asarray([(r - system_com) for r in res_com])

transform_per_atom = np.concatenate([np.tile(x, (n_atoms,1)) for x, n_atoms in zip(scaled_com, [i.n_atoms for i in bilayer.topology.residues])])

# Translate residue coordinates to new COM
# Add to atomic positions the difference in res COM

scaled_coords = np.dstack((bilayer.xyz[:,:,:2] + transform_per_atom[np.newaxis, ...], bilayer.xyz[0,:,2]))

bilayer.xyz = scaled_coords

bilayer.save_gro(output_file)


