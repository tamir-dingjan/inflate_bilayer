"""
A tool to inflate/deflate bilayers in preparation for backmapping, to avoid stereoclashes in the result.

"""


import mdtraj as md
import numpy as np

scale_factor = 4

# Load in bilayer
bilayer = md.load("eq4_lipids.gro")

# Get COM of system
system_com = np.mean(bilayer.xyz[0, :, :2], axis=0)

# Get COM of each residue
res_com = np.asarray([np.mean(coords, axis=0) for coords in [bilayer.xyz[0, [i.index for i in res.atoms], :2] for res in bilayer.topology.residues]])

# Get inflated position of each residue COM
# Scale the magnitude of the XY vector from the system COM to the residue COM

scaled_com = scale_factor * np.asarray([(r - system_com) for r in res_com])

transform_per_atom = np.concatenate([np.tile(x, (n_atoms,1)) for x, n_atoms in zip(scaled_com, [i.n_atoms for i in bilayer.topology.residues])])

# Translate residue coordinates to new COM
# Add to atomic positions the difference in res COM

scaled_coords = np.dstack((bilayer.xyz[:,:,:2] + transform_per_atom[np.newaxis, ...], bilayer.xyz[0,:,2]))

bilayer.xyz = scaled_coords

bilayer.save_gro("scaled_%s.gro" % scale_factor)


