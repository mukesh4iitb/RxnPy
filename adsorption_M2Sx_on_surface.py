#!/usr/bin/env python
# coding: utf-8

# In[4]:


import re
from ase import Atoms
import numpy as np
from typing import Union, Optional
import ase.io
from scipy.optimize import minimize_scalar


# In[28]:


def principle_axis_coordinates(filename):
    """
    Here it has not been used.
    """
    mol = ase.io.read("{}".format(filename), format='xyz')
    I, eigvecs = mol.get_moments_of_inertia(vectors=True)
    
    r = mol.positions - mol.get_center_of_mass()
    rotated = (eigvecs @ r.T).T
    mol.positions = rotated
    ase.io.write("{}_principle_axis.xyz".format(filename.split(".xyz")[0]), mol, format='xyz')
    return mol


# In[3]:


def rotation_matrix(axis: list, theta: float):
    """Rodrigues’ rotation formula."""
    axis = axis / np.linalg.norm(axis)
    a = np.cos(theta/2)
    b, c, d = -axis*np.sin(theta/2)
    return np.array([
        [a*a+b*b-c*c-d*d, 2*(b*c-a*d),     2*(b*d+a*c)],
        [2*(b*c+a*d),     a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
        [2*(b*d-a*c),     2*(c*d+a*b),     a*a+d*d-b*b-c*c]
    ])



def rotation_matrix_vector_to_vector(mol: Atoms, a: int, b: int, target_vector: list):
    """Find rotation matrix that aligns v to u."""
    v = mol.positions[b]-mol.positions[a]
    u = target_vector
    v = v / np.linalg.norm(v)
    u = u / np.linalg.norm(u)
    cross = np.cross(v, u)
    dot = np.dot(v, u)
    if np.isclose(dot, 1.0):  # already aligned
        return np.eye(3)
    if np.isclose(dot, -1.0):  # opposite direction
        # 180° rotation around any axis perpendicular to v
        perp = np.array([1,0,0]) if not np.allclose(v,[1,0,0]) else np.array([0,1,0])
        axis = np.cross(v, perp)
        axis /= np.linalg.norm(axis)
        return rotation_matrix(axis, np.pi)

    s = np.linalg.norm(cross)
    kmat = np.array([[0, -cross[2], cross[1]],
                     [cross[2], 0, -cross[0]],
                     [-cross[1], cross[0], 0]])
    R = np.eye(3) + kmat + kmat @ kmat * ((1 - dot) / (s**2))
    return R


# In[14]:


#def rotated_point(t_pt0, u, theta):
#    """Rotate t_pt0 around axis u by angle theta (Rodrigues formula)."""
#    return (
#        t_pt0 * np.cos(theta)
#        + np.cross(u, t_pt0) * np.sin(theta)
#        + u * np.dot(u, t_pt0) * (1 - np.cos(theta)))
#
#def objective(theta, t_pt0, u, slab):
#    """Negative of z-distance, since we maximize."""
#    t_pt0_rot = rotated_point(t_pt0, u, theta)
#    t_pt0_dist = np.dot(t_pt0_rot - slab.positions[0], [0, 0, 1])
#    return -t_pt0_dist  # minimize negative → maximize positive

# Optimize over theta in [0, 2π]
def find_best_rotation_such_that_tpt_at_top(t_pt0, u, slab):
    res = minimize_scalar(objective, args=(t_pt0, u, slab), bounds=(0, 2*np.pi), method='bounded')
    print("Optimal θ (radians):", res.x)
    return res.x 


# In[15]:


# the following is generalization of one point having maximum distance. Done by GPT. So, check each line befor implimenting. 
# Here every points will have maximum distance.

# import numpy as np
# from scipy.optimize import minimize_scalar

# def rotate_molecule(mol_positions, a_idx, b_idx, theta):
#     """
#     Rotate all atoms in mol_positions around axis defined by atoms a_idx-b_idx
#     by angle theta (Rodrigues formula).
#     """
#     A = mol_positions[a_idx]
#     B = mol_positions[b_idx]
#     u = B - A
#     u = u / np.linalg.norm(u)  # normalize axis

#     # build rotation matrix (Rodrigues)
#     K = np.array([[0, -u[2], u[1]],
#                   [u[2], 0, -u[0]],
#                   [-u[1], u[0], 0]])
#     R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)

#     # translate so axis passes through origin, rotate, then translate back
#     shifted = mol_positions - A
#     rotated = (R @ shifted.T).T + A
#     return rotated

# def objective(theta, mol_positions, a_idx, b_idx, slab_positions):
#     """
#     Objective: maximize distance from surface (z=0 plane).
#     Equivalent: minimize negative of min z-coordinate.
#     """
#     rotated = rotate_molecule(mol_positions, a_idx, b_idx, theta)
#     z_coords = rotated[:, 2]

#     # Example objective: maximize minimum height of all atoms above slab
#     min_dist = np.min(z_coords - np.max(slab_positions[:, 2]))
#     return -min_dist  # negate for minimization

# def find_best_rotation(mol_positions, a_idx, b_idx, slab_positions):
#     """
#     Find optimal angle θ that maximizes distance from slab surface.
#     """
#     res = minimize_scalar(objective, args=(mol_positions, a_idx, b_idx, slab_positions),
#                           bounds=(0, 2*np.pi), method='bounded')
#     print("Optimal θ (radians):", res.x)
#     return res.x


# In[73]:


def orient_and_adsorption0(slab: Atoms, absorption_site: Union[int, list], target: list, mol: Atoms, a: int, b: int, c: Optional[int] = None, h: Optional[float] = 1, slab_name: str = "slab", mol_name: str = "mol"):

    mol.positions = mol.positions -mol.get_center_of_mass()
    if c is None:
        R = rotation_matrix_vector_to_vector(mol, a, b, target)
    else:
        R = rotation_matrix_plane_to_plane(mol, a, b, c, target)
    
    new_cart = np.dot(R, mol.positions.T).T
    mol.positions = new_cart

    if c is None:
        lowest_atom = np.min(mol[[a, b]].positions[:, 2])
    else:
        lowest_atom = np.min(mol[[a, b, c]].positions[:, 2])

    print(lowest_atom) 

    if isinstance(absorption_site, int):
        mol.positions = mol.positions + slab.positions[absorption_site] -lowest_atom + [0, 0, h]
    else:
        mol.positions = mol.positions + absorption_site -lowest_atom + [0, 0, h]

    M21=mol.positions[b]-mol.positions[a]
    u = M21/np.linalg.norm(M21)
    S_atoms = mol.positions[2:]
    S0 = S_atoms -mol.positions[a]
    t_pt = (mol.positions[5]+mol.positions[7])/2
    print(t_pt)
    t_pt0 = t_pt - mol.positions[a]
    
    theta = find_best_rotation_such_that_tpt_at_top(t_pt0, u, slab)
    print("theta:", theta)
    
    # Optimize over theta in [0, 2π]
    #res = minimize_scalar(objective, bounds=(0, 2*np.pi), method='bounded')
    #print("theta:", res.x)
    #theta = res.x #5.3347972047602505
    t_pt0_rot = t_pt0*np.cos(theta) + np.cross(u, t_pt0)*np.sin(theta) + u*np.dot(u, t_pt0)*(1-np.cos(theta))
    #t_pt0_dist = np.dot(t_pt0_rot - slab.positions[44], [0, 0, 1])
    print('t_pt:',t_pt0_rot+mol.positions[a])
    
    S0_new_pos = []
    for s0 in S0:
        s0_rot = s0*np.cos(theta) + np.cross(u, s0)*np.sin(theta) + u*np.dot(u, s0)*(1-np.cos(theta))
        print('S atom:',s0_rot+mol.positions[a])
        S0_new_pos.append(s0_rot+mol.positions[a])
    #print('Distance from ref point to t_pt after rotation:', t_pt0_dist)

    # S_atoms_rot = S0*np.cos(theta) + np.cross(u, S0)*np.sin(theta) + u*np.dot(u, S0)*(1-np.cos(theta))
    mol.positions[2:] = S0_new_pos

    print("new rotated mol's positions:", mol.positions)

    slab_with_doped_mol = slab + mol 
 #ase.io.write("slab_with_doped_mol.vasp", slab_with_doped_mol, format='vasp')
    ase.io.write("{}_{}.vasp".format(slab_name, mol_name), slab_with_doped_mol, format='vasp')


# In[74]:


#for poly in [ "K2S6.xyz", "K2S8.xyz"]:
#    #mol = ase.io.read(poly, format='xyz')
#    mol = principle_axis_coordinates(poly)
#    slab = ase.io.read('gCN_melon.vasp')
#    slab.cell[2][2]=20
#    slab.center(axis=(2))
#    pt=(slab.positions[44]+slab.positions[45])/2
#    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, h=1.75, slab_name="gCN", mol_name=poly.split('.xyz')[0])
#



def rotated_point(theta, t_pt0, u):
    """Rotate t_pt0 around axis u by angle theta (Rodrigues formula)."""
    u = u / np.linalg.norm(u)  # make sure u is normalized
    return (
        t_pt0 * np.cos(theta)
        + np.cross(u, t_pt0) * np.sin(theta)
        + u * np.dot(u, t_pt0) * (1 - np.cos(theta))
    )


def objective(theta, t_pt0, u, slab, mode):
    """Negative of z-distance, since we maximize."""
    t_pt0_rot = rotated_point(theta, t_pt0, u)
    if mode == 'perpendicular':
        t_pt0_dist = np.dot(t_pt0_rot - slab.positions[0], [0, 0, 1])
    elif mode == 'parallel':
        t_pt0_dist = np.linalg.norm(np.dot(t_pt0_rot - slab.positions[0], [1, 0, 0]) + np.dot(t_pt0_rot - slab.positions[0], [0, 1, 0]))
    else:
        print("mode is not specified.")
    #np.linalg.norm(np.dot(pt - ref, [1, 0, 0]) + np.dot(pt - ref, [0, 1, 0]))
    return -t_pt0_dist  # minimize negative → maximize positive


#def objective(theta, t_pt0, u, slab):
#    """Negative of z-distance, since we maximize."""
#    t_pt0_rot = rotated_point(theta, t_pt0, u)
#    t_pt0_dist = np.dot(t_pt0_rot - slab.positions[0], [0, 0, 1])
#    return t_pt0_dist  # minimize negative → maximize positive


def orient_and_adsorption1(slab: Atoms, absorption_site: Union[int, list], target: list, mol: Atoms, a: int, b: int, c: Optional[int] = None, h: Optional[float] = 1, slab_name: str = "slab", mol_name: str = "mol"):

    mol.positions = mol.positions -mol.get_center_of_mass()
    if c is None:
        R = rotation_matrix_vector_to_vector(mol, a, b, target)
    else:
        R = rotation_matrix_plane_to_plane(mol, a, b, c, target)
    
    new_cart = np.dot(R, mol.positions.T).T
    mol.positions = new_cart

#     if c is None:
#         lowest_atom = np.min(mol[[a, b]].positions[:, 2])
#     else:
#         lowest_atom = np.min(mol[[a, b, c]].positions[:, 2])

#     print("lowest pt:", lowest_atom) 

   

    M21=mol.positions[b]-mol.positions[a]
    u = M21/np.linalg.norm(M21)
    S_atoms = mol.positions[2:]
    S0 = S_atoms -mol.positions[a]
    
    # if isinstance(t_pt, int):
    #     t_pt = mol.positions[t_pt]
    # else:
    #     t_pt = mol.positions[t_pt].mean()
        
    t_pt = (mol.positions[7]+mol.positions[9])/2
    print(t_pt)
    t_pt0 = t_pt - mol.positions[a]

    # Optimize over theta in [0, 2π]
    res = minimize_scalar(objective,  args=(t_pt0, u, slab), bounds=(0, 2*np.pi), method='bounded')
    print("theta:", res.x)
    theta = res.x #5.3347972047602505
    t_pt0_rot = t_pt0*np.cos(theta) + np.cross(u, t_pt0)*np.sin(theta) + u*np.dot(u, t_pt0)*(1-np.cos(theta))
    #t_pt0_dist = np.dot(t_pt0_rot - slab.positions[44], [0, 0, 1])
    print('t_pt:',t_pt0_rot+mol.positions[a])
    
    S0_new_pos = []
    for s0 in S0:
        s0_rot = s0*np.cos(theta) + np.cross(u, s0)*np.sin(theta) + u*np.dot(u, s0)*(1-np.cos(theta))
        print('S atom:', s0_rot+mol.positions[a])
        S0_new_pos.append(s0_rot+mol.positions[a])
    #print('Distance from ref point to t_pt after rotation:', t_pt0_dist)

    # S_atoms_rot = S0*np.cos(theta) + np.cross(u, S0)*np.sin(theta) + u*np.dot(u, S0)*(1-np.cos(theta))
    mol.positions[2:] = S0_new_pos

    print("min-height:", np.min(mol.positions[:, 2]))
    
    min_h = np.min(mol.positions[:, 2])

    if isinstance(absorption_site, int):
        mol.positions = mol.positions - mol.get_center_of_mass() + slab.positions[absorption_site] - [0, 0, min_h] + [0, 0, h]
    else:
        mol.positions = mol.positions - mol.get_center_of_mass() + absorption_site  - [0, 0, min_h] + [0, 0, h]
    
    print("COM", mol.get_center_of_mass())
    slab_with_doped_mol = slab + mol 
 #ase.io.write("slab_with_doped_mol.vasp", slab_with_doped_mol, format='vasp')
    ase.io.write("{}_{}.vasp".format(slab_name, mol_name), slab_with_doped_mol, format='vasp')


def orient_and_adsorption(slab: Atoms, absorption_site: Union[int, list], target: list, mol: Atoms, a: int, b: int, t_pt: Union[int, list], c: Optional[int] = None, h: Optional[float] = 1, slab_name: str = "slab", mol_name: str = "mol", mode: str = 'perpendicular'):

    mol.positions = mol.positions - mol.get_center_of_mass()
    if c is None:
        R = rotation_matrix_vector_to_vector(mol, a, b, target)
    else:
        R = rotation_matrix_plane_to_plane(mol, a, b, c, target)
    
    new_cart = np.dot(R, mol.positions.T).T
    mol.positions = new_cart

#     if c is None:
#         lowest_atom = np.min(mol[[a, b]].positions[:, 2])
#     else:
#         lowest_atom = np.min(mol[[a, b, c]].positions[:, 2])
#     print("lowest pt:", lowest_atom)

    M21=mol.positions[b]-mol.positions[a]
    u = M21/np.linalg.norm(M21)
    S_atoms = mol.positions[2:]
    S0 = S_atoms -mol.positions[a]
    
    if isinstance(t_pt, int):
        t_pt = mol.positions[t_pt]
    else:
        t_pt = mol.positions[t_pt].mean(axis=0)
        
    #t_pt = (mol.positions[7]+mol.positions[9])/2
    print(t_pt)
    t_pt0 = t_pt - mol.positions[a]

    # Optimize over theta in [0, 2π]
    res = minimize_scalar(objective,  args=(t_pt0, u, slab, mode), bounds=(0, 2*np.pi), method='bounded')
    print("theta:", res.x)
    theta = res.x #5.3347972047602505
    t_pt0_rot = t_pt0*np.cos(theta) + np.cross(u, t_pt0)*np.sin(theta) + u*np.dot(u, t_pt0)*(1-np.cos(theta))
    #t_pt0_dist = np.dot(t_pt0_rot - slab.positions[44], [0, 0, 1])
    print('t_pt:',t_pt0_rot+mol.positions[a])
    
    S0_new_pos = []
    for s0 in S0:
        s0_rot = s0*np.cos(theta) + np.cross(u, s0)*np.sin(theta) + u*np.dot(u, s0)*(1-np.cos(theta))
        print('S atom:', s0_rot+mol.positions[a])
        S0_new_pos.append(s0_rot+mol.positions[a])
    #print('Distance from ref point to t_pt after rotation:', t_pt0_dist)

    # S_atoms_rot = S0*np.cos(theta) + np.cross(u, S0)*np.sin(theta) + u*np.dot(u, S0)*(1-np.cos(theta))
    mol.positions[2:] = S0_new_pos

    #print("min-height:", np.min(mol.positions[:, 2]))
    
    #min_h = np.min(mol.positions[:, 2])

    if isinstance(absorption_site, int):
        mol.positions = mol.positions - mol.get_center_of_mass() + slab.positions[absorption_site] #- [0, 0, min_h] + [0, 0, h]
    else:
        mol.positions = mol.positions - mol.get_center_of_mass() + absorption_site # - [0, 0, min_h] + [0, 0, h]
    
    min_h = np.min(mol.positions[:, 2]-slab.positions[0][2])
    print("min_h:", min_h)
    mol.positions = mol.positions - [0, 0, min_h] + [0, 0, h]
    print("COM", mol.get_center_of_mass())
    slab_with_doped_mol = slab + mol 
 #ase.io.write("slab_with_doped_mol.vasp", slab_with_doped_mol, format='vasp')
    ase.io.write("{}_{}.vasp".format(slab_name, mol_name), slab_with_doped_mol, format='vasp')






## For Li2Sx

def doping_Li2S(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_Li2S.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, 2, h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="Li2S")



def doping_Li2S2(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_Li2S2.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, [2, 3], h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="Li2S2")

def doping_Li2S4(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_Li2S4.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, [2, 3], h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="Li2S4")

def doping_Li2S6(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_Li2S6.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, [3, 4], h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="Li2S6")



def doping_Li2S8(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_Li2S8.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, [3, 4], h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="Li2S8")

### For Na2Sx


def doping_Na2S(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_Na2S.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, 2, h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="Na2S")



def doping_Na2S2(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_Na2S2.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, [2, 3], h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="Na2S2")

def doping_Na2S4(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_Na2S4.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, [2, 3], h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="Na2S4")

def doping_Na2S6(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_Na2S6.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, [3, 4], h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="Na2S6")



def doping_Na2S8(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_Na2S8.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, [3, 4], h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="Na2S8")

## For K2Sx


def doping_K2S(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_K2S.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, 2, h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="K2S")



def doping_K2S2(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_K2S2.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, [2, 3], h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="K2S2")

def doping_K2S4(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_K2S4.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, [2, 5], h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="K2S4")

def doping_K2S6(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_K2S6.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, [6, 7], h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="K2S6")

def doping_K2S8(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_K2S8.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, [7, 9], h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="K2S8")

def doping_S8(surface_filename, absorption_site, target_vector, h=1.5):
    mol  = ase.io.read("CONTCAR_S8.vasp", format='vasp')
    slab = ase.io.read(surface_filename)
    slab.cell[2][2]=20
    slab.center(axis=(2))
    pt=absorption_site
    target_vector = target_vector
    print("pt", pt)
    orient_and_adsorption(slab, pt, target_vector, mol, 0, 1, 5, h=h, slab_name=surface_filename.split(".vasp")[0], mol_name="S8", mode='parallel')



Li_func_map = {
     "Li2S": doping_Li2S,
    "Li2S2": doping_Li2S2,
    "Li2S4": doping_Li2S4,
    "Li2S6": doping_Li2S6,
    "Li2S8": doping_Li2S8,
      "S8": doping_S8
}

for i in ["Li2S", "Li2S2", "Li2S4", "Li2S6", "Li2S8", "S8"]:

    func = Li_func_map[str(i)]
#
    filename='melon.vasp'
    slab = ase.io.read(filename)
    absorption_site = (slab.positions[54] + slab.positions[55])/2
    target_vector = slab.positions[54] - slab.positions[55]
    func(filename, absorption_site, target_vector, h=1.75)
#
#    filename='gCN.vasp'
#    slab = ase.io.read(filename)
#    absorption_site = (slab.positions[2] + slab.positions[39]) / 2
#    target_vector = slab.positions[2] - slab.positions[39]
#    func(filename, absorption_site, target_vector)
#

#    filename='graphene_551.vasp'
#    slab = ase.io.read(filename)
#    absorption_site = slab.positions[25]
#    target_vector = slab.positions[34] - slab.positions[26]
#    func(filename, absorption_site, target_vector, h=1.8)
#    
#    filename='N3_graphene_551.vasp'
#    slab = ase.io.read(filename)
#    absorption_site = (slab.positions[46] + slab.positions[47]+ slab.positions[48]) / 3
#    target_vector = slab.positions[47] - slab.positions[31]
#    func(filename, absorption_site, target_vector)
#    
#    filename='N4_graphene_551.vasp'
#    slab = ase.io.read(filename)
#    absorption_site = (slab.positions[44] + slab.positions[45]+ slab.positions[46] + slab.positions[47]) / 4
#    target_vector = slab.positions[47] - slab.positions[44]
#    func(filename, absorption_site, target_vector)
#
#
#
#
#Na_func_map = {
#     "Na2S": doping_Na2S,
#    "Na2S2": doping_Na2S2,
#    "Na2S4": doping_Na2S4,
#    "Na2S6": doping_Na2S6,
#    "Na2S8": doping_Na2S8,
#      "S8": doping_S8
#}
#
#for i in ["Na2S", "Na2S2", "Na2S4", "Na2S6", "Na2S8", "S8"]:
#    func = Na_func_map[str(i)]
#
#    filename='gCN_melon.vasp'
#    slab = ase.io.read(filename)
#    absorption_site = (slab.positions[7] + slab.positions[29]) / 2
#    target_vector = slab.positions[7] - slab.positions[29]
#    func(filename, absorption_site, target_vector)
#
#    filename='gCN.vasp'
#    slab = ase.io.read(filename)
#    absorption_site = (slab.positions[2] + slab.positions[39]) / 2
#    target_vector = slab.positions[2] - slab.positions[39]
#    func(filename, absorption_site, target_vector)
#
#    filename='graphene_551.vasp'
#    slab = ase.io.read(filename)
#    absorption_site = slab.positions[25]
#    target_vector = slab.positions[34] - slab.positions[26]
#    func(filename, absorption_site, target_vector, h=2.1)
#    
#    filename='N3_graphene_551.vasp'
#    slab = ase.io.read(filename)
#    absorption_site = (slab.positions[46] + slab.positions[47]+ slab.positions[48]) / 3
#    target_vector = slab.positions[47] - slab.positions[31]
#    func(filename, absorption_site, target_vector, h=1.7)
#    
#    filename='N4_graphene_551.vasp'
#    slab = ase.io.read(filename)
#    absorption_site = (slab.positions[44] + slab.positions[45]+ slab.positions[46] + slab.positions[47]) / 4
#    target_vector = slab.positions[47] - slab.positions[44]
#    func(filename, absorption_site, target_vector)
#
#K_func_map = {
#     "K2S": doping_K2S,
#    "K2S2": doping_K2S2,
#    "K2S4": doping_K2S4,
#    "K2S6": doping_K2S6,
#    "K2S8": doping_K2S8,
#      "S8": doping_S8
#}
#
#for i in ["K2S", "K2S2", "K2S4", "K2S6", "K2S8", "S8"]:
##for i in ['']:
#    func = K_func_map[str(i)]
#
#    #filename='gCN_melon.vasp'
#    #slab = ase.io.read(filename)
#    #absorption_site = (slab.positions[7] + slab.positions[29]) / 2
#    #target_vector = slab.positions[7] - slab.positions[29]
#    #func(filename, absorption_site, target_vector)
#
#    #filename='gCN.vasp'
#    #slab = ase.io.read(filename)
#    #absorption_site = (slab.positions[2] + slab.positions[39]) / 2
#    #target_vector = slab.positions[2] - slab.positions[39]
#    #func(filename, absorption_site, target_vector)
#
#    filename='graphene_551.vasp'
#    slab = ase.io.read(filename)
#    absorption_site = slab.positions[25]
#    target_vector = slab.positions[34] - slab.positions[26]
#    func(filename, absorption_site, target_vector, h=2.6)
#    
#    filename='N3_graphene_551.vasp'
#    slab = ase.io.read(filename)
#    absorption_site = (slab.positions[46] + slab.positions[47]+ slab.positions[48]) / 3
#    target_vector = slab.positions[47] - slab.positions[31]
#    func(filename, absorption_site, target_vector, h=2.3)
#    
#    filename='N4_graphene_551.vasp'
#    slab = ase.io.read(filename)
#    absorption_site = (slab.positions[44] + slab.positions[45]+ slab.positions[46] + slab.positions[47]) / 4
#    target_vector = slab.positions[47] - slab.positions[44]
#    func(filename, absorption_site, target_vector, h=1.8)

