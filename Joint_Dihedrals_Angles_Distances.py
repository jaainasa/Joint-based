#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  10 23:01:08 2020

@author: jaainasa
"""
"Program to calculate bond and Theta angles, bond distances and Phi, Psi, Omega, & Tau dihedral angles from PDB structures"
"The program uses PDB structures (PDB format) as input and generates angles, distances and dihedrals output files for each chain of any given PDB structure"
import sys
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.vectors import calc_angle
from Bio.PDB.vectors import calc_dihedral

# Read PDB structure 

structure_name = "1a0t_f.B99990001"
parser=PDBParser(PERMISSIVE=True)
structure = parser.get_structure(structure_name, "%s.pdb" % structure_name)
standard_aa_residue_names = [
        "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
        "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"
] 

#For measuring the bonds angles N(i)-CA(i)-C(i); CA(i)-C(i)-N(i+1) and C(i)-N(i+1)-CA(i+1)
def get_angles(structure):
    result = []
    chain_result = []
    for model in structure:
        for chain in model:
            if chain.id != ' ':
                vector = []
                current_residues = []
                for residue in chain:
                    if residue.get_resname() in standard_aa_residue_names:
                        for atom in residue.get_atoms():
                                if atom.get_name() == 'N' or atom.get_name() == 'CA' or atom.get_name() == 'C':
                                    vector.append(atom.get_vector())
                                    current_residues.append(atom.get_parent())
                                    if len(vector) == 3:
                                        angle = (calc_angle(vector[0], vector[1], vector[2]))*180/np.pi
                                        current_chain = atom.get_parent().get_parent()
                                        current_residue = current_residues[0]
                                        chain_result.append((current_chain.id, current_residue.id[1], current_residue.get_resname(), round(angle,2)))
                                        del vector[0]
                                        del current_residues[0]
                if not len(chain_result) == 0:
                    chain_result.append((chain_result[-1][0], chain_result[-1][1], chain_result[-1][2], 'NA'))
                    chain_result.append((chain_result[-1][0], chain_result[-1][1], chain_result[-1][2], 'NA'))
                    result.append(chain_result)
                    chain_result = []
    return result

#For measuring the phi:C(i-1)-N(i)-CA(i)-C(i); psi:N(i)-CA(i)-C(i)-N(i+1) and omega:CA(i)-C(i)-N(i+1)-CA(i+1) dihedral angles
def get_dihedral_angles(structure):
    result = []
    chain_result = []
    for model in structure:
        for chain in model:
            if chain.id != ' ':
                vector = []
                current_residues = []
                is_new_chain = True
                for residue in chain:
                    if residue.get_resname() in standard_aa_residue_names:
                        for atom in residue.get_atoms():
                                if atom.get_name() == 'N' or atom.get_name() == 'CA' or atom.get_name() == 'C':
                                    vector.append(atom.get_vector())
                                    if is_new_chain == True:
                                        chain_result.append((chain.id, residue.id[1], residue.get_resname(), 'NA'))
                                        is_new_chain = False
                                    else:
                                        current_residues.append(atom.get_parent())
                                    if len(vector) == 4:
                                        angle = (calc_dihedral(vector[0], vector[1], vector[2], vector[3]))*180/np.pi
                                        current_chain = atom.get_parent().get_parent()
                                        current_residue = current_residues[0]
                                        chain_result.append((current_chain.id, current_residue.id[1], current_residue.get_resname(), round(angle,2)))
                                        del vector[0]
                                        del current_residues[0]
                if not len(chain_result) == 0:
                    chain_result.append((chain_result[-1][0], chain_result[-1][1], chain_result[-1][2], 'NA'))
                    chain_result.append((chain_result[-1][0], chain_result[-1][1], chain_result[-1][2], 'NA'))
                    result.append(chain_result)
                    chain_result = []
    return result

#For measuring the theta angles of all 3 consecutive Calpha atoms
def get_theta(structure):
    result = []
    chain_result = []
    for model in structure:
        for chain in model:
            if chain.id != ' ':
                i = 1
                vector = []
                for residue in chain:
                    if residue.get_resname() in standard_aa_residue_names:
                        for atom in residue.get_atoms():
                                if atom.get_name() == 'CA':
                                    vector.append(atom.get_vector())
                                    if len(vector) == 3:
                                        theta = (calc_angle(vector[0], vector[1], vector[2]))*180/np.pi
                                        chain_result.append((chain.id, i, round(theta,2)))
                                        i = i+1
                                        del vector[0]
                if not len(chain_result) == 0:
                    result.append(chain_result)
                    chain_result = []
    return result

#For measuring the tau dihedral angles of all 4 consecutive Calpha atoms
def get_tau(structure):
    result = []
    chain_result = []
    for model in structure:
        for chain in model:
            if chain.id != ' ':
                i = 1
                vector = []
                for residue in chain:
                    if residue.get_resname() in standard_aa_residue_names:
                        for atom in residue.get_atoms():
                                if atom.get_name() == 'CA':
                                    vector.append(atom.get_vector())
                                    if len(vector) == 4:
                                        theta = (calc_dihedral(vector[0], vector[1], vector[2], vector[3]))*180/np.pi
                                        chain_result.append((chain.id, i, round(theta,2)))
                                        i = i+1
                                        del vector[0]
                if not len(chain_result) == 0:
                    result.append(chain_result)
                    chain_result = []
    return result

#For measuring the bond distances between N(i)-CA(i); CA(i)-C(i) and C(i)-N(i+1) atoms
def get_distance(structure):
    result = []
    chain_result = []
    for model in structure:
        for chain in model:
            if chain.id != ' ':
                atoms = []
                current_residues = []
                for residue in chain:
                    if residue.get_resname() in standard_aa_residue_names:
                        for atom in residue.get_atoms():
                                if atom.get_name() == 'N' or atom.get_name() == 'CA' or atom.get_name() == 'C':
                                    atoms.append(atom)
                                    current_residues.append(atom.get_parent())
                                    if len(atoms) == 2:
                                        distance = atoms[0] - atoms[1]
                                        current_chain = atom.get_parent().get_parent()
                                        current_residue = current_residues[0]
                                        chain_result.append((current_chain.id, current_residue.id[1], current_residue.get_resname(), round(distance,2)))
                                        del atoms[0]
                                        del current_residues[0]
                if not len(chain_result) == 0:
                    chain_result.append((chain_result[-1][0], chain_result[-1][1], chain_result[-1][2], 'NA'))
                    result.append(chain_result)
                    chain_result = []
    return result
 
# Write angles
def print_angles(structure, structure_name):
    for angle in get_angles(structure):
        i = 0
        new_chain = True
        for angle_chain in angle:
            if new_chain == True:
                file_name = structure_name + "_" + angle_chain[0] + "_bond_angles.txt"
                sys.stdout = open(file_name,"w")
                print 'Chain' + '\t' + 'ResNum' + '\t' + 'ResName' +  '\t' + 'NCAC' + '\t' + 'CACN' + '\t' + 'CNCA'
                new_chain = False
            if i == 3:
                print '\n'
                i=0
            if i == 0:
                print str(angle_chain[0]) + '\t' + str(angle_chain[1]) + '\t' + str(angle_chain[2]) + '\t' + str(angle_chain[3]) + '\t',
                i=i+1
            else:
                print str(angle_chain[3]) + '\t',
                i=i+1

# Write dihedral_angles
def print_dihedral_angles(structure, structure_name):
    for angle in get_dihedral_angles(structure):
        i = 0
        new_chain = True
        for angle_chain in angle:
            if new_chain == True:
                file_name = structure_name + "_" + angle_chain[0] + "_dihedral_angles.txt"
                sys.stdout = open(file_name,"w")
                print 'Chain' + '\t' + 'ResNum' + '\t' + 'ResName' +  '\t' + 'Phi' + '\t' + 'Psi' + '\t' + 'Omega'
                new_chain = False
            if i == 3:
                print '\n'
                i=0
            if i == 0:
                print str(angle_chain[0]) + '\t' + str(angle_chain[1]) + '\t' + str(angle_chain[2]) + '\t' + str(angle_chain[3]) + '\t',
                i=i+1
            else:
                print str(angle_chain[3]) + '\t',
                i=i+1
            
# Write CA_dihedral_angle
def print_tau(structure, structure_name):
    for tau in get_tau(structure):
        new_chain = True
        for tau_chain in tau:
            if new_chain == True:
                file_name = structure_name + "_" + tau_chain[0] + "_tau_dihedral.txt"
                sys.stdout = open(file_name,"w")
                print 'Chain' + '\t' + 'Num' + '\t' + 'Tau'
                new_chain = False
            print str(tau_chain[0]) + '\t' + str(tau_chain[1]) + '\t' + str(tau_chain[2])

# Write CA_bond_angle   
def print_theta(structure, structure_name):
    for theta in get_theta(structure):
        new_chain = True
        for theta_chain in theta:
            if new_chain == True:
                file_name = structure_name + "_" + theta_chain[0] + "_theta_angle.txt"
                sys.stdout = open(file_name,"w")
                print 'Chain' + '\t' + 'Num' + '\t' + 'Theta'
                new_chain = False
            print str(theta_chain[0]) + '\t' + str(theta_chain[1]) + '\t' + str(theta_chain[2])
        
# Write_bond_distances        
def print_distance(structure,structure_name):
    for distance in get_distance(structure):
        new_chain = True
        i = 0
        for distance_chain in distance:
            if new_chain == True:
                file_name = structure_name + "_" + distance_chain[0] + "_bond_distances.txt"
                sys.stdout = open(file_name,"w")
                print 'Chain' + '\t' + 'ResNum' + '\t' + 'ResName'  + '\t' + 'D:N-CA' + '\t' + 'D:CA-C' + '\t' + 'D:C-N'
                new_chain = False
            if i == 3:
                print '\n'
                i = 0
            if i == 0:
                print str(distance_chain[0]) + '\t' + str(distance_chain[1]) + '\t' + str(distance_chain[2]) + '\t' + str(distance_chain[3]) + '\t',
                i=i+1
            else:
                print str(distance_chain[3]) + '\t',
                i=i+1

def run():
        print_angles(structure, structure_name)
        
        print_dihedral_angles(structure, structure_name)
                
        print_tau(structure, structure_name)
        
        print_theta(structure, structure_name)
               
        print_distance(structure, structure_name)

run()

