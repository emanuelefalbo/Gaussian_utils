#!/usr/bin/env python3

"""
Example of Json
{
    "trjprm": {"sphradius": 30,
               "solres": 1},
    "files": {"topdata": "frame0_mod.pdb",
              "trjdata": "weighted_frames_100_349.pdb",
              "solute": "solute.dat",
              "solvent": "solvent.dat"},
    "qmdata": {"functional": "PW6B95D3",
               "basis": "gen 6D 10F",
               "molchr": 0,
               "molspn": 2,
               "addroot": "TD=(nstates=11)"
               "addline": "@/home/m.fuse/basis/SNSD.gbs"}
}
Example of the solvent dat file
O -0.848440 15.9949146
H  0.424220 1.0078250
H  0.424220 1.0078250
"""


import os
import sys
import argparse
import copy
import json
import warnings

import numpy as np
from numpy.linalg import norm

import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
import MDAnalysis.analysis.distances as ds

DEBUG=False

## Generic functions
def read_optfile(fname):
    """
    Read the option file. expected in json style
    """
    if not os.path.exists(fname):
        raise FileNotFoundError('Missing JSON file')
    with open(fname, 'r', encoding='UTF-8') as fopen:
        data = json.load(fopen)
    for i in ['topdata', 'trjdata', 'solute', 'solvent']:
        if not os.path.exists(data['files'][i]):
            raise FileNotFoundError(f'Missing {i} file')
    return data

def read_molinfo(fname):
    """
    Small function to parse the molecule dat file
    """
    data = {'atm':[], 'charge': [], 'mass':[]}
    with open(fname, 'r', encoding='UTF-8') as fopen:
        for line in fopen:
            tokens = line.split()
            if len(tokens) > 1:
                data['atm'].append(tokens[0])
                data['charge'].append(float(tokens[1]))
                data['mass'].append(float(tokens[2]))
    return data

##Functions
def get_close_solvent(traj, frame, atom1, atom2type, nsolv, guessradius=4):
    """Returns the indices of the N(nsolv) residues closest to the selected atom
    Args:
        traj (mda.Universe): mdanalysis trajectory
        frame (int): The frame number
        atom1 (int): atomic index of the selected atom
        atom2type (str): atom type of the atom to search (e.g. HO)
        nsolv (int): number of solvent molecules to return
        guessradius (int, optional): [description]. Defaults to 4.
    Returns:
        [list(int)]: residue indices
    """
    traj.universe.trajectory[frame]
    atm_pos = traj.select_atoms('bynum {}'.format(atom1))
    nsolute = 0
    srad = guessradius
    while not nsolute >= nsolv:
        whs = traj.select_atoms('name {:s} and around {:f} (bynum {:d})'.format(atom2type, srad, int(atom1)))
        nsolute = len(list(set(whs.ids.tolist())))
        #nsolute = whs.ids.shape[0]
        srad += 0.5

    distances = ds.distance_array(atm_pos.positions, whs.atoms.positions).flatten()
    resind = whs.resids
    # Checks for duplicates
    dist_2 = []
    resind_2 = []
    for i, val in enumerate(resind):
        if val not in resind_2:
            dist_2.append(distances[i])
            resind_2.append(val)
    sortresid = np.array([x for _,x in sorted(zip(dist_2,resind_2))])[:nsolv]

    return sortresid

def get_drop_chrg(traj, drop=30, solute_res=1):
    tmp_sel = traj.select_atoms("sphzone {drop:.3f} (resid {solute:d})".format(drop=drop, solute=solute_res))
    solv_indx = list(set(tmp_sel.resids))
    return solv_indx

def write_charge_frame(traj, frame, qmindex, solchrs, 
                       qmdata, drop=30,
                       out_file='sel_frame', where=""):
    """Write a gaussian input file
    Args:
        traj ([type]): [description]
        frame ([type]): [description]
        qmindex ([type]): [description]
        solchrs ([type]): [description]
        out_file (str, optional): [description]. Defaults to 'sel_frame'.
    """
    template ="""%mem=40GB
%nprocshared=16
%chk={CHK}
#p {FUN}/{BASIS} Charge
 nosymm {ADDROOT}

CQ cumulative structure Charges

{MOLCHR} {MOLSPN}
"""
    # select the frame
    traj.universe.trajectory[frame]
    geomstr = ""
    chrgstr = ""
    line_a = '{a:4s}{b[0]:12.6f}{b[1]:12.6f}{b[2]:12.6f}\n'
    line_c = '{a[0]:12.6f}{a[1]:12.6f}{a[2]:12.6f}{b:12.8f}\n'
    solv_drop = get_drop_chrg(traj, drop=drop)
    for i in qmindex:
        solv_drop.remove(i)
    for res in traj.residues:
        sel = traj.select_atoms('resid {} and (not name VS*)'.format(res.resid))
        for i, pos in enumerate(sel.atoms.positions):
            if res.resid in qmindex:
                geomstr += line_a.format(a=sel.elements[i], b=pos)
            elif solv_drop:
                chrgstr += line_c.format(a=pos, b=solchrs[i])
    geomstr += '\n'
    chrgstr += '\n'
    # write the molecule
    fout = out_file+'_frame{:03d}_charge.gjf'.format(frame)
    with open(where+fout, 'w') as fopen:
        fopen.write(template.format(CHK=fout[:-3]+'chk',
                                    FUN=qmdata['functional'],
                                    BASIS=qmdata['basis'],
                                    MOLCHR=qmdata['molchr'],
                                    MOLSPN=qmdata['molspn'],
                                    ADDROOT=qmdata['addroot']))
        fopen.write(geomstr)
        fopen.write(chrgstr)
        fopen.write(f'{qmdata["addline"]}\n\n')
        #fopen.write('@/home/m.fuse/basis/SNSD.gbs\n\n')

def build_parser():
    """Builds options parser.
    Builds the full option parser.
    Returns
    -------
    :obj:`ArgumentParser`
        `ArgumentParser` object
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('optfile', nargs='?', help='option file in json')
    parser.add_argument('patom', default=1, type=int,
                        help='Pivot atom on the solute around which look for solvent molecules')
    parser.add_argument('solattype', default='OW', type=str,
                        help='atom type of the solvent to search (es. "OW")')
    parser.add_argument('-ns', '--nsolvent', default=0, type=int,
                        help='number of solvent molecule to be included in QM')
    txt = """Reference frame type:
        single: only one frame used usually the cluster centroid
        all: all the frame in the trajectory
        """
    parser.add_argument('-m', '--mode', choices=['single', 'all'],
                        default='single', help=txt)
    parser.add_argument('--frame', default=-1, type=int, help='the frame to be\
                        used as reference in single mode')
    parser.add_argument('-p', '--prefix', default='fromtraj',
                        type=str, help='prefix name for the output files')
    parser.add_argument('--silent', action='store_true', help="Suppres almost all the printing")
    return parser


def main():
    """
    The main
    """
    parser = build_parser()
    opts = parser.parse_args()
    jopts = read_optfile(opts.optfile)
    # Name of the topology and of the trajectory
    top_name = jopts['files']['topdata']
    trj_name = jopts['files']['trjdata']
    # load it with mdanalysis
    utraj = mda.Universe(top_name , trj_name)
    # User defined files containing for the solute and the solvent
    # AtomName Charge Mass
    solvent = read_molinfo(jopts['files']['solvent'])

    # Variabiles to inquire md analysis
    residue = jopts['trjprm']['solres']
    sphradius = jopts['trjprm']['sphradius']

    atom = int(opts.patom)
    atomtyp = opts.solattype
    sframe = int(opts.frame)
    nsolv = int(opts.nsolvent)
    nframe = utraj.universe.trajectory.n_frames
    sel = utraj.select_atoms('resid {} and (not name VS*)'.format(residue))
    natoms = sel.n_atoms
    if atom >= natoms:
        print("Atom not in the solute")
        sys.exit()
    if -1 < sframe or sframe >= nframe:
        print("Frame not in the traj")
        sys.exit()
    if opts.mode == "all":
        framestd = list(range(nframe))
    else:
        framestd = [sframe]

    for frm in framestd:
        quant = [residue]
        if nsolv > 0:
            quant += get_close_solvent(utraj, frm,
                                       atom, atomtyp,
                                       nsolv, guessradius=4).tolist()
        # writes the ref in gm traj the centroid is the last
        write_charge_frame(utraj, frm, quant,
                           solvent['charge'], jopts["qmdata"],
                           drop=30,
                           out_file='_nframe{:d}_nsolv{:02d}'.format(frm, nsolv))



if __name__ == '__main__':
   if not DEBUG:
        # Hide MDAnalysis warnings about the virtualsites
        warnings.filterwarnings("ignore")
        main()

