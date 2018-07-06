#!/usr/bin/env python
"""
@author: Dazhi Tan
@created: July 5th, 2018
"""
import sys
import os
import argparse
from glob import glob

import numpy as np
import pandas as pd

import mdtraj as mtj

def extract_frame(parm_file, traj_file, num_frames, outpdb):  

    traj = mtj.load(traj_file, top=parm_file)
    if num_frames == 1:
        sel_frames = [traj.n_frames - 1]
    elif num_frames == 2:
        sel_frames = [0, traj.n_frames - 1]
    elif num_frames > 2:
        sel_frames = [0, traj.n_frames - 1]
        rand_frames = np.random.randint(1, traj.n_frames - 2, num_frames - 2)
        sel_frames.extend(list(rand_frames))

    sel_string = "not (water or resname 'K+' or resname 'Na+' or resname 'Cl-')"
    com_traj = traj.atom_slice(traj.top.select(sel_string)).slice(sel_frames)
    com_traj.save_pdb(outpdb)

def read_footprint(footprint_file):

    df = pd.read_csv(footprint_file, sep='\s+')

    return df

def assign_footprint_values(inpdb, outpdb, data_df):

    lines = open(inpdb, 'r').readlines()
    with open(outpdb, 'w') as outfh:
        for l in lines:
            l = l.rstrip()
            if 'ATOM' in l:
                _, atomid, atomname, resname, chainid, resid, x, y, z, _, _, elem = l.split()

                data_idx = data_df.index[data_df.ResNum == int(resid)]
                assert resname == data_df.loc[data_idx, 'Residue'], \
                    raise ValueError("Please check the footprint file.")
                new_occu = data_df.loc[data_idx, 'Total-Energy']
                print >> outfh, 'ATOM  %5s%5s%4s%2s%4s%12.3f%8.3f%8.3f%6.1f  0.00%12s' \
                                % (atomid, atomname, resname, chainid, resid,
                                   float(x), float(y), float(z), float(new_occu), elem.upper())
            else:
                print >> outfh, l

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                description=__doc__,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--rbfe_dir')
    parser.add_argument('-l', '--ligname')
    parser.add_argument('-n', '--num_frames', type=int, default=1) 
    args = parser.parse_args()

    if args.n < 1:
        raise ValueError("--num_frames has to be a positive integer.")
