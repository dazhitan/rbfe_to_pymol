#!/software/anaconda2/bin/python
"""
@author: Dazhi Tan
@created: March 15th, 2018
"""
import sys
import os
import argparse
import shutil
import itertools
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
