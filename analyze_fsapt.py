#!/software/anaconda2/bin/python
"""
@author: Dazhi Tan
@created: March 15th, 2018
"""
import sys
import os
import re
sys.path.insert(0, '/software/anaconda2/lib/python2.7/site-packages')
pkg_dir = '/software/ffgen/ffgen_v3'
sys.path.insert(0, os.path.dirname(pkg_dir))
import argparse
import cPickle
import shutil
import tempfile
from glob import glob

import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.rdFMCS as rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry.rdGeometry import Point3D

import mdtraj as mtj
import ffgen_v3.chemsystems as cs
import ffgen_v3.geometry as geometry
import ffgen_v3.const as const

import seaborn
import svgutils.compose as sc

import getpass

def fragment_mol(mol, frag_label):
    excl_smarts = ['[CX3](=[OX1])[OX1]',
                   '[SX4](=[OX1])=[OX1]',
                   '[NX3](=[OX1])[OX1]',
                   '[OX1]=C[NX3][*]',
                   '[#6r5]~[#7r6]~[#6r6]~[#7r6]~[#6r6]~[#6r5]', # Special for aspire
                   '[SX3](=[OX1])']
    excl_substructs = [Chem.MolFromSmarts(smart) for smart in excl_smarts]
    excl_idx = []
    for es in excl_substructs:
        excl_idx.extend(list(mol.GetSubstructMatches(es)))
    excl_idx = [i for sublist in excl_idx for i in sublist]

    sigma_bond_idx = []
    frag_dict = dict()
    for bnd in mol.GetBonds():
        if bnd.GetIsAromatic():
            continue
        elif bnd.GetBeginAtom().GetAtomicNum() == 1:
            continue
        elif bnd.GetEndAtom().GetAtomicNum() == 1:
            continue
        elif bnd.GetEndAtomIdx() in excl_idx and bnd.GetBeginAtomIdx() in excl_idx:
            continue
        else:
            sigma_bond_idx.append(bnd.GetIdx())

    fragmentation = Chem.FragmentOnBonds(mol, sigma_bond_idx, addDummies=False)
    frags_as_idx = Chem.GetMolFrags(fragmentation)

    for n, indices in enumerate(frags_as_idx):
        resname = frag_label
        resid = n + 1
        for idx in indices:
            frag_dict[idx] = {'resname':resname, 'resid':resid}

    return frag_dict, frags_as_idx
def _get_ene_matrix(result_file, ene_type, ignore_chg=True):
    
    elabels = {'total': 'EsaptAB',
               'dispersion': 'EdispAB',
               'electrostatic': 'EelstAB',
               'induction_AB': 'EindAB_AB',
               'induction_BA': 'EindBA_AB',
               'exchange': 'EexchAB'}

    try:
        blocks = open(result_file, 'r').read().split(elabels[ene_type])[-1].split('\n\n\n')[0]
        blocks = blocks.split('\n\n')[1:]
        data_dfs = []
        for b in blocks:
            tmpdir = tempfile.mkdtemp(dir='/dev/shm')
            csv_path = os.path.join(tmpdir, 'tmp.csv')
            with open(csv_path, 'w') as fh:
                print >> fh, b
            df = pd.read_csv(csv_path, index_col=0, sep='\s+')
            dummy_cols = [c for c in df.columns if 'dummy' in c]
            for dc in dummy_cols:
                df.drop(dc, axis=1, inplace=True)
            chg_cols = [c for c in df.columns if re.match('ASP|GLU|LYS|CYM|ARG', c)]
            for cc in chg_cols:
                df[cc] = np.nan
            data_dfs.append(df)
            os.system('rm -f %s' % csv_path)
        matrix_df = pd.concat(data_dfs, axis=1)
        matrix_df.loc[:, 'Total'] = matrix_df.sum(axis=1)
        matrix_df.loc['Total', :] = matrix_df.sum(axis=0)
        return matrix_df
    except:
        return None

def fsapt_analyze(lig_dir, mode, ene_type):
    lig_name = os.path.basename(os.path.abspath(lig_dir))
    matrix_dfs = []
    outfiles = glob('%s/FSAPT*out' % lig_dir) 
    for of in outfiles:
        df = _get_ene_matrix(of, ene_type)
        if not df is None:
            matrix_dfs.append(df)

    all_df = pd.concat(matrix_dfs, axis=1)
    mean_df = all_df.stack().groupby(level=[0,1]).mean().unstack()
    std_df = all_df.stack().groupby(level=[0,1]).std().unstack()

    if mode in ['prolig', 'proliglig']:
        old_columns = mean_df.columns[:]
        new_labels = []
        numbering = []
        for old_label in old_columns:
            if old_label == 'Total':
                new_labels.append('Total')
                numbering.append(100000)
            else:
                labels = old_label.split('-')
                if len(labels) == 2:
                    new_labels.append(''.join(labels))
                    numbering.append(float(labels[-1]))
                elif len(labels) == 3:
                    new_labels.append('-'.join(labels[1:]))
                    numbering.append(0.5*(float(labels[-1]) + float(labels[-2])))
        new_columns = [nl for _, nl in sorted(zip(numbering, new_labels))]
        old_columns = [ol for _, ol in sorted(zip(numbering, old_columns))]

        new_mean_df = pd.DataFrame()
        new_std_df = pd.DataFrame()
        for nc, oc in zip(new_columns, old_columns):
            new_mean_df[nc] = mean_df[oc]
            new_std_df[nc] = std_df[oc]
        mean_df = new_mean_df 
        std_df = new_std_df

    mean_anno = mean_df.applymap(lambda x: '%+.2f\n' % x)
    std_anno = std_df.applymap(lambda x: r'+/-%.2f' % x)
    all_anno = mean_anno + std_anno

    matrix_svg = '%s/ene_matrix_%s.svg' % (lig_dir, ene_type)
    plot_matrix(mean_df, all_anno, matrix_svg, mode, ene_type)

    mean_df.to_csv('%s/ene_mean_%s_%s_%s.csv' % (lig_dir, lig_name, mode, ene_type))
    std_df.to_csv('%s/ene_std_%s_%s_%s.csv' % (lig_dir, lig_name, mode, ene_type))

    # Plot the ligand
    dpi = 96
    width = len(mean_df.columns) + 2
    height = 4

    ligmol = cs._RdkitMolBase.from_file('MD/%s/cmp_sybyl.mol2' % lig_name)
    ligmol._init_atominfo(reset=False)
    ligmol.charged_mol2file = 'MD/%s/cmp_sybyl.mol2' % lig_name
    ligmol.get_noh_mol()
    AllChem.Compute2DCoords(ligmol.noh_mol, canonOrient=True, bondLength=1.5)
    drawer =rdMolDraw2D.MolDraw2DSVG(width*dpi, height*dpi)
    opts = drawer.drawOptions()
    opts.additionalAtomLabelPadding = 0.1

    frag_dict, _ = fragment_mol(ligmol, 'L1')

    for noha in ligmol.noh_mol.GetAtoms():
        noh_idx = noha.GetIdx()
        h_idx = ligmol.noh_to_h_atom_mapping[noh_idx]
        frag_label = str(frag_dict[h_idx]['resid'])
        if not 'L1-%02d' % int(frag_label) in mean_df.index:
            continue
        if noha.GetAtomicNum() == 6:
            opts.atomLabels[noh_idx] = '%02d' % int(frag_label)
        else:
            elem = ligmol.GetAtomWithIdx(h_idx).GetProp('_TriposAtomType').split('.')[0]
            opts.atomLabels[noh_idx] = '%s/%02d' % (elem, int(frag_label))
    drawer.DrawMolecule(ligmol.noh_mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    struct_svg = '%s/lig_frag_%s.svg' % (lig_dir, ene_type)
    with open(struct_svg, 'w') as fh: 
        fh.writelines(svg)

    # Consolidate the panels
    if mode == 'prolig':
        mat_title = 'Protein-Ligand %s Interaction' % ene_type.capitalize()
    else:
        mat_title = 'Ligand-Ligand %s Interaction' % ene_type.capitalize()

    mat_title = sc.Panel(sc.Text(mat_title, size=24)).move(20, 20)
    mat_panel = sc.Panel(sc.SVG(matrix_svg).scale(1.4)).move(0, 20)
    struct_title = sc.Panel(sc.Text('Ligand %s' % lig_name, size=24)).move(20, dpi*len(mean_df)+20)
    struct_panel = sc.Panel(sc.SVG(struct_svg)).move(0, dpi*len(mean_df)+20)
    final_figure = sc.Figure(dpi*width, dpi*(len(mean_df)+height)+40,
                             mat_panel, mat_title, struct_panel, struct_title)
    final_name = '%s/%s_%s_%s' % (lig_dir, lig_name, mode, ene_type)
    final_figure.save('%s.svg' % final_name)
    os.system('convert -density 100 %s.svg %s.pdf' % (final_name, final_name))
    os.system('rm -f %s %s' % (matrix_svg, struct_svg))

    # Write pdb for pymol
    inpdb = '%s/frame0/fsapt.pdb' % lig_dir
    outpdb = '%s_pymol.pdb' % final_name
    write_pymol_pdb(inpdb, outpdb, mean_df)

def write_pymol_pdb(inpdb, outpdb, data_df):

    lines = open(inpdb, 'r').readlines()
    with open(outpdb, 'w') as outfh:
        for l in lines:
            l = l.rstrip()
            if 'ATOM' in l:
                _, atomid, atomname, resname, chainid, resid, x, y, z, _, _, elem = l.split()
                if resname == 'L1':
                    row_label = '%s-%02d' % (resname, int(resid))
                    try:
                        new_occu = data_df.loc[row_label, 'Total']
                    except:
                        new_occu = 0.0
                elif resname == 'L2':
                    col_label = '%s-%02d' % (resname, int(resid))
                    try:
                        new_occu = data_df.loc['Total', col_label]
                    except:
                        new_occu = 0.0
                elif atomname in ['N', 'H']:
                    col_label = '%d-%d' % (int(resid)-1, int(resid))
                    try:
                        new_occu = data_df.loc['Total', col_label]
                    except:
                        new_occu = 0.0
                elif atomname in ['C', 'O']:
                    col_label = '%d-%d' % (int(resid), int(resid)+1)
                    try:
                        new_occu = data_df.loc['Total', col_label]
                    except:
                        new_occu = 0.0
                else:
                    col_label = '%s%d' % (resname, int(resid))
                    try:
                        new_occu = data_df.loc['Total', col_label]
                    except:
                        new_occu = 0.0
                print >> outfh, 'ATOM  %5s%5s%4s%2s%4s%12.3f%8.3f%8.3f%6.1f  0.00%12s' \
                                % (atomid, atomname, resname, chainid, resid,
                                   float(x), float(y), float(z), float(new_occu), elem.upper())
            else:
                print >> outfh, l

def plot_matrix(mean_df, anno_df, outfile, mode, ene_type):
    colors = [(0, 0.4, 0), (1, 1, 1), (0.8, 0, 0)] # 'G-W-R'
    n_bins = 100
    cmap_name = 'GWR'
    cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

    width = len(mean_df.columns) + 2
    height = len(mean_df)
    #fontsize = 1.0*np.mean([width, height])
    fontsize = 13
    
    fig = plt.figure(figsize=(width, height))
    ax = seaborn.heatmap(mean_df, vmin=-6, vmax=6, cmap=cm, 
                         annot=anno_df, fmt='', linewidths=1.0, linecolor='k',
                         annot_kws={'fontsize': fontsize},
                         )
    ax.tick_params(axis='both', labelsize=4+np.mean([width, height]))
    ax.tick_params(axis='x', labelbottom=False, labeltop=True)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=fontsize)
    plt.setp(ax.get_xticklabels(), rotation=20, fontsize=fontsize,
             ha='left', rotation_mode='anchor')
    cbar = fig.axes[-1]
    cbar.set_yticks(np.arange(-6, 7, 2))
    cbar.tick_params(labelsize=fontsize)
    cbar.set_ylabel(r'E_interaction (%s) (kcal/mol)' % ene_type.capitalize(), fontsize=fontsize)
    plt.savefig(outfile, dpi=96)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                description=__doc__,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dir')
    parser.add_argument('-m', '--mode')
    parser.add_argument('-e', '--ene_type', 
                        choices=['total', 'dispersion', 'exchange', 'electrostatic',
                                 'induction_AB', 'induction_BA'],
                        default='total')
    args = parser.parse_args()
    fsapt_analyze(args.dir, args.mode, args.ene_type)
