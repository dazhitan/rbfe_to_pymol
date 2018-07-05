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
import pandas as pd
import argparse
import parmed
from glob import glob

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                description=__doc__,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-m', '--mode')
    parser.add_argument('-e', '--ene_type', 
                        choices=['total', 'dispersion', 'exchange', 'electrostatic',
                                 'induction_AB', 'induction_BA'],
                        default='total')
    args = parser.parse_args()

    mean_csvs = glob('fsapt_%s/*/*mean*%s*csv' % (args.mode, args.ene_type))
    std_csvs = glob('fsapt_%s/*/*std*%s*csv' % (args.mode, args.ene_type))
    pymol_pdbs = glob('fsapt_%s/*/*%s*pymol*pdb' % (args.mode, args.ene_type))

    mean_csvs.sort()
    std_csvs.sort()

    post_df = pd.DataFrame(columns=['Mean_ene', 'Std_ene', 'Normalized_ene'])

    for mc in mean_csvs:
        ligname = os.path.basename(os.path.dirname(mc))
        ligmol = parmed.load_file(os.path.join('MD', ligname, 'cmp_sybyl.mol2'),
                                  structure=True)
        num_Cs = len([a for a in ligmol.atoms if a.atomic_number == 6])
        value = pd.read_csv(mc, index_col=0).loc['Total', 'Total']
        post_df.loc[ligname, 'Mean_ene'] = value
        post_df.loc[ligname, 'Normalized_ene'] = value / num_Cs

    for sc in std_csvs:
        ligname = os.path.basename(os.path.dirname(sc))
        value = pd.read_csv(sc, index_col=0).loc['Total', 'Total']
        post_df.loc[ligname, 'Std_ene'] = value

    post_df.sort_values(by='Normalized_ene', inplace=True)
    post_df.to_csv('sorted_%s_%s.csv' % (args.mode, args.ene_type))

    with open('%s_%s.pml' % (args.mode, args.ene_type), 'w') as pmlfh:
        for pp in pymol_pdbs:
            objname = os.path.basename(os.path.dirname(pp))
            print >> pmlfh, 'load %s, %s' % (pp, objname)

        print >> pmlfh, 'hide everything'
        print >> pmlfh, 'show sticks'
        print >> pmlfh, 'label resname L1+L2 and not hydrogens, "%s/%+.1f" % (name, q)'
        print >> pmlfh, 'label not resname L1+L2 and name ca, "%s%s/%+.1f" % (resn, resi, q)'
        print >> pmlfh, 'label not resname L1+L2 and name n, "%s-%s/%+.1f" % (resi-1, resi, q)'
        print >> pmlfh, 'label not resname L1+L2 and name c, "%s-%s/%+.1f" % (resi, resi+1, q)'
        print >> pmlfh, 'spectrum q, green_white_red, minimum=-10, maximum=3.0'
        print >> pmlfh, 'color white, q = 0.0'
