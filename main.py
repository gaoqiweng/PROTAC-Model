#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os
import argparse
import utils.frodock as fro
import utils.rosetta as ros


def parse_args():
    """Parses arguments from cmd"""
    parser = argparse.ArgumentParser(description="Modeling of PROTAC-mediated ternary complexes.")

    parser.add_argument('-irec',"--input-receptor-pdb",
                        help='PDB file of receptor protein should include the small molecular binder and '
                              'exclude other heteroatoms. Please submit the larger protein as the receptor.',
                        type=str, required=True, metavar='<string>')
    parser.add_argument('-ilig',"--input-target-pdb",
                        help='PDB file of target protein should include the small molecular binder and '
                              'exclude other heteroatoms. Please submit the smaller protein as the target.',
                        type=str, required=True, metavar='<string>')
    parser.add_argument('-site', '--input-docking-site',
                        help='The docking site in the receptor for protein-protein docking. '
                             'e.g. -site=-35.73,13.75,-27.92',
                        type=str, required=True, metavar='X,Y,Z')
    parser.add_argument('-ismi', '--input-protac-smi', help='Smiles file of PROTAC.',
                        type=str, required=True, metavar='<string>')
    parser.add_argument('-o', '--output-filepath', help='The filepath for storing files.',
                        type=str, required=True, metavar='<string>')
    parser.add_argument('-cpu', help='Number of cpu to calculate. Default value: 1',
                        type=int, default='1', metavar='<int>')
    parser.add_argument('-ie3lig1', '--input-e3-ligand-sdf1',
                        help='First sdf file of the e3 ligand which contains two possible locations '
                             'for the attachment of the atoms, e.g. thalidomide.',
                        type=str, default='none', metavar='<string>')
    parser.add_argument('-ie3lig2', '--input-e3-ligand-sdf2',
                        help='Second sdf file of the e3 ligand which contains two possible locations '
                             'for the attachment of the atoms, e.g. thalidomide.',
                        type=str, default='none', metavar='<string>')
    parser.add_argument('-refine', '--rosettadock-refinement',
                        help='Use the RosettaDock-based refinement to improve the prediction performance. '
                             'But it will cost much more time.', action="store_true", default=False)

    return parser.parse_args()

def main():
    """Main function"""
    args = parse_args()

    filepath_rec = args.input_receptor_pdb
    filepath_lig = args.input_target_pdb
    site = args.input_docking_site
    filepath_smiles = args.input_protac_smi
    filepath_out = args.output_filepath
    cpu = args.cpu
    filepath_e3sdf_1 = args.input_e3_ligand_sdf1
    filepath_e3sdf_2 = args.input_e3_ligand_sdf2
    refine = args.rosettadock_refinement
    lig_locate_num = 1
    filepath_frodock = filepath_out + '/frodock'
    if os.path.exists(filepath_out) == False:
        os.makedirs(filepath_out)
    if os.path.exists(filepath_frodock) == False:
        os.makedirs(filepath_frodock)
    if filepath_e3sdf_1 != 'none' and filepath_e3sdf_2 != 'none':
        os.system('cp %s %s/rec_lig_1.sdf' % (filepath_e3sdf_1, filepath_frodock))
        os.system('cp %s %s/rec_lig_2.sdf' % (filepath_e3sdf_2, filepath_frodock))
        lig_locate_num = 2
    os.system('cp %s %s/receptor.pdb' % (filepath_rec, filepath_frodock))
    os.system('cp %s %s/target.pdb' % (filepath_lig, filepath_frodock))
    os.system('cp %s %s/protac.smi' % (filepath_smiles, filepath_frodock))
    os.chdir(filepath_frodock) #working directory is in filepath_frodock

    #FRODOCK docking
    fro.frodock(site)
    fro.filter_frodock(cpu, lig_locate_num)

    #Rosetta docking
    filepath_rosetta = 'rosetta'
    if refine:
        if os.path.exists('../%s' % filepath_rosetta) == False:
            os.makedirs('../%s' % filepath_rosetta)
        os.chdir('../%s' % filepath_rosetta)
        ros.rosetta(cpu, lig_locate_num)

if __name__ == '__main__':
    main()
