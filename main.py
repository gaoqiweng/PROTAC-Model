#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os
import argparse
import utils.frodock as fro
import utils.rosetta as ros
import utils.preprocess as pre

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
    parser.add_argument('-itsmi', '--input-target-smi', help='Smiles file of target ligand in target protein.',
                        type=str, default='none', metavar='<string>')
    parser.add_argument('-irsmi', '--input-receptor-smi', help='Smiles file of receptor ligand in receptor protein.',
                        type=str, default='none', metavar='<string>')

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
    filepath_target_smi = args.input_target_smi
    filepath_rec_smi = args.input_receptor_smi
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
    target_smi='none'
    if filepath_target_smi != 'none':
        with open(filepath_target_smi, 'r') as file_target_smi:
            target_smi = file_target_smi.read().splitlines()[0]
    rec_smi='none'
    if filepath_rec_smi != 'none':
        with open(filepath_rec_smi, 'r') as file_rec_smi:
            rec_smi = file_rec_smi.read().splitlines()[0]
    #alter the chain id of proteins to avoid the same chain id between receptor and target for rosetta docking
    pre.alter_pro_chain(filepath_rec, filepath_lig,
                        '%s/receptor.pdb' % filepath_frodock, '%s/target.pdb' % filepath_frodock)
    os.system('cp %s %s/protac.smi' % (filepath_smiles, filepath_frodock))
    os.chdir(filepath_frodock) #working directory is in filepath_frodock

    #FRODOCK docking
    fro.frodock(site)
    fro.filter_frodock(cpu, lig_locate_num, target_smi, rec_smi)

    #Rosetta docking
    filepath_rosetta = 'rosetta'
    if refine:
        if os.path.exists('../%s' % filepath_rosetta) == False:
            os.makedirs('../%s' % filepath_rosetta)
        os.chdir('../%s' % filepath_rosetta)
        ros.rosetta(cpu, lig_locate_num, target_smi, rec_smi)

if __name__ == '__main__':
    main()
