#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os,re
from string import digits
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

FRODOCK = os.environ["FRODOCK"]
ADFRSUITE = os.environ['ADFRSUITE']
VINA = os.environ['VINA']
VOROMQA = os.environ['VOROMQA']
FCC= os.environ['FCC']

#Obenergy and vina
def obenergy_vina(pdb_num, filepath_vina, filepath_out):
    #obenergy
    os.system(ADFRSUITE + '/bin/obenergy -h -ff GAFF %s/protac_%s.mol2 | grep "TOTAL ENERGY" '
                   '| awk \'{print $4}\' > %s/obenergy_%s' % (filepath_out, pdb_num, filepath_vina, pdb_num))
    os.system('nl %s/obenergy_%s | awk \'{if($2<10000) print $0}\' > %s/obenergy_process_%s'
              % (filepath_vina, pdb_num, filepath_vina, pdb_num))
    split_mol2_obenergy('%s/protac_%s.mol2' % (filepath_out, pdb_num), '.', '%s/obenergy_process_%s' % (filepath_vina, pdb_num))
    #vina
    os.system(ADFRSUITE + '/bin/prepare_receptor -r model_nolig.%s.pdb -o %s/model_nolig.%s.pdbqt '
                          '-A checkhydrogens -U nphs_lps_waters' % (pdb_num, filepath_out, pdb_num))
    #filter score
    conetnt_score = ''
    for i in range(1,101):
        if os.path.exists('protac_%s_%s.mol2' % (pdb_num, i)):
            #os.path.basename is used for "-l" and hence we can't specify the path of protac_%s_%s.mol2
            os.system(ADFRSUITE + '/bin/prepare_ligand -l protac_%s_%s.mol2 -o %s/protac_%s_%s.pdbqt -A checkhydrogens'
                      % (pdb_num, i, filepath_out, pdb_num, i))
            with os.popen(VINA+ '/bin/vina --score_only --receptor %s/model_nolig.%s.pdbqt '
                                '--ligand %s/protac_%s_%s.pdbqt | grep Affinity | cut -d\" \" -f2'
                          % (filepath_out, pdb_num, filepath_out, pdb_num, i), 'rb') as file_score:
                score = file_score.read().splitlines()[0]
            conetnt_score += '%s_%s %s\n' % (pdb_num, i, score)
    num = 0
    if conetnt_score != '':
        with open('%s/score_%s' % (filepath_vina, pdb_num), 'wb') as score_out:
            score_out.write(conetnt_score)
        os.system('paste %s/score_%s %s/obenergy_process_%s > %s/obenergy_merge_%s'
                  % (filepath_vina, pdb_num, filepath_vina, pdb_num, filepath_vina, pdb_num))
        os.system('awk \'{if($2<0 && $4<10000) print $1" "$2" "$4}\' %s/obenergy_merge_%s > %s/obenergy_filter_%s'
                  % (filepath_vina, pdb_num, filepath_vina, pdb_num))
        with os.popen('awk \'{if($2<0) print $0}\' %s/obenergy_filter_%s | wc -l'
                      % (filepath_vina, pdb_num), 'rb') as file_num:
            num = file_num.read().splitlines()[0]
        if int(num) > 0:
            os.system('echo \"%s %s\" >> %s/score_filter' % (pdb_num, num, filepath_vina))
            os.system('sort -k 2 -n  %s/obenergy_filter_%s | head -n 1 >> %s/score_all_top1' % (filepath_vina, pdb_num, filepath_vina))
    if filepath_out != '.':
        os.system('mv protac_%s_*.mol2 %s' % (pdb_num, filepath_out))
    return num
    #os.system('rm model_nolig.%s.* protac_%s.mol2 protac_%s_* %s/score_%s %s/obenergy_process_%s %s/obenergy_%s'
    #          % (pdb_num, pdb_num, pdb_num, filepath_vina, pdb_num, filepath_vina, pdb_num, filepath_vina, pdb_num))


#preprocess pdb files for clustering
def preprocess_cluster(file_input_pdb, file_out_pdb):
    with open(file_input_pdb, 'rb') as pdb_file:
        pdb_lines = pdb_file.read().splitlines()
    coord_re = re.compile('^(ATOM|HETATM)')
    space = 18 * ' '
    content = ''
    for pdb_line in pdb_lines:
        if coord_re.match(pdb_line):
            content += pdb_line[:54] + space + pdb_line[21].ljust(4) + '\n'
        else:
            content += pdb_line + '\n'
    with open(file_out_pdb, 'wb') as file_out:
        file_out.write(content)

#alter the chain id
def alter_chain(file_input_pdb, file_out_pdb, chain_id):
    with open(file_input_pdb, 'rb') as file_in:
        pdb_lines = file_in.read().splitlines()
    coord_re = re.compile('^(ATOM|HETATM)')
    content = ''
    for line in pdb_lines:
        line = line.strip()
        if coord_re.match(line):
            content += line[:21] + chain_id + line[22:] + '\n'
        else:
            content += line + '\n'
    with open(file_out_pdb, 'wb') as file_out:
        file_out.write(content)

#alter the chain id of protein in order
def alter_pro_chain(rec_input, target_input, rec_out, target_out):
    chain_list = []
    for i in range(65, 88): #exclude X,Y which is the id for small molecules
        chain_list.append('%c' % i)
    for i in range(97,123):
        chain_list.append('%c' % i)
    with open(rec_input, 'rb') as file_in:
        rec_lines = file_in.read().splitlines()
    with open(target_input, 'rb') as file_in:
        target_lines = file_in.read().splitlines()
    j = 0
    flag_model = 0
    first_chain_id = ''
    content_rec = ''
    content_target = ''
    for rec_line in rec_lines:
        if rec_line[:5] == 'MODEL':  # Only keep the first model
            if flag_model == 0:
                flag_model = 1
            else:
                break
        if rec_line[0:4] == 'ATOM':
            chain_id = rec_line[21]
            if chain_id != first_chain_id:
                first_chain_id = chain_id
                new_chain_id = chain_list[j]
                j += 1
            content_rec += rec_line[:21] + new_chain_id + rec_line[22:] + '\n'
        else:
            content_rec += rec_line + '\n'
    first_chain_id = ''
    flag_model = 0
    for target_line in target_lines:
        if target_line[:5] == 'MODEL':  # Only keep the first model
            if flag_model == 0:
                flag_model = 1
            else:
                break
        if target_line[0:4] == 'ATOM':
            chain_id = target_line[21]
            if chain_id != first_chain_id:
                first_chain_id = chain_id
                new_chain_id = chain_list[j]
                j += 1
            content_target += target_line[:21] + new_chain_id + target_line[22:] + '\n'
        else:
            content_target += target_line + '\n'
    with open(rec_out, 'wb') as file_out:
        file_out.write(content_rec)
    with open(target_out, 'wb') as file_out:
        file_out.write(content_target)

#get the chain id
def get_chain_id(file_input_pdb):
    with open(file_input_pdb, 'rb') as file_in:
        pdb_lines = file_in.read().splitlines()
    id_list = []
    for pdb_line in pdb_lines:
        if pdb_line[:4] == 'ATOM':
            chain_id = pdb_line[21]
            if chain_id not in id_list:
                id_list.append(chain_id)
    return id_list

#Cluster
def cluster(pdb_model_list, filepath_cluster, cpu):
    os.chdir(filepath_cluster)
    os.system('python %s/scripts/make_contacts.py -f %s -n %s -e %s/src/contact_fcc' % (FCC, pdb_model_list, cpu, FCC))
    os.system('sed -e \'s/pdb/contacts/\' pdb_model.list | sed -e \'/^$/d\' > pdb_model.contacts')
    os.system('python %s/scripts/calc_fcc_matrix.py -f pdb_model.contacts -o fcc_model_matrix.out' % FCC)
    os.system('python %s/scripts/cluster_fcc.py fcc_model_matrix.out 0.5 -o clusters_model_0.5_3.out -c 3' % FCC)
    os.system('python %s/scripts/ppretty_clusters.py clusters_model_0.5_3.out pdb_model.list > cluster_5_3_model' % FCC)
    get_best_cluster('cluster_5_3_model', 'cluster_5_3_results_top1')
    os.system(' sed \'s/model.//g\' cluster_5_3_model | sed \'s/_chainxseg//g\' > cluster_5_3_model_process')
    os.chdir('..')

#get the best one from each cluster
def get_best_cluster(file_input,file_output):
    re_cluster = re.compile('Cluster (.*?) -> (.*?)\n', re.S)
    re_number = re.compile('model.(\d*)')
    cluster_best_list = []
    all_number_list = []
    with open(file_input, 'rb') as cluster_file:
        cluster_lines = cluster_file.readlines()
    for cluster_line in cluster_lines:
        re_line = re.search(re_cluster, cluster_line)
        model = re_line.group(2)
        model_list = model.split()
        # get the number of model
        for i in range(len(model_list)):
            model_list[i] = int(re.search(re_number, model_list[i]).group(1))
            all_number_list.append(model_list[i])
        # get the bets model from cluster
        model_list.sort()
        cluster_best_list.append(model_list[0])
    cluster_best_list.sort()
    all_number_list.sort()
    content = ''
    #for i in range(1,1001):
    #    if i in cluster_best_list:
    #        content += '%s\n' % i
    #    elif i not in all_number_list:
    #        content += '%s\n' % i
    for i in cluster_best_list:
        content += '%s\n' % i
    with open(file_output, 'wb') as cluster_file_output:
        cluster_file_output.write(content)

#extract top 3 conformations from each cluster
def extract_top3_cluster(cluster_input, top1_intput, result_input, file_output):
    with open(cluster_input, 'rb') as cluster_file:
        cluster_lines = cluster_file.read().splitlines()
    with open(top1_intput, 'rb') as top1_file:
        top1_lines = top1_file.read().splitlines()
    top1_list = []
    for top1_line in top1_lines:
        top1_list.append(int(top1_line))
    cluster_dict = {}
    all_top3_list = []
    all_cluster_list = []
    for cluster_line in cluster_lines:
        cluster_items = cluster_line.split()
        # get the number of model
        cluster_list = []
        for i in range(3,len(cluster_items)):
            cluster_list.append(int(cluster_items[i]))
            all_cluster_list.append(int(cluster_items[i]))
        # get the best model from cluster
        cluster_list.sort()
        if cluster_list[0] in top1_list:
            for j in range(0,3):
                cluster_dict[cluster_list[j]] = cluster_list[0]
                all_top3_list.append(cluster_list[j])
    with open(result_input, 'rb') as result_file:
        result_lines = result_file.read().splitlines()
    content = 'Rank Interface energy Best rank from the same cluster\n'
    for result_line in result_lines:
        number = int(result_line.split()[0])
        if number in all_top3_list:
            content += '%s %s\n' % (result_line, cluster_dict[number])
    # For cluster less than 15，extract the models from the conformations which has not been clusterd
    if len(all_top3_list) < 45:
        cluster_num = len(all_top3_list) / 3
        for result_line in result_lines:
            number = int(result_line.split()[0])
            if number not in all_cluster_list and cluster_num < 15:
                content += '%s %s\n' % (result_line, number)
                cluster_num += 1

    with open(file_output, 'wb') as cluster_file_output:
        cluster_file_output.write(content)

#voromqa
def voromqa(file_input_pdb):
    os.system(VOROMQA + '/bin/voronota-voromqa --input %s --contacts-query \'--no-same-chain --no-solvent\' '
                        '--print-energy-of-contacts-selection >> results_voromqa' % (file_input_pdb))

#rank voromqa score according to the interface score
def rank_voromqa(file_voromqa_score, file_out):
    os.system('sort -k10 -n %s | awk \'{if($10<0 && ! (! $10)) print $0}\' '
              '| nl > %s' % (file_voromqa_score, file_out))

#split mol2 file
def split_mol2_obenergy(input_file, output_filepath, obenergy_file):
    with open(input_file, 'rb') as mol2_file:
        mol2_lines = mol2_file.read().splitlines()
    with open(obenergy_file, 'rb') as ob_file:
        obenergy_lines = ob_file.read().splitlines()
    list_obenergy = []
    for obenergy_line in obenergy_lines:
        obenergy_num = int(obenergy_line.split()[0])
        list_obenergy.append(obenergy_num)
    ligand_name = input_file.split('/')[-1][:-5]
    content = ''
    i=0
    for mol2_line in mol2_lines:
        if mol2_line == '@<TRIPOS>MOLECULE':
            i += 1
            content = ''
        content += mol2_line + '\n'
        if i in list_obenergy:
            with open('%s/%s_%s.mol2' % (output_filepath, ligand_name, i), 'wb') as filter_file:
                filter_file.write(content)

#Convert the file format by openbabel
def obabel_convert_format(iformat, file_input, oformat, file_out, addH = False):
    if addH:
        os.system(ADFRSUITE + '/bin/obabel -h -i%s %s -o%s -O %s' % (iformat, file_input, oformat, file_out))
    else:
        os.system(ADFRSUITE + '/bin/obabel -i%s %s -o%s -O %s' % (iformat, file_input, oformat, file_out))

#Assign bond order for ligands in pdb file according to the smiles.
#Ligands in pdb file usually meets wrong bond orders.
def pdb_AssignBondOrder(file_input_pdb, file_smi, file_out):
    ref_smi = Chem.MolFromSmiles(file_smi)
    pdb = Chem.MolFromPDBFile(file_input_pdb)
    pdb = AllChem.AssignBondOrdersFromTemplate(ref_smi, pdb)
    writer = Chem.SDWriter(file_out)
    writer.write(pdb)

#get the small molecule from PDB file
def preprocess_pdb_element(file_input_pdb, file_output_pdb):
    with open(file_input_pdb, 'rb') as file_pdb:
        pdb_lines = file_pdb.read().splitlines()
    content = ''
    for pdb_line in pdb_lines:
        if pdb_line[:6] == 'HETATM' and pdb_line[12:14]:
            element = pdb_line[12:14].translate(None, digits)
            if element[0] == 'H':
                element = ' H'
            elif len(element) == 1:
                element = ' %s' % element
            length = len(pdb_line)
            content += pdb_line[:76]
            diff = 76 - length
            blank_length = ' ' * diff
            content += '%s%s\n' % (blank_length, element)
    with open(file_output_pdb, 'wb') as file_output:
        file_output.write(content)

#adding hydrogens to the protein
def addH_protein(file_input_pdb, file_output_pdb):
    os.system(ADFRSUITE + '/bin/reduce -OH -HIS -NOADjust -NUClear %s 1> %s 2>> addH_log'
              % (file_input_pdb, file_output_pdb))

#get the number of residues around the small molecule
def lig_around_residue(file_input_pdb, out_site):
    with open(file_input_pdb, 'rb') as pdb_file:
        pdb_lines = pdb_file.read().splitlines()
    lig_list = []
    protein_list = []
    lig_dict = {}
    protein_dict = {}
    site_list = []
    for pdb_line in pdb_lines:
        if pdb_line[:4] == 'ATOM':
            chain_id = pdb_line[21]
            res_number = pdb_line[22:27].strip()
            atom_number = pdb_line[12:16].strip()
            x = pdb_line[30:38].strip()
            y = pdb_line[38:46].strip()
            z = pdb_line[46:54].strip()
            protein_dict['chain_id'] = chain_id
            protein_dict['res_number'] = res_number
            protein_dict['atom_number'] = atom_number
            protein_dict['x'] = float(x)
            protein_dict['y'] = float(y)
            protein_dict['z'] = float(z)
            protein_list.append(protein_dict.copy())
        elif pdb_line[:6] == 'HETATM':
            x = float(pdb_line[30:38].strip())
            y = float(pdb_line[38:46].strip())
            z = float(pdb_line[46:54].strip())
            lig_dict['x'] = x
            lig_dict['y'] = y
            lig_dict['z'] = z
            lig_list.append(lig_dict.copy())
    for pdb_line in pdb_lines:
        if pdb_line[:4] == 'ATOM':
            chain_id = pdb_line[21]
            res_number = pdb_line[22:27].strip()
            res = chain_id+res_number
            atom_number = pdb_line[12:16].strip()
            x = float(pdb_line[30:38].strip())
            y = float(pdb_line[38:46].strip())
            z = float(pdb_line[46:54].strip())
            if res not in site_list:
                for i in range(len(lig_list)):
                    lig_x = lig_list[i]['x']
                    lig_y = lig_list[i]['y']
                    lig_z = lig_list[i]['z']
                    length = (x - lig_x) * (x - lig_x) + (y - lig_y) * (y - lig_y) + (z - lig_z) * (z - lig_z)
                    if length <= 25:
                        if res not in site_list:
                            site_list.append(res)
    content = ''
    for item in site_list:
        content += "%s\n" % item
    with open(out_site, 'wb') as file_out:
        file_out.write(content)
    return site_list

#judge whether there is contact
def filter_interface_residue(input_info):
    items = input_info.split()
    content = "none"
    if items[1] != 'No' and items[2] != 'No' and int(items[1]) > 1 and int(items[2]) > 1:
        content = items[0]
    return content

#get the PROTAC conformations
def getConformers(file_rec_lig_sdf, file_target_sdf, protac_smi, file_docked,
                  file_out, target_smi='none', rec_smi='none'):
    rmsList = []
    #rdkit might meet some errors for some ligands
    try:
        rec_ligand = Chem.SDMolSupplier(file_rec_lig_sdf)[0]
        target = Chem.SDMolSupplier(file_target_sdf)[0]
        if rec_smi !='none':
            rec_smi = rec_smi
        else:
            rec_smi = Chem.MolToSmiles(rec_ligand)
        if target_smi != 'none':
            target_smi = target_smi
        else:
            target_smi = Chem.MolToSmiles(target)
        ref_docked = Chem.MolFromSmiles('%s.%s' % (rec_smi, target_smi))
        docked_head = Chem.MolFromPDBFile(file_docked)
        docked_head = AllChem.AssignBondOrdersFromTemplate(ref_docked,
                                                           docked_head)
        with open(protac_smi,'rb') as protac_smi_input:
            protac_smi = protac_smi_input.read().splitlines()[0]
        protac = Chem.MolFromSmiles(protac_smi)
        Chem.AddHs(protac)
        Chem.AddHs(docked_head)
        docked_rec = docked_head.GetSubstructMatch(rec_ligand)
        docked_target = docked_head.GetSubstructMatch(target)
        protac_rec = protac.GetSubstructMatch(rec_ligand)
        protac_target = protac.GetSubstructMatch(target)
        if len(docked_target) == 0 or len(docked_rec) == 0 or len(protac_rec) == 0 or len(protac_target) == 0:
            print "The smiles of PROTAC doesn't match the structures of ligands of target or receptor proteins."
        #print docked_rec
        #print docked_target
        #print protac_rec
        #print protac_target
        protac_align_id = list(protac_rec)+list(protac_target)
        docked_align_id = list(docked_rec)+list(docked_target)
        if not (len(docked_rec) == 0 or len(docked_target) == 0):
            cmap = {protac_rec[j]: docked_head.GetConformer().GetAtomPosition(docked_rec[j]) for j in range(len(docked_rec))}
            cmap.update({protac_target[j]: docked_head.GetConformer().GetAtomPosition(docked_target[j]) for j in range(len(docked_target))})
            cids = AllChem.EmbedMultipleConfs(protac, numConfs=100, coordMap=cmap, maxAttempts=1000, numThreads=1)
            if len(cids) > 0:
                writer = Chem.SDWriter(file_out)
                for i in range(len(cids)):
                    rms = rdMolAlign.AlignMol(protac, docked_head, prbCid=i,atomMap=zip(protac_align_id,docked_align_id))
                    #rms_rec = rdMolAlign.AlignMol(protac, docked_rec, prbCid=i,atomMap=zip(list(protac_rec),list(docked_rec)))
                    #rms_target = rdMolAlign.AlignMol(protac, docked_target, prbCid=i,atomMap=zip(list(protac_target),list(docked_target)))
                    #print rms
                    if rms < 0.5:
                    #if rms_rec < 0.5 and rms_target < 0.5:
                        rmsList.append(rms)
                        writer.write(protac, confId=i)
                    #rmsList.append(rms)
                    #writer.write(protac, confId=i)
        #content = '%s %s' % (file_docked, len(rmsList))
        #print '%s %s' % (file_docked, len(rmsList))
        return len(rmsList)
    except:
        return len(rmsList)

#extract the frodock results for rosetta
def extract_frodock_result(voromqa_input, cluster_input, model_input):
    with open(voromqa_input, 'rb') as voromqa_input_file:
        voromqa_lines = voromqa_input_file.read().splitlines()
    with open(cluster_input, 'rb') as cluster_input_file:
        cluster_lines = cluster_input_file.read().splitlines()
    cluster_list = []
    content_list = []
    model_list = []
    for cluster_line in cluster_lines:
        cluster_list.append(cluster_line)
    cluster_num = len(cluster_list)
    for voromqa_line in voromqa_lines:
        items = voromqa_line.split()
        rank = items[0]
        model_num = items[1].split('.')[1]
        if rank in cluster_list:
            content_list.append('%s %s' % (rank, model_num))
    # For cluster less than 15，extract the models from the conformations which has not been clusterd
    if cluster_num < 15:
        with open(model_input, 'rb') as model_input_file:
            model_lines = model_input_file.read().splitlines()
        for model_line in model_lines:
            model_items = model_line.split()
            for i in range(3,len(model_items)):
                model_list.append(model_items[i])
        for voromqa_line in voromqa_lines:
            items = voromqa_line.split()
            rank = items[0]
            model_num = items[1].split('.')[1]
            if rank not in model_list and cluster_num < 15:
                content_list.append('%s %s' % (rank, model_num))
                cluster_num += 1
    return content_list
