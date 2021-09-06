#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os, sys
import utils.preprocess as pre
from threading import Thread
from Queue import Queue

ROSETTA = os.environ["ROSETTA"]
ADFRSUITE = os.environ['ADFRSUITE']

#run FRODOCK
def rosetta(cpu, lig_locate_num):
    print ('ROSETTA docking')
    #working directory is in filepath_rosetta

    #preprocess files for rosetta docking
    filepath_frodock = '../frodock'
    filepath_frodock_results = '../frodock_results'
    filepath_frodock_cluster = '../frodock/cluster'
    frodock_voromqa = '%s/results_voromqa_interface_rank' % filepath_frodock
    frodock_cluster = '%s/cluster/cluster_5_3_results_top1' % filepath_frodock
    frodock_model = '%s/cluster/cluster_5_3_model_process' % filepath_frodock
    frodock_list = extract_frodock_result(frodock_voromqa, frodock_cluster, frodock_model)
    os.system('cp %s/{rec_lig_*.sdf,rec_lig.sdf,target_lig.sdf,protac.smi} .' % (filepath_frodock))
    molfile_to_params('LG1', 'rec_lig', 'rec_lig_H.sdf') #generate the parameters of rec_lig for rosetta
    pre.alter_chain('rec_lig.pdb', 'rec_lig_Y.pdb', 'Y')
    rec_id_list = pre.get_chain_id('%s/receptor.pdb' % filepath_frodock)
    rec_id = ''.join(rec_id_list) + 'Y'
    target_id_list = pre.get_chain_id('%s/target.pdb' % filepath_frodock)
    target_id = ''.join(target_id_list) + 'X'
    rosetta_partner = rec_id + '_' + target_id
    file_docking_flags = 'docking.flags'
    generate_rosetta_para(file_docking_flags)
    rosetta_mpi_icc = "%s/main/source/bin/docking_protocol.mpi.linuxiccrelease" % ROSETTA
    rosetta_mpi_gcc = "%s/main/source/bin/docking_protocol.mpi.linuxgccrelease" % ROSETTA
    rosetta_icc = "%s/main/source/bin/docking_protocol.linuxiccrelease" % ROSETTA
    rosetta_gcc = "%s/main/source/bin/docking_protocol.linuxgccrelease" % ROSETTA
    flag_mpi = 0
    if os.path.exists(rosetta_mpi_icc):
        rosetta_docking = rosetta_mpi_icc
        flag_mpi = 1
    elif os.path.exists(rosetta_mpi_gcc):
        rosetta_docking = rosetta_mpi_gcc
        flag_mpi = 1
    elif os.path.exists(rosetta_icc):
        rosetta_docking = rosetta_icc
    elif os.path.exists(rosetta_gcc):
        rosetta_docking = rosetta_gcc
    else:
        print "No RosettaDock is found"
        sys.exit()
    for j in range(len(frodock_list)):
        new_rank = frodock_list[j].split()[0]
        old_rank = frodock_list[j].split()[1]
        pre.preprocess_pdb_element('%s/target_%s_addH.pdb'% (filepath_frodock, old_rank),
                                   'target_lig_pre.%s.pdb' % (new_rank))
        pre.obabel_convert_format('pdb', 'target_lig_pre.%s.pdb' % (new_rank),
                                  'sdf', 'target_lig_pre.%s.sdf' % (new_rank), addH=True)
        molfile_to_params('LG2', 'lig.%s' % new_rank, 'target_lig_pre.%s.sdf' % (new_rank))
        os.system('cp %s/cluster/model.%s.pdb .' % (filepath_frodock, new_rank))
        os.system('cat rec_lig_Y.pdb lig.%s.pdb >> model.%s.pdb' % (new_rank, new_rank))

        #rosetta docking
        if flag_mpi == 1:
            os.system('mpirun -np %s %s @%s -s model.%s.pdb '
                      '-extra_res_fa lig.%s.params -partners %s -out:file:scorefile model_%s.sc >/dev/null 2>&1'
                      % (cpu, rosetta_docking, file_docking_flags, new_rank, new_rank, rosetta_partner, new_rank))
        elif flag_mpi == 0:
            os.system('%s @%s -s model.%s.pdb '
                      '-extra_res_fa lig.%s.params -partners %s -out:file:scorefile model_%s.sc >/dev/null 2>&1'
                      % (rosetta_docking, file_docking_flags, new_rank, new_rank, rosetta_partner, new_rank))


    #preprocess for filtering
    with open('%s/rec_site' % filepath_frodock, 'rb') as rec_file:
        rec_site = rec_file.read().splitlines()
    with open('%s/lig_site' % filepath_frodock, 'rb') as lig_file:
        lig_site = lig_file.read().splitlines()
    rosetta_list = 'rosetta_list'
    os.system('ls model.*_*.pdb > %s' % rosetta_list) #get the conformation list
    with open(rosetta_list, 'rb') as rosetta_list_file:
        rosetta_list_lines = rosetta_list_file.read().splitlines()
    filepath_vina = 'vina'
    #preprocess the ligand with two possible locations for linker
    filepath_rec_lig_1 = 'rec_lig_1'
    filepath_rec_lig_2 = 'rec_lig_2'
    if lig_locate_num > 1:
        pre.obabel_convert_format('sdf', 'rec_lig_1.sdf', 'pdb', 'rec_lig_1.pdb')
        pre.obabel_convert_format('sdf', 'rec_lig_2.sdf', 'pdb', 'rec_lig_2.pdb')
        if os.path.exists('%s/vina' % filepath_rec_lig_1) == False:
            os.makedirs('%s/vina' % filepath_rec_lig_1)
        if os.path.exists('%s/vina' % filepath_rec_lig_2) == False:
            os.makedirs('%s/vina' % filepath_rec_lig_2)
    else:
        if os.path.exists(filepath_vina) == False:
            os.makedirs(filepath_vina)
        else:
            os.system('rm %s/*' % filepath_vina)

    #Run interface_residue, obenergy, vina, voromqa
    filtering(cpu, rec_site, rec_id_list, lig_site, target_id_list, rosetta_list_lines, filepath_vina,
              filepath_rec_lig_1, filepath_rec_lig_2, lig_locate_num)

    #preprocess for clustering
    results_voromqa_rosetta = 'results_voromqa'
    results_voromqa_frodock = '%s/all/results_frodock.txt' % filepath_frodock_results
    results_rosetta_top50 = 'results_rosetta_top50'
    rosetta_cluster_topN(results_voromqa_rosetta, results_voromqa_frodock, results_rosetta_top50, 50)
    with open(results_rosetta_top50, 'rb') as voromqa_score:
       score_lines = voromqa_score.read().splitlines()
    filepath_cluster = 'cluster'
    if os.path.exists(filepath_cluster) == False:
       os.makedirs(filepath_cluster)
    else:
       os.system('rm %s/*' % filepath_cluster)
    pdb_model = ''
    # The ligand with only one location for linker
    if lig_locate_num == 1:
        vina_dict = {}
        vina_score_filter = '%s/score_all_top1' % filepath_vina
        with open(vina_score_filter, 'rb') as file_score:
           vina_score_lines = file_score.read().splitlines()
        for vina_score_line in vina_score_lines:
           items = vina_score_line.split()[0].split('_')
           vina_dict['%s_%s' % (items[0], items[1])] = items[2]
    else:
        vina_dict_1 = {}
        vina_score_filter_1 = '%s/%s/score_all_top1' % (filepath_rec_lig_1, filepath_vina)
        with open(vina_score_filter_1, 'rb') as file_score:
            vina_score_lines_1 = file_score.read().splitlines()
        for vina_score_line_1 in vina_score_lines_1:
            items = vina_score_line_1.split()[0].split('_')
            vina_dict_1['%s_%s' % (items[0], items[1])] = items[2]
        vina_dict_2 = {}
        vina_score_filter_2 = '%s/%s/score_all_top1' % (filepath_rec_lig_2, filepath_vina)
        with open(vina_score_filter_2, 'rb') as file_score:
            vina_score_lines_2 = file_score.read().splitlines()
        for vina_score_line_2 in vina_score_lines_2:
            items = vina_score_line_2.split()[0].split('_')
            vina_dict_2['%s_%s' % (items[0], items[1])] = items[2]
    for score_line in score_lines:
       score_rank = score_line.split()[0]
       model_pdb_id = score_line.split()[1]
       os.system('grep ATOM model.%s.pdb > %s/model.%s.pdb' % (model_pdb_id, filepath_cluster, score_rank))
       pre.preprocess_cluster('%s/model.%s.pdb' % (filepath_cluster, score_rank),
                              '%s/model.%s_chainxseg.pdb' % (filepath_cluster, score_rank))
       pdb_model += 'model.%s_chainxseg.pdb\n' % score_rank
       #merge teh best protac and protein
       # Judge whether this pdb_id comes from rosetta results.
       # If it is from frodock, we should copy the protac from frodock results
       # The ligand with only one location for linker
       if lig_locate_num == 1:
           if vina_dict.has_key('%s' % model_pdb_id):
               protac_best_num = vina_dict['%s' % model_pdb_id]
               protac_best_mol2 = 'protac_%s_%s.mol2' % (model_pdb_id, protac_best_num)
               os.system('cp %s %s/protac_%s.mol2' % (protac_best_mol2, filepath_cluster, score_rank))
           else:
               os.system('cp %s/protac_%s.mol2 %s/protac_%s.mol2' %
                         (filepath_frodock_cluster, model_pdb_id, filepath_cluster, score_rank))
           pre.obabel_convert_format('mol2', '%s/protac_%s.mol2' % (filepath_cluster, score_rank),
                                     'pdb', '%s/protac_%s.pdb' % (filepath_cluster, score_rank))
           os.system('cat %s/model.%s.pdb %s/protac_%s.pdb > %s/model_merge_%s.pdb' %
                     (filepath_cluster, score_rank, filepath_cluster, score_rank, filepath_cluster, score_rank))
       else:
           os.system('cat %s/model.%s.pdb > %s/model_merge_%s.pdb' %
                     (filepath_cluster, score_rank, filepath_cluster, score_rank))
           flag = 0
           if vina_dict_1.has_key(model_pdb_id):
               protac_best_num_1 = vina_dict_1['%s' % model_pdb_id]
               protac_best_mol2_1 = '%s/protac_%s_%s.mol2' % (filepath_rec_lig_1, model_pdb_id, protac_best_num_1)
               os.system('cp %s %s/protac_%s_1.mol2' % (protac_best_mol2_1, filepath_cluster, score_rank))
               pre.obabel_convert_format('mol2', '%s/protac_%s_1.mol2' % (filepath_cluster, score_rank),
                                         'pdb', '%s/protac_%s_1.pdb' % (filepath_cluster, score_rank))
               os.system('cat %s/protac_%s_1.pdb >> %s/model_merge_%s.pdb' %
                         (filepath_cluster, score_rank, filepath_cluster, score_rank))
               flag = 1
           if vina_dict_2.has_key(model_pdb_id):
               protac_best_num_2 = vina_dict_2['%s' % model_pdb_id]
               protac_best_mol2_2 = '%s/protac_%s_%s.mol2' % (filepath_rec_lig_2, model_pdb_id, protac_best_num_2)
               os.system('cp %s %s/protac_%s_2.mol2' % (protac_best_mol2_2, filepath_cluster, score_rank))
               pre.obabel_convert_format('mol2', '%s/protac_%s_2.mol2' % (filepath_cluster, score_rank),
                                         'pdb', '%s/protac_%s_2.pdb' % (filepath_cluster, score_rank))
               os.system('cat %s/protac_%s_2.pdb >> %s/model_merge_%s.pdb' %
                         (filepath_cluster, score_rank, filepath_cluster, score_rank))
               flag = 1
           if flag == 0:
               os.system('cp %s/model_merge_%s.pdb %s/model_merge_%s.pdb'  %
                         (filepath_frodock_cluster, model_pdb_id, filepath_cluster, score_rank))

    with open('%s/pdb_model.list' % filepath_cluster, 'wb') as file_out:
       file_out.write(pdb_model)

    #cluster
    pre.cluster('pdb_model.list', filepath_cluster, cpu)

    #Output results
    results_rosetta= 'results_rosetta.txt'
    cluster_top3 = 'cluster_top3.txt'
    filepath_rosetta_results = '../rosetta_results'
    if os.path.exists(filepath_rosetta_results) == False:
       os.makedirs(filepath_rosetta_results)
    os.system('awk \'{print $1\" \"$3}\' %s > %s' % (results_rosetta_top50, results_rosetta))
    pre.extract_top3_cluster('%s/cluster_5_3_model_process' % filepath_cluster,
                            '%s/cluster_5_3_results_top1' % filepath_cluster, results_rosetta,
                            '%s/%s' % (filepath_rosetta_results, cluster_top3))
    with open('%s/%s' % (filepath_rosetta_results, cluster_top3), 'rb') as file:
       cluster_lines = file.read().splitlines()
    for cluster_line in cluster_lines:
       rank = cluster_line.split()[0]
       cluster_num = cluster_line.split()[-1]
       os.system('cp %s/model_merge_%s.pdb %s/cluster_%s_%s.pdb' % (
           filepath_cluster, rank, filepath_rosetta_results, cluster_num, rank))
    filepath_rosetta_results_all = '%s/all' % filepath_rosetta_results
    if os.path.exists(filepath_rosetta_results_all) == False:
       os.makedirs(filepath_rosetta_results_all)
    os.system('cp %s/model_merge_*.pdb %s %s' % (filepath_cluster, results_rosetta, filepath_rosetta_results_all))
    os.system('rm target_lig_pre.* target_lig.*_*.pdb target_lig.*_*.sdf protac_*_*.mol2 *.pdbqt model_nolig.*.pdb')

#get the topN conformations from each ensemble of Rosetta and merge the frodock result
def rosetta_cluster_topN(input_voro_score, input_frodock_score, output_file, number):
    content = ''
    with open(input_frodock_score, 'rb') as frodock_file:
        frodock_lines = frodock_file.read().splitlines()
    for frodock_line in frodock_lines:
        frodock_model = frodock_line.split()[0]
        if os.path.exists('model.%s.pdb' % frodock_model):
            content += '%s\n' % frodock_line
    voro_score_sort = '%s_sort' % input_voro_score
    os.system('sort -k10 -n %s | awk \'{if($10<0 && ! (! $10)) print $1\" \"$10}\' > %s'
              % (input_voro_score, voro_score_sort))
    with open(voro_score_sort, 'rb') as rank_file:
        rank_lines = rank_file.read().splitlines()
    rank_dict = {}
    i = 0
    for rank_line in rank_lines:
        rank_items = rank_line.split()
        score = rank_items[1]
        rank_model_num = rank_items[0].split('.')[1]
        rank_model = rank_items[0].split('.')[1].split('_')[0]
        if rank_model in rank_dict:
            if rank_dict[rank_model] < int(number):
                i += 1
                rank_dict[rank_model] = rank_dict[rank_model] + 1
                #content += '%s %s\n' % (i, rank_line)
                content += '%s %s\n' % (rank_model_num, score)
        else:
            i += 1
            rank_dict[rank_model] = 1
            #content += '%s %s\n' % (i, rank_line)
            content += '%s %s\n' % (rank_model_num, score)
    output_file_temp = '%s_merge' % voro_score_sort
    with open(output_file_temp, 'wb') as file_out:
        file_out.write(content)
    os.system('sort -n -k2 %s | nl > %s ' % (output_file_temp, output_file))

#get the interface residue
def interface_rosetta(rec_site, lig_site, pdb, rec_chain_list, lig_chain_list):
    with open(pdb, 'rb') as pdb_file:
        pdb_lines = pdb_file.read().splitlines()
    rec_list = []
    rec_site_info_list = []
    rec_dict = {}
    lig_dict = {}
    lig_list = []
    lig_site_info_list = []
    rec_interface_list = []
    lig_interface_list = []
    for pdb_line in pdb_lines:
        if pdb_line[:4] == 'ATOM':
            chain_id = pdb_line[21]
            if chain_id in rec_chain_list:
                res_number = pdb_line[22:27].strip()
                res = chain_id + res_number
                atom_number = pdb_line[12:16].strip()
                x = pdb_line[30:38].strip()
                y = pdb_line[38:46].strip()
                z = pdb_line[46:54].strip()
                rec_dict['res'] = res
                rec_dict['atom_number'] = atom_number
                rec_dict['x'] = float(x)
                rec_dict['y'] = float(y)
                rec_dict['z'] = float(z)
                if res in rec_site:
                    rec_site_info_list.append(rec_dict.copy())
                rec_list.append(rec_dict.copy())
            elif chain_id in lig_chain_list:
                res_number = pdb_line[22:27].strip()
                res = chain_id + res_number
                atom_number = pdb_line[12:16].strip()
                x = pdb_line[30:38].strip()
                y = pdb_line[38:46].strip()
                z = pdb_line[46:54].strip()
                lig_dict['res'] = res
                lig_dict['atom_number'] = atom_number
                lig_dict['x'] = float(x)
                lig_dict['y'] = float(y)
                lig_dict['z'] = float(z)
                if res in lig_site:
                    lig_site_info_list.append(lig_dict.copy())
                lig_list.append(lig_dict.copy())
    for i in range(len(lig_site_info_list)):
        lig_res = lig_site_info_list[i]['res']
        lig_x = lig_site_info_list[i]['x']
        lig_y = lig_site_info_list[i]['y']
        lig_z = lig_site_info_list[i]['z']
        for j in range(len(rec_site_info_list)):
            rec_res = rec_site_info_list[j]['res']
            rec_x = rec_site_info_list[j]['x']
            rec_y = rec_site_info_list[j]['y']
            rec_z = rec_site_info_list[j]['z']
            length = (rec_x - lig_x) * (rec_x - lig_x) + (rec_y - lig_y) * (rec_y - lig_y) + (rec_z - lig_z) * (
                    rec_z - lig_z)
            if length <= 25:
                if rec_res not in rec_interface_list:
                    rec_interface_list.append(rec_res)
                if lig_res not in lig_interface_list:
                    lig_interface_list.append(lig_res)
    if len(rec_interface_list) > 0:
        content = '%s %s %s' % (pdb, len(lig_interface_list), len(rec_interface_list))
        #print('%s %s %s' % (pdb, len(lig_interface_list), len(rec_interface_list)))
    else:
        content = '%s No No' % (pdb)
        #print('%s No No' % (pdb))
    return content

#Use the molfile_to_params.py from Rosetta
def molfile_to_params(res_name, para_name, sdf_name):
    os.system(ROSETTA + '/main/source/scripts/python/public/molfile_to_params.py -n %s -p %s '
                        '--conformers-in-one-file %s --clobber' % (res_name, para_name, sdf_name))

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

#generate docking parametes for rosetta
def generate_rosetta_para(file_out):
    content = '-in:file:extra_res_fa rec_lig.params\n'
    content += '-nstruct 400\n'
    content += '-use_input_sc\n'
    content += '-load_PDB_components false\n'
    content += '-dock_pert 3 8\n'
    content += '-dock_mcm_trans_magnitude 0.1\n' #refinement translational perturbation
    content += '-dock_mcm_rot_magnitude 5.0\n' #refinement rotational perturbation
    content += '-ex1\n'
    content += '-ex2aro\n'
    content += '-overwrite\n'
    with open(file_out, 'wb') as file:
        file.write(content)

#Run addH, interface_residue, obenergy, vina, voromqa in multithreading
def filtering(cpu, rec_site, rec_id_list, lig_site, target_id_list, rosetta_list_lines, filepath_vina,
              filepath_rec_lig_1, filepath_rec_lig_2, lig_locate_num):
    filtering_queue = Filtering_queue(rec_site, rec_id_list, lig_site, target_id_list, rosetta_list_lines,
                                      filepath_vina, filepath_rec_lig_1, filepath_rec_lig_2, lig_locate_num)
    queue = Queue(maxsize=int(cpu*2))
    producer_thread = Thread(target=filtering_queue.producer, args=(queue,))
    producer_thread.daemon = True
    producer_thread.start()

    for index in range(int(cpu)):
        consumer_thread = Thread(target=filtering_queue.consumer, args=(queue,))
        consumer_thread.daemon = True
        consumer_thread.start()
    queue.join()

#The queue class for addH, interface_residue, obenergy, vina, voromqa
class Filtering_queue:
    def __init__(self, rec_site, rec_id_list, lig_site, target_id_list, rosetta_list_lines, filepath_vina,
                 filepath_rec_lig_1, filepath_rec_lig_2, lig_locate_num):
        self.rec_site = rec_site
        self.rec_id_list = rec_id_list
        self.lig_site = lig_site
        self.target_id_list = target_id_list
        self.rosetta_list_lines = rosetta_list_lines
        self.filepath_vina = filepath_vina
        self.lig_locate_num = lig_locate_num
        self.filepath_rec_lig_1 = filepath_rec_lig_1
        self.filepath_rec_lig_2 = filepath_rec_lig_2
        self.filepath_vina_1 = '%s/vina' % filepath_rec_lig_1
        self.filepath_vina_2 = '%s/vina' % filepath_rec_lig_2

    def producer(self,in_q):
        while in_q.full() is False:
            for line in self.rosetta_list_lines:
                in_q.put(line)

    def consumer(self, in_q):
        while in_q.empty() is not True:
            pdb = in_q.get()
            #get the interface residue
            content_interface = interface_rosetta(self.rec_site, self.lig_site, 
                                                  pdb, self.rec_id_list, self.target_id_list)
            filter_num = pre.filter_interface_residue(content_interface)

            #get the PROTAC conformations
            if filter_num != 'none':
                pdb_num = filter_num.split('.')[1]
                # The ligand with only one location for linker
                if self.lig_locate_num == 1:
                    target_lig_pdb = 'target_lig.%s.pdb' % pdb_num
                    target_lig_sdf = 'target_lig.%s.sdf' % pdb_num
                    os.system('awk \'{if($1=="HETATM") print $0}\' model.%s.pdb > %s' % (pdb_num, target_lig_pdb))
                    pre.obabel_convert_format('pdb', target_lig_pdb, 'sdf', target_lig_sdf)
                    protac_sdf = 'protac_%s.sdf' % pdb_num
                    num_confor = pre.getConformers('rec_lig.sdf', 'target_lig.sdf', 'protac.smi',
                                                   target_lig_sdf, protac_sdf)
                    #Vina and obenergy
                    if int(num_confor) > 0 :
                        #preprocess files for obenergy and vina
                        protac_mol2 = 'protac_%s.mol2' % (pdb_num)
                        pre.obabel_convert_format('sdf', protac_sdf, 'mol2', protac_mol2)
                        model_nolig_pdb = 'model_nolig.%s.pdb' % (pdb_num)
                        model_pdb = 'model.%s.pdb' % pdb_num
                        os.system('grep ATOM %s > %s' % (model_pdb, model_nolig_pdb))
                        #obenergy and vina
                        num_obenergy_vina = pre.obenergy_vina(pdb_num, self.filepath_vina, '.')
                        #voromqa
                        if int(num_obenergy_vina) > 0:
                            pre.voromqa(model_nolig_pdb)
                #The ligand with two possible location for linker
                else:
                    target_lig_pdb = 'target_lig.%s.pdb' % pdb_num
                    target_lig_sdf = 'target_lig.%s.sdf' % pdb_num
                    os.system('awk \'{if($1=="HETATM" && $5=="X") print $0}\' model.%s.pdb > %s' % (pdb_num, target_lig_pdb))
                    target_lig_pdb_1 = '%s/%s' % (self.filepath_rec_lig_1, target_lig_pdb)
                    target_lig_sdf_1 = '%s/%s' % (self.filepath_rec_lig_1, target_lig_sdf)
                    os.system('cat %s rec_lig_1.pdb > %s'
                              % (target_lig_pdb, target_lig_pdb_1))
                    pre.obabel_convert_format('pdb', target_lig_pdb_1, 'sdf', target_lig_sdf_1)
                    protac_sdf_1 = '%s/protac_%s.sdf' % (self.filepath_rec_lig_1, pdb_num)
                    num_confor_1 = pre.getConformers('rec_lig_1.sdf', 'target_lig.sdf', 'protac.smi',
                                                     target_lig_sdf_1, protac_sdf_1)
                    num_obenergy_vina_1 = 0
                    num_obenergy_vina_2 = 0
                    model_nolig_pdb = 'model_nolig.%s.pdb' % (pdb_num)
                    if int(num_confor_1) > 0 :
                        #preprocess files for obenergy and vina
                        protac_mol2_1 = '%s/protac_%s.mol2' % (self.filepath_rec_lig_1, pdb_num)
                        pre.obabel_convert_format('sdf', protac_sdf_1, 'mol2', protac_mol2_1)
                        model_pdb = 'model.%s.pdb' % pdb_num
                        os.system('grep ATOM %s > %s' % (model_pdb, model_nolig_pdb))
                        #obenergy and vina
                        num_obenergy_vina_1 = pre.obenergy_vina(pdb_num, self.filepath_vina_1, self.filepath_rec_lig_1)

                    target_lig_pdb_2 = '%s/%s' % (self.filepath_rec_lig_2, target_lig_pdb)
                    target_lig_sdf_2 = '%s/%s' % (self.filepath_rec_lig_2, target_lig_sdf)
                    os.system('cat %s rec_lig_2.pdb > %s'
                              % (target_lig_pdb, target_lig_pdb_2))
                    pre.obabel_convert_format('pdb', target_lig_pdb_2, 'sdf', target_lig_sdf_2)
                    protac_sdf_2 = '%s/protac_%s.sdf' % (self.filepath_rec_lig_2, pdb_num)
                    num_confor_2 = pre.getConformers('rec_lig_2.sdf', 'target_lig.sdf', 'protac.smi',
                                                     target_lig_sdf_2, protac_sdf_2)
                    if int(num_confor_2) > 0 :
                        #preprocess files for obenergy and vina
                        protac_mol2_2 = '%s/protac_%s.mol2' % (self.filepath_rec_lig_2, pdb_num)
                        pre.obabel_convert_format('sdf', protac_sdf_2, 'mol2', protac_mol2_2)
                        model_pdb = 'model.%s.pdb' % pdb_num
                        os.system('grep ATOM %s > %s' % (model_pdb, model_nolig_pdb))
                        #obenergy and vina
                        num_obenergy_vina_2 = pre.obenergy_vina(pdb_num, self.filepath_vina_2, self.filepath_rec_lig_2)
                    # voromqa
                    if int(num_obenergy_vina_1) > 0 or int(num_obenergy_vina_2) > 0:
                        pre.voromqa(model_nolig_pdb)

            in_q.task_done()

