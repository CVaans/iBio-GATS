# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'SearchTemplate.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from __future__ import print_function

import matplotlib.pyplot as plt

import pandas as pd
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio import SeqIO
from helper import Helper


class SearchTemplate(object):
    def __init__(self, seq_id):
        self.targ_id = seq_id
        self.hotspot_residue_positions = Helper.non_olfactory_hotspot_residue_positions

        self.template = []
        self.target = []

        self.ress = []
        self.hydrh = []
        self.summ = [0] * 7
        self.value = [0] * 7
        self.tm_value = [0] * 2
        self.plot_data = {}
        self.current_df = ''

        self.index_temp1 = 0
        self.index_temp2 = 0
        self.index_temp3 = 0
        self.state_selected_index = 0
        self.identities = []

        self.temp_df = pd.DataFrame()
        self.scoring_scheme()

    def scoring_scheme(self):
        global selected_template_str
        df4 = Helper.df_target.apply(lambda row: row.astype(str).str.contains(self.targ_id).any(),
                                     axis=1)  # searching entered sequence in target database
        df5, = df4[df4 == True].index  # getting row number of target sequence
        targ_seq = Helper.df_target.loc[df5, "Sequence"]  # getting target sequence

        current_identity = []
        current_resolution = []
        self.current_df = Helper.df_active

        # select resolution with values less than or equal to 2.5
        resolution_index = 0
        if resolution_index == 1:
            self.current_df = self.current_df[self.current_df['Resolution'] <= 2.5]

        # get number of rows in excel sheet (depending on selected value in dropdown i.e. active, inactive, or intermediate)
        number_of_rows_in_sheet = len(self.current_df.axes[0])
        tm1 = [None] * number_of_rows_in_sheet
        tm2 = [None] * number_of_rows_in_sheet
        tm3 = [None] * number_of_rows_in_sheet
        tm4 = [None] * number_of_rows_in_sheet
        tm5 = [None] * number_of_rows_in_sheet
        tm6 = [None] * number_of_rows_in_sheet
        tm7 = [None] * number_of_rows_in_sheet

        self.identities = [None] * number_of_rows_in_sheet
        tm_score_aggregate = [None] * number_of_rows_in_sheet

        similarity_scores = []

        for j in range(number_of_rows_in_sheet):
            current_pdbid = self.current_df['PDBID'][j]
            self.split_helix(self.current_df, current_pdbid, type='template', numbers=True)
            _, self.identities[j] = self.align(self.current_df, current_pdbid, self.targ_id)
            ssd_tm = self.ssd(self.current_df, current_pdbid, self.targ_id)
            self.tm_value[j - 1] = ssd_tm

            tm1[j], tm2[j], tm3[j], tm4[j], tm5[j], tm6[j], tm7[j] = ssd_tm
            iden = Helper.identity_for_whole_sequence(self.current_df['Sequence'][j], targ_seq)

            current_identity.append(iden)
            current_resolution.append(self.current_df['Resolution'][j])
        self.tm_value.reverse()

        
            # similarity score based on hotspot residues not calculated for iBio-GATS
            #sum_value = 0
            #for hotspot_residue_position in self.hotspot_residue_positions:
            #    sum_value += Helper.calculate_hotspot_residue_similarity(targ_seq, self.current_df['Sequence'][j],
                                                                         #hotspot_residue_position)
            #similarity_scores.append(sum_value)

        # calculate SSD score
        self.temp_df = pd.DataFrame()
        self.temp_df['PDBID'] = self.current_df['PDBID']
        self.temp_df['Receptor'] = self.current_df['Receptor']
        self.temp_df['Sequence'] = self.current_df['Sequence']

        for i in range(1, 8):
            self.temp_df[f'TM{i} start'] = self.current_df[f'TM{i} start']
            self.temp_df[f'TM{i} end'] = self.current_df[f'TM{i} end']
            self.temp_df[f'BW{i}.50'] = self.current_df[f'BW{i}.50']

        self.temp_df['TM1_score'] = self.calculate_score(number_of_rows_in_sheet, tm1, 1)
        self.temp_df['TM2_score'] = self.calculate_score(number_of_rows_in_sheet, tm2, 2)
        self.temp_df['TM3_score'] = self.calculate_score(number_of_rows_in_sheet, tm3, 2)
        self.temp_df['TM4_score'] = self.calculate_score(number_of_rows_in_sheet, tm4, 2)
        self.temp_df['TM5_score'] = self.calculate_score(number_of_rows_in_sheet, tm5, 2)
        self.temp_df['TM6_score'] = self.calculate_score(number_of_rows_in_sheet, tm6, 2)
        self.temp_df['TM7_score'] = self.calculate_score(number_of_rows_in_sheet, tm7, 2)

        # calculate normalized value of all TM#_score
        tm_columns = ['TM1_score', 'TM2_score', 'TM3_score', 'TM4_score', 'TM5_score', 'TM6_score', 'TM7_score']
        self.temp_df['tm_score_aggregate'] = self.temp_df[tm_columns].sum(axis=1)
        self.temp_df['tm_score_normalized'] = (self.temp_df['tm_score_aggregate'] - self.temp_df[
            'tm_score_aggregate'].min()) / (self.temp_df['tm_score_aggregate'].max() - self.temp_df[
            'tm_score_aggregate'].min())

        # calculate identity score
        # self.temp_df['identity_score'] = self.calculate_score(number_of_rows_in_sheet, current_identity, reverse=True)

        # calculate resolution score
        self.temp_df['resolution'] = current_resolution
        self.temp_df['resolution_score'] = self.calculate_resolution_score(number_of_rows_in_sheet, current_resolution)
        #doesnot calculate similarity scores for iORs
        # calculate similarity score
        #self.temp_df['similarity_score'] = similarity_scores
        #similarity_scores_df = pd.DataFrame(similarity_scores)
        #self.temp_df['similarity_score_normalized'] = (similarity_scores_df - similarity_scores_df.min()) / (
        #            similarity_scores_df.max() - similarity_scores_df.min())

        # doesnot take hotspot residues in score calculation for iORs
        # calculate normalized score (normalized_score score = tm_score_normalized + similarity_score_normalized + resolution_score)
        columns_to_aggregate = ['tm_score_normalized', 'resolution_score']
        self.temp_df['normalized_score'] = self.temp_df[columns_to_aggregate].sum(axis=1)

        # arrange aggregate scores in descending order
        aggregate_score_sorted_indexes = sorted(range(len(self.temp_df['normalized_score'])),
                                                key=lambda k: self.temp_df['normalized_score'][k], reverse=True)

        # TODO: comment in production
        # print(self.temp_df)
        # print(aggregate_score_sorted_indexes)  # use first three index values ONLY

        # Printing resolution, position and identity of top 3 templates
        if len(aggregate_score_sorted_indexes) > 0:
            self.index_temp1 = aggregate_score_sorted_indexes[0]  # getting the index of best template

            #print("Template 1 Identity (%)")
            #[print("{}".format(str(self.identities[self.index_temp1][i]))) for i in range(7)]
            print('iBio-GATS showing results of the first-best and second-best template for the given target sequence')
            print('*' * 80)

            res1 = self.current_df['Resolution'].iloc[self.index_temp1]  # getting resolution of the best template
            print("Resolution of the first best template : {}".format(str(res1)))

            ident1 = current_identity[self.index_temp1]  # getting identity of the best template
            #print("Sequence Identity of the first best template : {}".format(str(ident1)))

            pos1 = self.current_df['Position'].iloc[self.index_temp1]  # getting position of the best template
            print("Position of the first best template : {}".format(str(pos1)))

           # sim1 = self.temp_df['similarity_score'].iloc[
            #    self.index_temp1]  # getting similarity score of the best template
            #print("Similarity score of the best template 1: {}".format(str(sim1)))

            temp1_id = self.current_df['PDBID'].iloc[self.index_temp1]  # getting pdbid of the top template
            temp1_name = self.current_df['Receptor'].iloc[self.index_temp1]
            #print("Displaying PDBID of the best template 1" + " " + temp1_name + " " + "(PDBID" + ":" + temp1_id + ")")
            print("Displaying the name of first best template: " + " " + temp1_name)
        if len(aggregate_score_sorted_indexes) > 1:
            self.index_temp2 = aggregate_score_sorted_indexes[1]  # getting the index of 2nd best template

            #print("*" * 80 + "\n" + "Template 2 Identity (%)")
            #[print("{}".format(str(self.identities[self.index_temp2][i]))) for i in range(7)]
            print('*' * 80)

            res2 = self.current_df['Resolution'].iloc[self.index_temp2]  # getting resolution of the 2nd best template
            print("Resolution of the second best template : {}".format(str(res2)))

            ident2 = current_identity[self.index_temp2]  # getting identity of the 2nd best template
            #print("Sequence Identity of the second best template : {}".format(str(ident2)))

            pos2 = self.current_df['Position'].iloc[self.index_temp2]  # getting position of the 2nd best template
            print("Position of the second best template : {}".format(str(pos2)))

            #sim2 = self.temp_df['similarity_score'].iloc[
            #    self.index_temp2]  # getting similarity score of the best template
            #print("Similarity score of the best template 2: {}".format(str(sim2)))

            temp2_id = self.current_df['PDBID'].iloc[self.index_temp2]  # getting pdbid of the top template
            temp2_name = self.current_df['Receptor'].iloc[self.index_temp2]
            print("Displaying the name of second best template:  " + " " + temp2_name + "\n" + "*" * 80)

        for idx, current_index in enumerate(aggregate_score_sorted_indexes[:2]):
            print("Template {} SSD values: ".format(idx + 1) + str(tm1[current_index]) + " " + str(tm2[current_index]) + " " + str(tm3[current_index]) + " " + str(tm4[current_index]) + " " + str(tm5[current_index]) + " " + str(tm6[current_index]) + " " + str(tm7[current_index]))
            print("*" * 80)

        for idx, current_index in enumerate(aggregate_score_sorted_indexes[:2]):
            self.download_result_summary(current_index)
            self.download_full_alignment_with_index(current_index)

        #for idx, current_index in enumerate(aggregate_score_sorted_indexes[:2]):
         #   self.download_full_alignment(current_index)

        for idx, current_index in enumerate(aggregate_score_sorted_indexes[:1]):
           # print("*" * 80 + "\n" + "Please enter 1 or 2 to select one of the two templates for model building:" + " " + temp1_name + " " + "(PDBID" + ":" + temp1_id + ")" + " " + "or" + " " + temp2_name + " " + "(PDBID" + ":" + temp2_id + ")")
            print("*" * 80 + "\n" + "Please enter 1 or 2 to select one of the two templates for model building:" + " " +temp1_name + " " + "or" + " " + temp2_name)
           # print("1. " + temp1_name + "  " + "(PDBID" + ":" + temp1_id + ")")
            print("1. " + temp1_name)
            print("2. " + temp2_name)
            user_template_selection = int(input())
            if user_template_selection > 2:
                print("Please enter either 1 or 2")
            elif user_template_selection == 1:
                selected_template_str = temp1_id + "-" + self.targ_id
                lidString = "7LID"
                print("*" * 80 + "\n" + "Please select any one of structures of MhOR5, 7LIG or 7LID for using MhOR5 as template for model building")
                print("1. " + temp1_id)
                print("2. " + lidString)
                user_structure_selection = int(input())
                if user_structure_selection > 2:
                    print("Please enter either 1 or 2")
                elif user_structure_selection == 1:
                    selected_template_str = temp1_id + "-" + self.targ_id
                elif user_structure_selection == 2:
                    selected_template_str = lidString + "-" + self.targ_id

            elif user_template_selection == 2:
                selected_template_str = temp2_id + "-" + self.targ_id

            selected_template_file = open("selected-template.txt", "w")
            # selected_template_file.write(temp1_name + "-" + self.targ_id)
            selected_template_file.write(selected_template_str)
            selected_template_file.close()

    def split_helix(self, current_df, id, type, numbers):
        ref = ''
        df = pd.DataFrame()
        if type == 'template':
            ref = 'PDBID'
            df = current_df
        elif type == 'target':
            ref = 'Uniprot_ID'
            df = Helper.df_target

        prot = [[None for x in range(7)], [None for y in range(7)]]
        numb = [[None for x in range(7)], [None for y in range(7)]]
        seq = df.loc[df[ref] == id, 'Sequence'].iloc[0]

        for i in range(1, 8):
            beg = int(df.loc[df[ref] == id, 'TM' + str(i) + ' start'].iloc[0])
            end = int(df.loc[df[ref] == id, 'TM' + str(i) + ' end'].iloc[0])
            numb[0][i - 1] = [beg, end]
            prot[0][i - 1] = seq[(int(beg) - 1):int(end)]
            bw = (df.loc[df[ref] == id, 'BW' + str(i) + '.50'].iloc[0])
            numb[1][i - 1] = int(bw[1:])
            prot[1][i - 1] = int(bw[1:]) - int(beg)

        return (numb) if numbers else (prot)

    def align(self, df, template_id, target_id):
        self.template = self.split_helix(df, template_id, type='template', numbers=False)
        self.target = self.split_helix(Helper.df_target, target_id, type='target', numbers=False)

        identities = []
        for i in range(7):
            ident = Helper.identity_for_sequence_chunks(self.template[0][i], self.target[0][i])
            identities.append(ident)

        align = [[None for x in range(7)], [None for y in range(7)]]

        for i in range(7):
            if self.template[1][i] > self.target[1][i]:
                align[0][i] = self.template[0][i]
                align[1][i] = ''.join(['-'] * (self.template[1][i] - self.target[1][i])) + self.target[0][i]
            elif self.template[1][i] < self.target[1][i]:
                align[0][i] = ''.join(['-'] * (self.target[1][i] - self.template[1][i])) + self.template[0][i]
                align[1][i] = self.target[0][i]
            else:
                align[0][i] = self.template[0][i]
                align[1][i] = self.target[0][i]

            if len(align[0][i]) > len(align[1][i]):
                align[1][i] = align[1][i] + ''.join(['-'] * (len(align[0][i]) - len(align[1][i])))
            elif len(align[1][i]) > len(align[0][i]):
                align[0][i] = align[0][i] + ''.join(['-'] * (len(align[1][i]) - len(align[0][i])))

        return align, identities

    def show_alignment(self, current_df, index):
        df4 = Helper.df_target.apply(lambda row: row.astype(str).str.contains(self.targ_id).any(),
                                     axis=1)  # searching entered sequence in target database
        df5, = df4[df4 == True].index  # getting row number of target sequence
        seq_id = Helper.df_target.loc[df5, "Uniprot_ID"]
        target_receptor = Helper.df_target.loc[df5, 'Receptor']
        target_sequence = Helper.df_target.loc[df5, 'Sequence']

        pdb_id = current_df['PDBID'].iloc[index]  # getting pdbid of template name
        template_receptor = current_df['Receptor'].iloc[index]
        template_sequence = current_df['Sequence'].iloc[index]

        target_number = []  # list of tuples like ('TM1 start', 'TM1 end'), ('TM2 start', 'TM2 end'), ...
        template_number = []  # list of tuples like ('TM1 start', 'TM1 end'), ('TM2 start', 'TM2 end'), ...

        for i in range(1, 8):
            target_number.append((Helper.df_target.loc[df5, f'TM{i} start'], Helper.df_target.loc[df5, f'TM{i} end']))
            template_number.append((current_df[f'TM{i} start'].iloc[index], current_df[f'TM{i} end'].iloc[index]))

        alignment, _ = self.align(current_df, pdb_id, seq_id)
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_Alignment(self.window, alignment, template_receptor, target_receptor, template_sequence,
                               target_sequence, template_number, target_number)
        self.ui.setup_ui()
        self.window.setWindowTitle('Alignment - Bio-GATS')
        self.window.show()

    def identity(self, template_seq, target_seq):
        matrix = Helper.gpcr_tm
        for a in pairwise2.align.globaldx(template_seq, target_seq, matrix):
            align_temp = a[0]
            align_targ = a[1]
            seq_len = len(align_temp)
            matches = [align_temp[i] == align_targ[i] for i in range(seq_len)]
            iden = (100 * sum(matches)) / seq_len
            return round(iden, 2)
        return None

    def calculate_score(self, number_of_rows_in_sheet, current_list, add_value=1, reverse=False):
        score = [0] * number_of_rows_in_sheet

        for index, value in enumerate(current_list):
            if value > 2.5:
                score[index] = 0
            else:
                score[index] = 1

            if 0 <= value < 0.1:
                score[index] = 2
            elif value >= 0.1:
                score[index] = -1

        return score

    def calculate_resolution_score(self, number_of_rows_in_sheet, current_list):
        score = [0] * number_of_rows_in_sheet

        for index, value in enumerate(current_list):
            if value > 2.5:
                score[index] = 0
            else:
                score[index] = 1

        return score

    def show_hydrophobicity_plot(self, template_text):
        pdbid = template_text.rstrip(')').split(':')[1]
        current_plot_data = self.plot_data[pdbid]

        Helper.show_hydrophobicity_plot(current_plot_data['ress'], current_plot_data['hydrh'],
                                        current_plot_data['summ'], current_plot_data['value'],
                                        current_plot_data['template_receptor'], current_plot_data['target_receptor'],
                                        current_plot_data['name'], download=False)

    def show_helical_plot(self, template_text):
        pdbid = template_text.rstrip(')').split(':')[1]
        current_plot_data = self.plot_data[pdbid]
        current_df = self.plot_data[pdbid]['current_df']

        template_receptor = current_plot_data['template_receptor']
        target_receptor = current_plot_data['target_receptor']

        template = self.split_helix(current_df, pdbid, type='template', numbers=False)
        target = self.split_helix(Helper.df_target, self.targ_id, type='target', numbers=False)

        Helper.show_helical_plot(template, target, template_receptor, target_receptor, current_plot_data['name'],
                                 download=False)


    def ssd(self, current_df, template_id, target_id):
        template_sequence = current_df.loc[current_df['PDBID'] == template_id, 'Sequence'].iloc[0]
        template_receptor = current_df.loc[current_df['PDBID'] == template_id, 'Receptor'].iloc[0]
        target_sequence = Helper.df_target.loc[Helper.df_target['Uniprot_ID'] == target_id, 'Sequence'].iloc[0]
        target_receptor = Helper.df_target.loc[Helper.df_target['Uniprot_ID'] == target_id, 'Receptor'].iloc[0]

        window_size = 11
        hyd = [None, None]
        hyd[0] = [Helper.hydscale[res] for res in template_sequence]
        hyd[1] = [Helper.hydscale[res] for res in target_sequence]
        lens = [len(template_sequence), len(target_sequence)]
        hydr = [[1] * int(window_size / 2) for i in range(2)]
        for k in range(2):
            for i in range(int(window_size / 2), lens[k] - int(window_size / 2)):
                sumb = 0
                for j in range(-int(window_size / 2), int(window_size / 2) + 1):
                    sumb += hyd[k][i + j]
                avg = sumb / window_size
                hydr[k].append(avg)
        for i in range(2):
            hydr[i] += [1] * int(window_size / 2)
        temp_list = list(template_sequence)
        targ_list = list(target_sequence)
        length = len(temp_list)
        targ_align = ['-'] * length
        temp_align = temp_list
        temp_vals = self.split_helix(current_df, template_id, type='template', numbers=True)
        targ_vals = self.split_helix(Helper.df_target, target_id, type='target', numbers=True)
        temp_helix = temp_vals[0]
        targ_helix = targ_vals[0]
        temp_bw = temp_vals[1]
        targ_bw = targ_vals[1]
        name = template_id + ',' + target_id
        hydiff = []
        self.hydrh = [[[] for x in range(2)] for y in range(7)]
        self.ress = [[] for y in range(7)]
        per_res_ssd = []
        for i in range(7):
            self.summ = [0] * 7
            for j in range(int(targ_helix[i][0]) - 1, int(targ_helix[i][1])):
                targ_align[j + temp_bw[i] - targ_bw[i]] = targ_list[j]
            for j in range(int(temp_helix[i][0]) - 1, int(temp_helix[i][1])):
                self.hydrh[i][0].append(hydr[0][j])
                self.hydrh[i][1].append(hydr[1][targ_bw[i] - temp_bw[i] + j])
                self.summ[i] += ((hydr[1][targ_bw[i] - temp_bw[i] + j]) - (hydr[0][j])) ** 2
            resno = (temp_bw[i] - targ_bw[i] + targ_helix[i][0]) - temp_helix[i][0]
            self.ress[i] = [x for x in range(int(temp_helix[i][0]), int(temp_helix[i][1] + 1))]
            self.value[i] = self.summ[i] / (temp_helix[i][1] - temp_helix[i][0] + 1)
            hydiff.append(round(self.summ[i], 2))
            per_res_ssd.append(round(self.value[i], 4))

            temp_pdbid = current_df.loc[current_df['PDBID'] == template_id, 'PDBID'].iloc[0]
            self.plot_data[temp_pdbid] = {
                'current_df': current_df,
                'name': name,
                'ress': self.ress,
                'hydrh': self.hydrh,
                'summ': self.summ,
                'value': self.value,
                'template_receptor': template_receptor,
                'target_receptor': target_receptor,
            }

        return per_res_ssd

    def download_result_summary(self, current_index):
        current_item_text = self.temp_df['PDBID'].iloc[current_index]  # getting PDBID of the top 3 templates
        target_id, target_receptor, target_sequence, target_number, template_id, template_receptor, template_sequence, template_number = Helper.calculate_alignment(
            self.targ_id, current_item_text)

        current_plot_result = self.plot_data[current_item_text]
        name = template_id + ',' + target_id
        alignment, _ = self.align(self.temp_df, template_id, target_id)

        hydrophobicity_plot_filenames = Helper.show_hydrophobicity_plot(current_plot_result['ress'], current_plot_result['hydrh'], self.summ, self.tm_value[current_index],
                                                                        template_receptor, target_receptor, name,
                                                                        download=True)
        helical_plot_filenames = Helper.show_helical_plot(self.template, self.target, template_receptor,
                                                          target_receptor, name, download=True)
        Helper.generate_result_summary(current_index, alignment, template_id, target_id, template_sequence,
                                       target_sequence, template_number, target_number, hydrophobicity_plot_filenames,
                                       helical_plot_filenames)

    def download_full_alignment_with_index(self, current_index):
        current_item_text = self.temp_df['PDBID'].iloc[current_index]  # getting PDBID of the top 3 templates
        target_id, target_receptor, target_sequence, target_number, template_id, template_receptor, template_sequence, template_number = Helper.calculate_alignment(
            self.targ_id, current_item_text)

        name = template_id + ',' + target_id
        alignment, _ = self.align(self.temp_df, template_id, target_id)
        alignment_text = Helper.download_full_alignment(alignment, template_id, target_id, template_sequence,
                                       target_sequence, template_number, target_number)
        # print to file
        #fasta_filename = Helper.full_alignment_filename_fasta.format(current_index)
        if template_id == '7LIG':
            for i in range(2):
                ali_filename = Helper.full_alignment_filename_ali.format(template_id, target_id)
                if i == 0:
                    with open(ali_filename, 'w') as f:
                        f.write(alignment_text)
                if i == 1:
                    alignment_text = Helper.download_full_alignment(alignment, '7LID', target_id,
                                                                    template_sequence,
                                                                    target_sequence, template_number, target_number)
                    ali_filename = Helper.full_alignment_filename_ali.format('7LID', target_id)

                    with open(ali_filename, 'w') as f:
                        f.write(alignment_text)
        else:
            ali_filename = Helper.full_alignment_filename_ali.format(template_id, target_id)
            with open(ali_filename, 'w') as f:
                f.write(alignment_text)





