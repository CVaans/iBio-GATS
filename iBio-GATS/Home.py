# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Home.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

import pandas as pd

from SearchTemplate import SearchTemplate
from helper import Helper

df_targ = pd.read_excel('or_data.xlsx', sheet_name="Targets")

class HomeWindow(object):
    def __init__(self, sequence_target):
        self.sequence_target_new = sequence_target
        self.search()

    def example(self):
        self.target_sequence.setText(
            ">OR1A2\n"
            "MKKENQSFNLDFILLGVTSQQEQNNVFFVIFLCIYPITLTGNLLIILAICADIRLHNPMY\n"
            "FLLANLSLVDIIFSSVTIPKVLANHLLGSKFISFGGCLMQMYFMIALAKADSYTLAAMAY\n"
            "DRAVAISCPLHYTTIMSPRSCILLIAGSWVIGNTSALPHTLLTASLSFCGNQEVANFYCD\n"
            "IMPLLKLSCSDVHFNVKMMYLGVGVFSLPLLCIIVSYVQVFSTVFQVPSTKSLFKAFCTC\n"
            "GSHLTVVFLYYGTTMGMYFRPLTSYSPKDAVITVMYVAVTPALNPFIYSLRNWDMKAALQ\n"
            "KLFSKRISS")

    def search(self):
        seq_id = self.preprocess_input()
        print("The given OR sequence is {}".format(seq_id) + "\n")
        if not seq_id:
            return

        SearchTemplate(seq_id)


    def preprocess_input(self):
        # read input from text editor
        sequence = self.sequence_target_new.strip()
        sequence = Helper.sanitize_input(sequence)
        #print('*' * 80)
        print("{}".format(sequence))
        #print('*' * 80)

        # check if protein sequence is valid
        if not Helper.validate_input(sequence):
            print("Warning", "Please enter a valid sequence")
            return False

        # get unitprod_id of protein sequence
        #seq_id = Helper.get_unitprot_id(sequence)
        seq_id = Helper.get_receptor(sequence)
        if not seq_id:
            print("Warning", "Sequence not found in target file")
            return False

        return seq_id

    def helix_numbers(self, tup, df):
        name = tup[0]
        helix = []
        BW = []
        for i in range(1, 8):
            beg = 'TM' + str(i) + ' start'
            end = 'TM' + str(i) + ' end'
            bw = 'BW' + str(i) + '.50'
            tupbeg = df.loc[df['Uniprot_ID'] == name, beg].iloc[0]
            tupend = df.loc[df['Uniprot_ID'] == name, end].iloc[0]
            bw1 = df.loc[df['Uniprot_ID'] == name, bw].iloc[0]
            BW.append(int(bw1[1:]))
            helix.append((tupbeg, tupend))
            return (BW, helix)

    def split_helix(self, id, type, numbers):
        if type == 'template':
            ref = 'PDBID'
        elif type == 'target':
            ref = 'Uniprot_ID'
        prot = [[None for x in range(7)], [None for y in range(7)]]
        numb = [[None for x in range(7)], [None for y in range(7)]]
        seq = df_targ.loc[df_targ['Uniprot_ID'] == id, 'Sequence'].iloc[0]
        for i in range(1, 8):
            beg = int(df_targ.loc[df_targ['Uniprot_ID'] == id, 'TM' + str(i) + ' start'].iloc[0])
            end = int(df_targ.loc[df_targ['Uniprot_ID'] == id, 'TM' + str(i) + ' end'].iloc[0])
            numb[0][i - 1] = [beg, end]
            prot[0][i - 1] = seq[(int(beg) - 1):int(end)]
            bw = (df_targ.loc[df_targ['Uniprot_ID'] == id, 'BW' + str(i) + '.50'].iloc[0])
            numb[1][i - 1] = int(bw[1:])
            prot[1][i - 1] = int(bw[1:]) - int(beg)

        TM1 = prot[0][0]
        TM2 = prot[0][1]
        TM3 = prot[0][2]
        TM4 = prot[0][3]
        TM5 = prot[0][4]
        TM6 = prot[0][5]
        TM7 = prot[0][6]

        # self.window = QtWidgets.QMainWindow()
        # self.ui = Ui_BrowseTemplates(self.window, TM1, TM2, TM3, TM4, TM5, TM6, TM7, id)
        # self.ui.setup_ui()
        # self.home_window.hide()
        # self.home_window.setWindowTitle('Browse Template - Bio-GATS')
        # self.window.show()
