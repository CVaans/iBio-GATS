# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Alignment.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

from helper import Helper


class Ui_Alignment(object):
    def __init__(self, Alignment, align, template_receptor, target_receptor, template_sequence, target_sequence, template_number, target_number):
        self.alignment = Alignment
        self.tm_align = align
        self.template_receptor = template_receptor
        self.target_receptor = target_receptor
        self.template_sequence = template_sequence
        self.target_sequence = target_sequence
        self.template_number = template_number
        self.target_number = target_number

    def setup_ui(self):
        self.alignment.setObjectName("Alignment")
        self.alignment.resize(850, 600)
        self.centralwidget = QtWidgets.QWidget(self.alignment)
        self.centralwidget.setObjectName("centralwidget")

        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(30, 60, 71, 31))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(30, 200, 71, 31))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(30, 130, 71, 31))
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(30, 270, 71, 31))
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(30, 340, 71, 31))
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(30, 410, 71, 31))
        self.label_6.setObjectName("label_6")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(30, 480, 71, 31))
        self.label_7.setObjectName("label_7")

        self.closeButton = QtWidgets.QPushButton(self.centralwidget)
        self.closeButton.setGeometry(QtCore.QRect(750, 10, 80, 23))
        self.closeButton.setObjectName("homeButton")
        self.closeButton.clicked.connect(self.hide_window)

        self.downloadTmWiseAlignmentButton = QtWidgets.QPushButton(self.centralwidget)
        self.downloadTmWiseAlignmentButton.setGeometry(QtCore.QRect(150, 550, 250, 40))
        self.downloadTmWiseAlignmentButton.setObjectName("downloadTmWiseAlignmentButton")
        self.downloadTmWiseAlignmentButton.clicked.connect(self.download_tm_wise_alignment)

        self.downloadFullAlignmentButton = QtWidgets.QPushButton(self.centralwidget)
        self.downloadFullAlignmentButton.setGeometry(QtCore.QRect(400, 550, 250, 40))
        self.downloadFullAlignmentButton.setObjectName("downloadFullAlignmentButton")
        self.downloadFullAlignmentButton.clicked.connect(self.download_full_alignment)

        font = QtGui.QFont()
        font.setPointSize(12)
        font.setFamily('Courier New')
        self.TM1_align = QtWidgets.QLabel(self.centralwidget)
        self.TM1_align.setGeometry(QtCore.QRect(110, 50, 700, 51))
        self.TM1_align.setFont(font)
        self.TM1_align.setObjectName("TM1_align")
        self.TM2_align = QtWidgets.QLabel(self.centralwidget)
        self.TM2_align.setGeometry(QtCore.QRect(110, 110, 700, 51))
        self.TM2_align.setFont(font)
        self.TM2_align.setObjectName("TM2_align")
        self.TM3_align = QtWidgets.QLabel(self.centralwidget)
        self.TM3_align.setGeometry(QtCore.QRect(110, 190, 700, 51))
        self.TM3_align.setFont(font)
        self.TM3_align.setObjectName("TM3_align")
        self.TM4_align = QtWidgets.QLabel(self.centralwidget)
        self.TM4_align.setGeometry(QtCore.QRect(110, 260, 700, 51))
        self.TM4_align.setFont(font)
        self.TM4_align.setObjectName("TM4_align")
        self.TM5_align = QtWidgets.QLabel(self.centralwidget)
        self.TM5_align.setGeometry(QtCore.QRect(110, 330, 700, 51))
        self.TM5_align.setFont(font)
        self.TM5_align.setObjectName("TM5_align")
        self.TM6_align = QtWidgets.QLabel(self.centralwidget)
        self.TM6_align.setGeometry(QtCore.QRect(110, 400, 700, 51))
        self.TM6_align.setFont(font)
        self.TM6_align.setObjectName("TM6_align")
        self.TM7_align = QtWidgets.QLabel(self.centralwidget)
        self.TM7_align.setGeometry(QtCore.QRect(110, 470, 700, 51))
        self.TM7_align.setFont(font)
        self.TM7_align.setObjectName("TM7_align")
        self.target_template = QtWidgets.QLabel(self.centralwidget)
        self.target_template.setGeometry(QtCore.QRect(280, 10, 311, 41))

        font = QtGui.QFont()
        font.setPointSize(12)
        self.target_template.setFont(font)
        self.target_template.setObjectName("target_template")
        self.alignment.setCentralWidget(self.centralwidget)

        self.menubar = QtWidgets.QMenuBar(self.alignment)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 20))
        self.menubar.setObjectName("menubar")

        self.alignment.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(self.alignment)
        self.statusbar.setObjectName("statusbar")
        self.alignment.setStatusBar(self.statusbar)

        self.translate_ui()
        QtCore.QMetaObject.connectSlotsByName(self.alignment)

    def translate_ui(self):
        _translate = QtCore.QCoreApplication.translate
        self.alignment.setWindowTitle(_translate("Alignment", "MainWindow"))
        self.closeButton.setText(_translate("Alignment", "Close"))
        self.downloadTmWiseAlignmentButton.setText(_translate("Alignment", "Download TM-wise Alignment"))
        self.downloadFullAlignmentButton.setText(_translate("Alignment", "Download Full Alignment"))

        self.label.setText(_translate("Alignment", '<html><head/><body><p><span style="font-size:12pt; font-weight:600;">TM1</span></p></body></html>'))
        self.label_2.setText(_translate("Alignment", '<html><head/><body><p><span style="font-size:12pt; font-weight:600;">TM3</span></p></body></html>'))
        self.label_3.setText(_translate("Alignment", '<html><head/><body><p><span style="font-size:12pt; font-weight:600;">TM2</span></p></body></html>'))
        self.label_4.setText(_translate("Alignment", '<html><head/><body><p><span style="font-size:12pt; font-weight:600;">TM4</span></p></body></html>'))
        self.label_5.setText(_translate("Alignment", '<html><head/><body><p><span style="font-size:12pt; font-weight:600;">TM5</span></p></body></html>'))
        self.label_6.setText(_translate("Alignment", '<html><head/><body><p><span style="font-size:12pt; font-weight:600;">TM6</span></p></body></html>'))
        self.label_7.setText(_translate("Alignment", '<html><head/><body><p><span style="font-size:12pt; font-weight:600;">TM7</span></p></body></html>'))

        # calculate string sizes for padding
        id_max_size = len(self.template_receptor) if len(self.template_receptor) > len(self.target_receptor) else len(self.target_receptor)
        align_max_size = 0
        for i in range(7):
            if align_max_size < len(self.tm_align[0][i]):
                align_max_size = len(self.tm_align[0][i])

        self.TM1_align.setText(_translate("Alignment", self.format_string(0, id_max_size, align_max_size)))
        self.TM2_align.setText(_translate("Alignment", self.format_string(1, id_max_size, align_max_size)))
        self.TM3_align.setText(_translate("Alignment", self.format_string(2, id_max_size, align_max_size)))
        self.TM4_align.setText(_translate("Alignment", self.format_string(3, id_max_size, align_max_size)))
        self.TM5_align.setText(_translate("Alignment", self.format_string(4, id_max_size, align_max_size)))
        self.TM6_align.setText(_translate("Alignment", self.format_string(5, id_max_size, align_max_size)))
        self.TM7_align.setText(_translate("Alignment", self.format_string(6, id_max_size, align_max_size)))

        self.target_template.setText(_translate("Alignment", self.template_receptor + "-" + self.target_receptor))

    def format_string(self, i, id_max_size, align_max_size):
        return self.template_receptor.ljust(id_max_size) + '\t' + str(self.template_number[i][0]) + '\t' + self.tm_align[0][i].ljust(align_max_size) + '\t' + str(self.template_number[i][1]) + '\n' \
               + self.target_receptor.ljust(id_max_size) + '\t' + str(self.target_number[i][0]) + '\t' + self.tm_align[1][i].ljust(align_max_size) + '\t' + str(self.target_number[i][1])

    def hide_window(self):
        self.alignment.hide()

    def download_tm_wise_alignment(self):
        filename = 'tm_wise_alignment.txt'
        with open(filename, 'w') as f:
            for i in range(7):
                f.write(self.tm_align[0][i] + '\n' + self.tm_align[1][i] + '\n\n')

        QtWidgets.QMessageBox.about(QtWidgets.QMessageBox(), 'Information', 'Downloaded file ' + filename)

    def download_full_alignment(self):
        alignment_text = Helper.download_full_alignment(self.tm_align, self.template_receptor, self.target_receptor, self.template_sequence, self.target_sequence, self.template_number, self.target_number)
        print(alignment_text)
        # print to file
        filename = 'full_alignment.txt'
        with open(filename, 'w') as f:
            f.write(alignment_text)

        QtWidgets.QMessageBox.about(QtWidgets.QMessageBox(), 'Information', 'Downloaded file ' + filename)
