# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'aPdb.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_aPdb(object):
    def __init__(self, aPdb, df4):
        self.aPdb = aPdb
        self.df_pdb = df4

    def setup_ui(self):
        self.aPdb.setObjectName("aPdb")
        self.aPdb.resize(800, 600)

        self.centralwidget = QtWidgets.QWidget(self.aPdb)
        self.centralwidget.setObjectName("centralwidget")

        self.table_widget = QtWidgets.QTableWidget(self.centralwidget)
        self.table_widget.setGeometry(QtCore.QRect(210, 100, 320, 411))
        self.table_widget.setObjectName("tableWidget")
        self.table_widget.setRowCount(0)
        self.table_widget.setColumnCount(3)
        for i in range(3):
            self.table_widget.setHorizontalHeaderItem(i, QtWidgets.QTableWidgetItem())

        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(270, 10, 221, 61))

        font = QtGui.QFont()
        font.setPointSize(13)
        self.label.setFont(font)
        self.label.setObjectName("label")

        self.home_button = QtWidgets.QPushButton(self.centralwidget)
        self.home_button.setGeometry(QtCore.QRect(700, 10, 80, 23))
        self.home_button.setObjectName("homeButton")
        self.home_button.clicked.connect(self.hide_window)

        self.aPdb.setCentralWidget(self.centralwidget)
        self.table_widget.insertRow(self.table_widget.rowCount())

        self.statusbar = QtWidgets.QStatusBar(self.aPdb)
        self.statusbar.setObjectName("statusbar")
        self.aPdb.setStatusBar(self.statusbar)

        self.translate_ui()
        QtCore.QMetaObject.connectSlotsByName(self.aPdb)

    def translate_ui(self):
        _translate = QtCore.QCoreApplication.translate
        self.aPdb.setWindowTitle(_translate("aPdb", "MainWindow"))
        item = self.table_widget.horizontalHeaderItem(0)
        item.setText(_translate("aPdb", "PDBID"))
        item = self.table_widget.horizontalHeaderItem(1)
        item.setText(_translate("aPdb", "Resolution"))
        item = self.table_widget.horizontalHeaderItem(2)
        item.setText(_translate("aPdb", "Position"))
        self.label.setText(_translate("aPdb", "<html><head/><body><p><span style=\" font-weight:600;\">Available PDBs</span></p></body></html>"))
        self.home_button.setText(_translate("aPdb", "Close"))

        rows = self.df_pdb.shape[0]
        self.table_widget.setRowCount(rows)

        for i in range(rows):
            pdbid = self.df_pdb['PDBID'].iloc[i]
            res = self.df_pdb['Resolution'].iloc[i]
            pos = self.df_pdb['Position'].iloc[i]

            self.table_widget.setItem(i, 0, QtWidgets.QTableWidgetItem(pdbid))
            self.table_widget.setItem(i, 1, QtWidgets.QTableWidgetItem(str(res)))
            self.table_widget.setItem(i, 2, QtWidgets.QTableWidgetItem(pos))

    def hide_window(self):
        self.aPdb.hide()
