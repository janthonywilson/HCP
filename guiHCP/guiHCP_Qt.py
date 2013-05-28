# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\Tony\Documents\Qt\guiHCP.ui'
#
# Created: Mon May 20 22:45:20 2013
#      by: PyQt4 UI code generator 4.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(983, 502)
        self.groupBox = QtGui.QGroupBox(Dialog)
        self.groupBox.setGeometry(QtCore.QRect(10, 20, 201, 171))
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.UserNameLineEdit = QtGui.QLineEdit(self.groupBox)
        self.UserNameLineEdit.setGeometry(QtCore.QRect(10, 20, 181, 20))
        self.UserNameLineEdit.setText(_fromUtf8(""))
        self.UserNameLineEdit.setObjectName(_fromUtf8("UserNameLineEdit"))
        self.PasswordLineEdit = QtGui.QLineEdit(self.groupBox)
        self.PasswordLineEdit.setGeometry(QtCore.QRect(10, 50, 181, 20))
        self.PasswordLineEdit.setText(_fromUtf8(""))
        self.PasswordLineEdit.setEchoMode(QtGui.QLineEdit.Password)
        self.PasswordLineEdit.setObjectName(_fromUtf8("PasswordLineEdit"))
        self.comboBox = QtGui.QComboBox(self.groupBox)
        self.comboBox.setGeometry(QtCore.QRect(10, 80, 181, 22))
        self.comboBox.setObjectName(_fromUtf8("comboBox"))
        self.comboBox.addItem(_fromUtf8(""))
        self.comboBox.addItem(_fromUtf8(""))
        self.comboBox.addItem(_fromUtf8(""))
        self.pushOK = QtGui.QPushButton(self.groupBox)
        self.pushOK.setGeometry(QtCore.QRect(10, 110, 81, 23))
        self.pushOK.setObjectName(_fromUtf8("pushOK"))
        self.pushCancel = QtGui.QPushButton(self.groupBox)
        self.pushCancel.setGeometry(QtCore.QRect(110, 110, 81, 23))
        self.pushCancel.setObjectName(_fromUtf8("pushCancel"))
        self.StatusLineEdit = QtGui.QLineEdit(self.groupBox)
        self.StatusLineEdit.setGeometry(QtCore.QRect(10, 140, 181, 20))
        self.StatusLineEdit.setObjectName(_fromUtf8("StatusLineEdit"))
        self.groupBox_2 = QtGui.QGroupBox(Dialog)
        self.groupBox_2.setGeometry(QtCore.QRect(220, 20, 181, 171))
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.subjectsButton = QtGui.QPushButton(self.groupBox_2)
        self.subjectsButton.setGeometry(QtCore.QRect(10, 140, 161, 23))
        self.subjectsButton.setObjectName(_fromUtf8("subjectsButton"))
        self.projectListWidget = QtGui.QListWidget(self.groupBox_2)
        self.projectListWidget.setGeometry(QtCore.QRect(10, 20, 161, 111))
        self.projectListWidget.setObjectName(_fromUtf8("projectListWidget"))
        self.groupBox_3 = QtGui.QGroupBox(Dialog)
        self.groupBox_3.setGeometry(QtCore.QRect(410, 20, 181, 171))
        self.groupBox_3.setObjectName(_fromUtf8("groupBox_3"))
        self.sessionsButton = QtGui.QPushButton(self.groupBox_3)
        self.sessionsButton.setGeometry(QtCore.QRect(10, 140, 161, 23))
        self.sessionsButton.setObjectName(_fromUtf8("sessionsButton"))
        self.subjectsListWidget = QtGui.QListWidget(self.groupBox_3)
        self.subjectsListWidget.setGeometry(QtCore.QRect(10, 20, 161, 111))
        self.subjectsListWidget.setObjectName(_fromUtf8("subjectsListWidget"))
        self.groupBox_4 = QtGui.QGroupBox(Dialog)
        self.groupBox_4.setGeometry(QtCore.QRect(600, 20, 181, 171))
        self.groupBox_4.setObjectName(_fromUtf8("groupBox_4"))
        self.getSessionMetaButton = QtGui.QPushButton(self.groupBox_4)
        self.getSessionMetaButton.setGeometry(QtCore.QRect(10, 140, 161, 23))
        self.getSessionMetaButton.setObjectName(_fromUtf8("getSessionMetaButton"))
        self.sessionsListWidget = QtGui.QListWidget(self.groupBox_4)
        self.sessionsListWidget.setGeometry(QtCore.QRect(10, 20, 161, 81))
        self.sessionsListWidget.setObjectName(_fromUtf8("sessionsListWidget"))
        self.getSessionResourcesButton = QtGui.QPushButton(self.groupBox_4)
        self.getSessionResourcesButton.setGeometry(QtCore.QRect(10, 110, 161, 23))
        self.getSessionResourcesButton.setObjectName(_fromUtf8("getSessionResourcesButton"))
        self.sessionMetaTableWidget = QtGui.QTableWidget(Dialog)
        self.sessionMetaTableWidget.setGeometry(QtCore.QRect(10, 200, 961, 291))
        self.sessionMetaTableWidget.viewport().setProperty("cursor", QtGui.QCursor(QtCore.Qt.ArrowCursor))
        self.sessionMetaTableWidget.setObjectName(_fromUtf8("sessionMetaTableWidget"))
        self.sessionMetaTableWidget.setColumnCount(3)
        self.sessionMetaTableWidget.setRowCount(0)
        item = QtGui.QTableWidgetItem()
        self.sessionMetaTableWidget.setHorizontalHeaderItem(0, item)
        item = QtGui.QTableWidgetItem()
        self.sessionMetaTableWidget.setHorizontalHeaderItem(1, item)
        item = QtGui.QTableWidgetItem()
        self.sessionMetaTableWidget.setHorizontalHeaderItem(2, item)
        self.sessionMetaTableWidget.horizontalHeader().setDefaultSectionSize(302)
        self.groupBox_5 = QtGui.QGroupBox(Dialog)
        self.groupBox_5.setGeometry(QtCore.QRect(790, 20, 181, 171))
        self.groupBox_5.setObjectName(_fromUtf8("groupBox_5"))
        self.resourcesListWidget = QtGui.QListWidget(self.groupBox_5)
        self.resourcesListWidget.setGeometry(QtCore.QRect(10, 20, 161, 111))
        self.resourcesListWidget.setObjectName(_fromUtf8("resourcesListWidget"))
        self.getResourcesMetaButton = QtGui.QPushButton(self.groupBox_5)
        self.getResourcesMetaButton.setGeometry(QtCore.QRect(10, 140, 161, 23))
        self.getResourcesMetaButton.setObjectName(_fromUtf8("getResourcesMetaButton"))

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.groupBox.setTitle(_translate("Dialog", "pyHCP Initializer", None))
        self.UserNameLineEdit.setToolTip(_translate("Dialog", "User Name", None))
        self.UserNameLineEdit.setStatusTip(_translate("Dialog", "User Name", None))
        self.PasswordLineEdit.setToolTip(_translate("Dialog", "Password", None))
        self.PasswordLineEdit.setStatusTip(_translate("Dialog", "Password", None))
        self.comboBox.setToolTip(_translate("Dialog", "Server", None))
        self.comboBox.setItemText(0, _translate("Dialog", "db.humanconnectome.org", None))
        self.comboBox.setItemText(1, _translate("Dialog", "intradb.humanconnectome.org", None))
        self.comboBox.setItemText(2, _translate("Dialog", "demo.humanconnectome.org", None))
        self.pushOK.setText(_translate("Dialog", "OK", None))
        self.pushCancel.setText(_translate("Dialog", "Cancel", None))
        self.groupBox_2.setTitle(_translate("Dialog", "Projects", None))
        self.subjectsButton.setText(_translate("Dialog", "Get Subjects", None))
        self.groupBox_3.setTitle(_translate("Dialog", "Subjects", None))
        self.sessionsButton.setText(_translate("Dialog", "Get Sessions", None))
        self.groupBox_4.setTitle(_translate("Dialog", "Sessions", None))
        self.getSessionMetaButton.setText(_translate("Dialog", "Get Session Metadata", None))
        self.getSessionResourcesButton.setText(_translate("Dialog", "Get Session Resources", None))
        item = self.sessionMetaTableWidget.horizontalHeaderItem(0)
        item.setText(_translate("Dialog", "Types", None))
        item = self.sessionMetaTableWidget.horizontalHeaderItem(1)
        item.setText(_translate("Dialog", "Series", None))
        item = self.sessionMetaTableWidget.horizontalHeaderItem(2)
        item.setText(_translate("Dialog", "Quality", None))
        self.groupBox_5.setTitle(_translate("Dialog", "Resources", None))
        self.getResourcesMetaButton.setText(_translate("Dialog", "Get Resources Meta", None))

