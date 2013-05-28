'''
Created on May 2, 2013

@author: Tony
'''
#===============================================================================
# C:\Python27\lib\site-packages\PyQt4>pyuic4.bat C:\Users\Tony\Documents\Qt\guiHCP.ui > C:\Users\Tony\workspace\guiHCP\guiHCP_Qt.py
#===============================================================================

import sys
import guiHCP_Qt
from pyHCP import pyHCP, getHCP
from PyQt4 import QtGui, QtCore
        
class StartHCP(QtGui.QDialog, guiHCP_Qt.Ui_Dialog):
    def __init__(self,parent=None):
        QtGui.QDialog.__init__(self,parent)
        self.setupUi(self)
#        super(pyHCP, self).__init__()
#        super(getHCP, self).__init__()
    #===============================================================================   
    def getValues(self):
        return self.isHidden()
    #===============================================================================
    def getConnect(self):
        global pyHCP, getHCP
#        text, ok = QtGui.QInputDialog.getText(self, 'Input Dialog', 'Enter your name:')
        pyHCP = pyHCP(str(self.UserNameLineEdit.text()), str(self.PasswordLineEdit.text()), str(self.comboBox.currentText()))
        getHCP = getHCP(pyHCP)
        Projects = getHCP.getProjects().get('Projects')
        
        self.StatusLineEdit.setText('Session ID: %s' % getHCP.SessionId)
#        self.textEdit.setText('\n'.join(Projects[0]))
        
#        self.projectsComboBox.addItems(Projects)
        
        self.projectListWidget.clear()
        for i in xrange(0, len(Projects)):
            self.projectListWidget.addItem(Projects[i])
    #===============================================================================    
    def getSubjects(self):
        global getHCP

#        self.subjectsComboBox.clear()
        
#        print self.projectListWidget.currentItem().text()
        getHCP.Project = str(self.projectListWidget.currentItem().text())
#        getHCP.Project = str(self.projectsComboBox.currentText())
        Subjects = sorted(getHCP.getSubjects())

#        self.subjectsTextEdit.setText('\n'.join(Subjects))
#        self.subjectsComboBox.addItems(Subjects)
        
        self.subjectsListWidget.clear()
        for i in xrange(0, len(Subjects)):
            self.subjectsListWidget.addItem(Subjects[i])
    #===============================================================================    
    def getSubjectSessions(self):
        global getHCP

#        self.sessionsComboBox.clear()
        getHCP.Subject = str(self.subjectsListWidget.currentItem().text())
#        getHCP.Subject = str(self.subjectsComboBox.currentText())
        Sessions = getHCP.getSubjectSessions()
        SessionsList = list(set(Sessions.get('Sessions')))
#        self.sessionsTextEdit.setText('\n'.join(Sessions.get('Sessions')))
#        self.sessionsComboBox.addItems(SessionsList)
        
        
        self.sessionsListWidget.clear()
        for i in xrange(0, len(SessionsList)):
            self.sessionsListWidget.addItem(SessionsList[i])
    #===============================================================================
    def getSessionResources(self):
        global getHCP

#        self.sessionsComboBox.clear()
        getHCP.Session = str(self.sessionsListWidget.currentItem().text())
#        getHCP.Subject = str(self.subjectsComboBox.currentText())
        Resources = getHCP.getSubjectResources()
        ResourcesList = list(set(Resources.get('Names')))
#        self.sessionsTextEdit.setText('\n'.join(Sessions.get('Sessions')))
#        self.sessionsComboBox.addItems(SessionsList)
        
        
        self.resourcesListWidget.clear()
        for i in xrange(0, len(ResourcesList)):
            self.resourcesListWidget.addItem(ResourcesList[i])
    #===============================================================================        
    def getSessionResourceMeta(self):
        global getHCP
        
        getHCP.Resource = str(self.resourcesListWidget.currentItem().text())
        ResourceMeta = getHCP.getSubjectResourceMeta()
        
        # 'Name', 'Bytes', 'URI', 'Path', 'Readable', 'RealPath', 'Label', 'Format', 'Contents'
        Names = ResourceMeta.get('Name')
        Bytes = ResourceMeta.get('Bytes')
        URIs = ResourceMeta.get('URI')
        
        self.sessionMetaTableWidget.clearSpans()
        self.sessionMetaTableWidget.setRowCount(0)
        
        self.sessionMetaTableWidget.setHorizontalHeaderItem(0, QtGui.QTableWidgetItem('File Name'))
        self.sessionMetaTableWidget.setHorizontalHeaderItem(1, QtGui.QTableWidgetItem('Bytes'))
        self.sessionMetaTableWidget.setHorizontalHeaderItem(2, QtGui.QTableWidgetItem('URI'))
        for i in xrange(0, len(Names)):

            self.sessionMetaTableWidget.removeRow(i)
            self.sessionMetaTableWidget.insertRow(i)
            self.sessionMetaTableWidget.setItem(i, 0, QtGui.QTableWidgetItem(Names[i]))
            self.sessionMetaTableWidget.setItem(i, 1, QtGui.QTableWidgetItem(Bytes[i]))
            self.sessionMetaTableWidget.setItem(i, 2, QtGui.QTableWidgetItem(URIs[i]))
    #===============================================================================
    def getSessionMeta(self):
#        SessionMeta = {'IDs':ScanIds, 'Types':ScanTypes, 'Series':ScanSeries, 'Quality':ScanQualty, 'XNATID':ScanXnatId }
        global getHCP

        getHCP.Session = str(self.sessionsListWidget.currentItem().text())
#        getHCP.Session = str(self.sessionsComboBox.currentText())
        SessionMeta = getHCP.getSessionMeta()
#        self.typeTextEdit.setText('\n'.join(SessionMeta.get('Types')))
#        self.seriesTextEdit.setText('\n'.join(SessionMeta.get('Series')))
#        self.qualityTextEdit.setText('\n'.join(SessionMeta.get('Quality')))
        
        Types = SessionMeta.get('Types')
        Series = SessionMeta.get('Series')
        Quality = SessionMeta.get('Quality')
        
#        self.sessionMetaTableWidget.clearContents()
#        self.sessionMetaTableWidget.clear()
        self.sessionMetaTableWidget.clearSpans()
        self.sessionMetaTableWidget.setRowCount(0)
        

        self.sessionMetaTableWidget.setHorizontalHeaderItem(0, QtGui.QTableWidgetItem('Types'))
        self.sessionMetaTableWidget.setHorizontalHeaderItem(1, QtGui.QTableWidgetItem('Series'))
        self.sessionMetaTableWidget.setHorizontalHeaderItem(2, QtGui.QTableWidgetItem('Quality'))
        for i in xrange(0, len(Types)):

            self.sessionMetaTableWidget.removeRow(i)
            self.sessionMetaTableWidget.insertRow(i)
            self.sessionMetaTableWidget.setItem(i, 0, QtGui.QTableWidgetItem(Types[i]))
            self.sessionMetaTableWidget.setItem(i, 1, QtGui.QTableWidgetItem(Series[i]))
            self.sessionMetaTableWidget.setItem(i, 2, QtGui.QTableWidgetItem(Quality[i]))
    #===============================================================================    
    def cancelConnect(self):
        exitAction = QtGui.QAction(QtGui.QIcon('exit24.png'), 'Exit', self)
#        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)

if __name__ == '__main__':
    
    app = QtGui.QApplication(sys.argv)
    
    dlg = StartHCP()

#    print dlg.getValues()
#    print dlg.getConnect()
    app.connect(dlg.pushOK, QtCore.SIGNAL('clicked()'), dlg.getConnect)
    app.connect(dlg.pushCancel, QtCore.SIGNAL('clicked()'), app.quit)
#    app.connect(dlg.pushCancel, QtCore.SIGNAL('clicked()'), dlg.cancelConnect)
    app.connect(dlg.subjectsButton, QtCore.SIGNAL('clicked()'), dlg.getSubjects)
    app.connect(dlg.sessionsButton, QtCore.SIGNAL('clicked()'), dlg.getSubjectSessions)
    app.connect(dlg.getSessionMetaButton, QtCore.SIGNAL('clicked()'), dlg.getSessionMeta)
    app.connect(dlg.getSessionResourcesButton, QtCore.SIGNAL('clicked()'), dlg.getSessionResources)
    app.connect(dlg.getResourcesMetaButton, QtCore.SIGNAL('clicked()'), dlg.getSessionResourceMeta)
    
    dlg.show()
    dlg.exec_()
    

    sys.exit(app.exec_())

#    app = QtGui.QApplication(sys.argv)
#    ex = guiHCP()
#    sys.exit(app.exec_())




















#===============================================================================
# class guiHCP(QtGui.QMainWindow):
#    
#    def __init__(self):
#        super(guiHCP, self).__init__()
#        
#        self.initUI()
#        
#    def initUI(self):               
#        
#        textEdit = QtGui.QTextEdit()
#        self.setCentralWidget(textEdit)
# 
#        exitAction = QtGui.QAction(QtGui.QIcon('exit24.png'), 'Exit', self)
#        exitAction.setShortcut('Ctrl+Q')
#        exitAction.setStatusTip('Exit application')
#        exitAction.triggered.connect(self.close)
# 
#        self.statusBar()
# 
#        menubar = self.menuBar()
#        fileMenu = menubar.addMenu('&File')
#        fileMenu.addAction(exitAction)
# 
#        toolbar = self.addToolBar('Exit')
#        toolbar.addAction(exitAction)
#        
#        self.setGeometry(300, 300, 350, 250)
#        self.setWindowTitle('Main window')    
#        self.show()
#===============================================================================
