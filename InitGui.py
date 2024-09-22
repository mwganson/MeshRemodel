# -*- coding: utf-8 -*-
###################################################################################
#
#  InitGui.py
#  
#  Copyright 2019 Mark Ganson <TheMarkster> mwganson at gmail
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
###################################################################################

import meshremodelwb_locator
meshremodelWBPath = os.path.dirname(meshremodelwb_locator.__file__)
meshremodelWB_icons_path = os.path.join(meshremodelWBPath,'Resources','icons')

global main_meshremodelWB_Icon

main_meshremodelWB_Icon = os.path.join(meshremodelWB_icons_path , 'CreatePointsObject.svg')

#def myFunc(string):
#    print (string)
#    global act
#    act.setVisible(True)

#mw=Gui.getMainWindow()
#bar=mw.menuBar()
#act=bar.addAction("MyCmd")
#mw.workbenchActivated.connect(myFunc)

####################################################################################
# Initialize the workbench 
class MeshRemodelWorkbench(Workbench):
 

    global main_meshremodelWB_Icon

    MenuText = "Mesh Remodel"
    ToolTip = "MeshRemodel workbench"
    Icon = main_meshremodelWB_Icon #defined in package.xml
    
    def __init__(self):
        pass

    def Initialize(self):
        "This function is executed when FreeCAD starts"
        import MeshRemodelCmd #needed files for FreeCAD commands
        self.list = [
                    "MeshRemodelGroupCommandPointsObjects",
                  #  "MeshRemodelCreatePointsObject",
                  #  "MeshRemodelCreateWireFrameObject",
                    "MeshRemodelMeshBoundaryWires",
                    "MeshRemodelCreateCrossSectionsObject",
                    "MeshRemodelMovePoint",
                    "MeshRemodelRemovePoint",
                    "MeshRemodelAddOrRemoveFacet",
                    "MeshRemodelExpandedMesh",
                    "MeshRemodelMeshSimpleCopy",
                    "MeshRemodelOffsetMesh",
                    "MeshRemodelDuplicateSelectedFacets",
                  #  "MeshRemodelCreateCoplanarPointsObject",
                    "MeshRemodelAddSelectionObserver",
                    "MeshRemodelPartSolid",
                    "MeshRemodelSubObjectLoft",
                    "MeshRemodelCreateGridSurface",
                    "MeshRemodelCreatePointObject",
                    "MeshRemodelCreateLine",
                    "MeshRemodelCreatePolygon",
                    "MeshRemodelCreateBSpline",
                    "MeshRemodelCreatePlane",
                    "MeshRemodelCreateCircle",
                    "MeshRemodelCreateArc",
                    "MeshRemodelCreateWire",
                    "MeshRemodelRotateObject",
                    "MeshRemodelMoveAxial",
                    "MeshRemodelGoBackSelection",
                    "MeshRemodelGroupCommandExtras",
                    # "MeshRemodelMergeSketches",
                    # "MeshRemodelValidateSketch",
                    # "MeshRemodelPartCheckGeometry",
                    # "MeshRemodelSubShapeBinder",
                    # "MeshRemodelSettings",
                    
                    ] # A list of command names created in the line above
        self.appendToolbar("MeshRemodel Commands",self.list) # leave settings, validate sketch and merge sketch off toolbar
        self.appendMenu("Mesh&Remodel",self.list) # creates a new menu
        #considered putting the menu inside the Edit menu, but decided against it
        #self.appendMenu(["&Edit","MeshRemodel"],self.list) # appends a submenu to an existing menu

 
    def Activated(self):
        "This function is executed when the workbench is activated"
        import requests
        import xml.etree.ElementTree as ET
        
        def get_remote_version(user, repo, branch='master'):
            # GitHub raw URL for package.xml
            url = f"https://raw.githubusercontent.com/{user}/{repo}/{branch}/package.xml"
            try:
                response = requests.get(url)
                if response.status_code == 200:
                    # Parse the XML content
                    xml_content = response.content
                    tree = ET.ElementTree(ET.fromstring(xml_content))
                    root = tree.getroot()
                    # Find the version element and return its text
                    version = root.find("version").text
                    return version
                else:
                    print(f"Failed to fetch package.xml: {response.status_code}")
                    return None
            except Exception as e:
                print(f"Error fetching or parsing package.xml: {e}")
                return None
        
        def check_for_update(current_version, user, repo, branch, callback):
            latest_version = get_remote_version(user, repo, branch)
            if latest_version and latest_version != current_version:
                callback(latest_version)
        
        # Example usage
        def update_callback(latest_version):
            FreeCAD.Console.PrintWarning(f"MeshRemodel {latest_version} is now available in the Addon Manager.\n")
        
        import MeshRemodelCmd
        current_version = MeshRemodelCmd.__version__
        user = "mwganson"
        repo = "MeshRemodel"
        branch = "master"
        
        # Check for updates
        pg = FreeCAD.ParamGet("User parameter:/Plugins/MeshRemodel")
        checkUpdates = pg.GetBool("CheckForUpdates", True)
        #print(f"checkUpdates = {checkUpdates}")
        if checkUpdates:
            check_for_update(current_version, user, repo, branch, update_callback)

        return
 
    def Deactivated(self):
        "This function is executed when the workbench is deactivated"

        #FreeCAD will hide our menu and toolbar upon exiting the wb, so we setup a singleshot
        #to unhide them once FreeCAD is finished, 2 seconds later
        from PySide import QtCore
        QtCore.QTimer.singleShot(2000, self.showMenu)
        return 
        
    def showMenu(self):
        from PySide import QtGui
        window = QtGui.QApplication.activeWindow()
        #freecad hides wb toolbars on leaving wb, we unhide ours here to keep it around
        #if the user has it set in parameters to do so
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        keep = pg.GetBool('KeepToolbar',True)
        if not keep:
            return
        tb = window.findChildren(QtGui.QToolBar) if window else []
        for bar in tb:
            if "MeshRemodel Commands" in bar.objectName():
                bar.setVisible(True)

    def ContextMenu(self, recipient):
        "This is executed whenever the user right-clicks on screen"
        # "recipient" will be either "view" or "tree"
        self.appendContextMenu("MeshRemodel",self.list) # add commands to the context menu
 
    def GetClassName(self): 
        # this function is mandatory if this is a full python workbench
        return "Gui::PythonWorkbench"
wb = MeshRemodelWorkbench()
Gui.addWorkbench(wb)





