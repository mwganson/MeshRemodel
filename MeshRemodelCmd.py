# -*- coding: utf-8 -*-
###################################################################################
#
#  MeshRemodelCmd.py
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

__title__   = "MeshRemodel"
__author__  = "Mark Ganson <TheMarkster>"
__url__     = "https://github.com/mwganson/MeshRemodel"
__date__    = "2020.08.30"
__version__ = "1.72"
version = 1.72

import FreeCAD, FreeCADGui, Part, os, math
from PySide import QtCore, QtGui
import Draft, DraftGeomUtils, DraftVecUtils


if FreeCAD.GuiUp:
    from FreeCAD import Gui

__dir__ = os.path.dirname(__file__)
iconPath = os.path.join( __dir__, 'Resources', 'icons' )
keepToolbar = False
windowFlags = QtCore.Qt.WindowTitleHint | QtCore.Qt.WindowCloseButtonHint #no ? in title bar
global_picked = [] #picked points list for use with selection by preselection observer

######################################################################################
# geometry utilities

class MeshRemodelGeomUtils(object):
    """Geometry Utilities"""

#source for this block of code: https://stackoverflow.com/questions/9866452/calculate-volume-of-any-tetrahedron-given-4-points
#4 points are coplanar if the tetrahedron defined by them has volume = 0
##################################################################
    def determinant_3x3(self,m):
        """helper for isCoplanar()"""
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                m[1][0] * (m[0][1] * m[2][2] - m[0][2] * m[2][1]) +
                m[2][0] * (m[0][1] * m[1][2] - m[0][2] * m[1][1]))


    def subtract(self, a, b):
        """ helper for isCoplanar()"""
        return (a[0] - b[0],
                a[1] - b[1],
                a[2] - b[2])

    def tetrahedron_calc_volume(self, a, b, c, d):
        """helper for isCoplanar()"""
        return (abs(self.determinant_3x3((self.subtract(a, b),
                                 self.subtract(b, c),
                                 self.subtract(c, d),
                                 ))) / 6.0)

#a = [0.0, 0.0, 0.0]
#d = [2.0, 0.0, 0.0]
#c = [0.0, 2.0, 0.0]
#b = [0.0, 0.0, 2.0]

#print(tetrahedron_calc_volume(a, b, c, d))
    def isCoplanar(self, trio, pt, tol=1e-3):
        """ isCoplanar(trio, pt, tol=1e-3)
            trio is a 3-element list of vectors, pt is a vector to test, tol is tolerance, return True if all 4 are coplanar
            test is done by creating a tetrahedron and testing its volume against tol
            a tetrahedron from 4 coplanar points should have volume ~= 0
        """
        if len(trio) != 3:
            raise Exception("MeshRemodel GeomUtils Error: isCoplanar() trio parameter must be list of 3 vectors")
        A,B,C,D = trio[0],trio[1],trio[2],pt
        vol = self.tetrahedron_calc_volume(A,B,C,D)

        if vol <= tol:
            return True
        return False

    def hasPoint(self,pt,lis,tol):
        """hasPoint(pt,lis,tol)"""
        for l in lis:
            if self.isSamePoint(pt,l,tol):
                return True
        return False

    def isSamePoint(self,A,B,tol):
        """isSamePoint(A,B,tol)"""
        dis = self.dist(A,B)
        if dis < tol:
            return True
        return False

    def midpoint(self, A, B):
        """ midpoint(A, B)
            A,B are vectors, return midpoint"""
        mid = FreeCAD.Base.Vector()
        mid.x = (A.x + B.x)/2.0
        mid.y = (A.y + B.y)/2.0
        mid.z = (A.z + B.z)/2.0
        return mid

    def dist(self, p1, p2):
        """ dist (p1, p2)
            3d distance between vectors p1 and p2"""
        return self.getDistance3d(p1[0],p1[1],p1[2],p2[0],p2[1],p2[2])

    def getDistance3d(self, x1, y1, z1, x2, y2, z2):
        """ getDistance3d(x1, y1, z1, x2, y2, z2)
            3d distance between x1,y1,z1 and x2,y2,z2 float parameters"""
        return math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)

    def sortPoints(self,pts):
        """ sortPoints(pts)
            sort pts, a list of vectors, according to distance from one point to the next
            pts[0] is taken first, then the nearest point to it is placed at pts[1]
            then pts[1] is used to find the nearest point to it and placed at pts[2], and so on
        """
        newList = [pts[0]]
        for ii in range(0, len(pts)):
            newList.extend([self.nearestPoint(newList[ii],pts,newList)])
        return newList

    def nearestPoint(self, pt, pts, exclude):
        """ nearestPoint(pt, pts, exclude)
            pt is a vector, pts a list of vectors
            exlude is a list of vectors to exclude from process
            return nearest point to pt in pts and not in exclude"""
        if len(pts) == 0: #should never happen
            raise Exception("MeshRemodel GeomUtils Error: nearestPoint() pts length = 0\n")
        nearest = pts[0]
        d = 10**100
        for p in pts:
            if p in exclude:
                continue
            dis = self.dist(pt, p)
            if dis < d:
                d = dis
                nearest = p
        return nearest

    def isColinear(self,A,B,C):
        """ isColinear(A, B, C)
            determine whether vectors A,B,C are colinear """
        return DraftVecUtils.isColinear([A,B,C])

    def incenter(A,B,C):
        """ incenter(A, B, C)
            return incenter (vector) of triangle at vectors A,B,C 
            incenter is center of circle fitting inside the triangle
            tangent to all 3 sides
            raises exception if A,B,C are colinear
        """

        if self.isColinear(A,B,C):
            raise Exception("MeshRemodel Error: incenter() A,B,C are colinear")

        Ax,Ay,Az = A[0],A[1],A[2] 
        Bx,By,Bz = B[0],B[1],B[2]
        Cx,Cy,Cz = C[0],C[1],C[2]

        a = self.dist(B,C)
        b = self.dist(C,A)
        c = self.dist(A,B)
        s = a+b+c

        Ix = (a*Ax+b*Bx+c*Cx)/s
        Iy = (a*Ay+b*By+c*Cy)/s
        Iz = (a*Az+b*Bz+c*Cz)/s
        I = FreeCAD.Base.Vector(Ix,Iy,Iz)
        return I

    def inradius(A,B,C):
        """ inradius(A, B, C)
            return inradius of triangle A,B,C 
            this is radius of incircle, the circle that
            fits inside the triangle, tangent to all 3 sides
        """
        return self.dist(A, self.incenter(A,B,C))

#python code below was adapted from this javascript code
#from here: https://gamedev.stackexchange.com/questions/60630/how-do-i-find-the-circumcenter-of-a-triangle-in-3d
#in a question answered by user greenthings

#function circumcenter(A,B,C) {
#    var z = crossprod(subv(C,B),subv(A,B));
#    var a=vlen(subv(A,B)),b=vlen(subv(B,C)),c=vlen(subv(C,A));
#    var r = ((b*b + c*c - a*a)/(2*b*c)) * outeradius(a,b,c);
#    return addv(midpoint(A,B),multv(normaliz(crossprod(subv(A,B),z)),r));
#}

#function outeradius(a,b,c) { /// 3 lens
#    return (a*b*c) / (4*sss_area(a,b,c));
#}

#function sss_area(a,b,c) {
#    var sp = (a+b+c)*0.5;
#    return Math.sqrt(sp*(sp-a)*(sp-b)*(sp-c));
#    //semi perimeter
#}

    def circumcenter(self,A,B,C):
        """ circumcenter(A, B, C)
            return the circumcenter of triangle A,B,C
            the circumcenter is the circle that passes through
            all 3 of the triangle's vertices
            raises exception if A,B,C are colinear
        """
        if self.isColinear(A,B,C):
            raise Exception("MeshRemodel Error: circumcenter() A,B,C are colinear")

        z = C.sub(B).cross(A.sub(B))
        a = A.sub(B).Length
        b = B.sub(C).Length
        c = C.sub(A).Length
        r = ((b*b + c*c - a*a)/(2*b*c)) * self.outerradius(a,b,c)
        return  A.sub(B).cross(z).normalize().multiply(r).add(self.midpoint(A,B))

    def outerradius(self, a, b, c):
        """ helper for circumcenter()"""
        return (a*b*c) / (4*self.sss_area(a,b,c))

    def sss_area(self,a,b,c): #semiperimeter
        """ helper for circumcenter()"""
        sp = (a+b+c)*0.5;
        return math.sqrt(sp*(sp-a)*(sp-b)*(sp-c))

    def circumradius(self, A,B,C):
        """ circumradius(A, B, C)
            returns the radius of circumcircle of triangle A,B,C
            A,B,C are vectors
            the circumcircle is the circle passing through A, B, and C
        """
        return self.dist(A, self.circumcenter(A,B,C))



gu = MeshRemodelGeomUtils()
#######################################################################################
# Settings

class MeshRemodelSettingsCommandClass(object):
    """Settings"""

    def __init__(self):
        pass

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'Settings.png') , # the name of an icon file available in the resources
            'MenuText': "&Settings" ,
            'ToolTip' : "Workbench settings dialog"}
 
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        from PySide import QtGui
        window = QtGui.QApplication.activeWindow()
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        keep = pg.GetBool('KeepToolbar',False)
        point_size = pg.GetFloat("PointSize", 4.0)
        line_width = pg.GetFloat("LineWidth", 5.0)
        prec = pg.GetInt("SketchRadiusPrecision", 1)
        coplanar_tol = pg.GetFloat("CoplanarTolerance",.001)
        items=[("","*")[keep]+"Keep the toolbar active",
            ("","*")[not keep]+"Do not keep the toolbar active",
            "Change point size ("+str(point_size)+")",
            "Change line width ("+str(line_width)+")",
            "Change sketch radius precision ("+str(prec)+")",
            "Change coplanar tolerance ("+str(coplanar_tol)+")",
            "Cancel"]
        item,ok = QtGui.QInputDialog.getItem(window,'Mesh Remodel v'+__version__,'Settings\n\nSelect the settings option\n',items,0,False,windowFlags)
        if ok and item == items[-1]:
            return
        elif ok and item == items[0]:
            keep = True
            pg.SetBool('KeepToolbar', keep)
        elif ok and item==items[1]:
            keep = False
            pg.SetBool('KeepToolbar', keep)
        elif ok and item==items[2]:
            new_point_size,ok = QtGui.QInputDialog.getDouble(window,"Point size", "Enter point size", point_size,1,50,2)
            if ok:
                pg.SetFloat("PointSize", new_point_size)
        elif ok and item==items[3]:
            new_line_width,ok = QtGui.QInputDialog.getDouble(window,"Line width", "Enter line width", line_width,1,50,2)
            if ok:
                pg.SetFloat("LineWidth", new_line_width)
        elif ok and item==items[4]:
            new_prec, ok = QtGui.QInputDialog.getInt(window,"Sketch Radius Precision", "\n\
-1 = no radius constraints\n\
0 = radius constraints\n\
1-12 = decimal points to round to when constraining\n\n\
Enter new sketch radius precision", prec, -1,12,1,flags=windowFlags)
            if ok:
                pg.SetInt("SketchRadiusPrecision", new_prec)
        elif ok and item==items[5]:
            new_coplanar_tol, ok = QtGui.QInputDialog.getDouble(window,"Coplanar tolerance", "Enter coplanar tolerance\n(Only applies to internal coplanar check using Alt+Click to create coplanar points)", coplanar_tol,.0000001,1,8)
            if ok:
                pg.SetFloat("CoplanarTolerance", new_coplanar_tol)
        return

    def IsActive(self):
        return True

#end settings class


####################################################################################
# Create the Mesh Remodel Points Object

class MeshRemodelCreatePointsObjectCommandClass(object):
    """Create Points Object command"""

    def __init__(self):
        self.mesh = None

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreatePointsObject.png') ,
            'MenuText': "Create points &object" ,
            'ToolTip' : "Create the points object"}
 
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        point_size = pg.GetFloat("PointSize",4.0)
        QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        doc.openTransaction("Create points object")

        pts=[]
        meshpts = self.mesh.Mesh.Points
        for m in meshpts:
            p = Part.Point(m.Vector)
            pts.append(p.toShape())
        Part.show(Part.makeCompound(pts),"MR_Points")
        doc.ActiveObject.ViewObject.PointSize = point_size
        doc.recompute()
        doc.commitTransaction()
        doc.recompute()
        QtGui.QApplication.restoreOverrideCursor()
        return

    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelectionEx()
        if len(sel) == 0:
            return False
        elif "Mesh.Feature" not in str(type(sel[0].Object)):
            return False
        else:
            self.mesh = sel[0].Object
        return True

# end create points class
####################################################################################
# Create the Mesh Remodel WireFrame Object

class MeshRemodelCreateWireFrameObjectCommandClass(object):
    """Create WireFrame Object command"""

    def __init__(self):
        self.mesh = None
        self.pb = None
        self.btn = None
        self.bar = None
        self.bCanceled = False

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateWireFrameObject.png') ,
            'MenuText': "Create Wire&Frame object" ,
            'ToolTip' : "Create the WireFrame object"}

    def makeId(self,a,b):
        """makeId(a,b)
           make a unique id from 2 integers using the Cantor Pairing Function
           see https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
        """
        if a >= b:
            return (a + b) * (a + b + 1) / 2 + a
        else:
            return (a + b) * (a + b + 1) / 2 + b

    def Activated(self):
        doc = FreeCAD.ActiveDocument
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        line_width = pg.GetFloat("LineWidth",5.0)
        point_size = pg.GetFloat("PointSize",4.0)
        tolerance = pg.GetFloat("CoplanarTolerance",.001)
        #QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        mids = []
        lines=[]
        meshfacets = self.mesh.Mesh.Facets
        total = len(meshfacets)
        ii = 0
        self.btn = QtGui.QPushButton("Cancel")
        self.btn.setToolTip("Cancel Create WireFrame Object command")
        self.btn.clicked.connect(self.on_clicked)
        self.pb = QtGui.QProgressBar()
        self.bar = Gui.getMainWindow().statusBar()
        self.bar.addWidget(self.pb)
        self.bar.addWidget(self.btn)
        self.btn.show()
        self.pb.show()
        self.pb.reset()
        self.pb.setMinimum(0)
        self.pb.setMaximum(total);
        self.pb.setFormat("%v/%m")
        idMap = {}
        for f in meshfacets:
            try:
                pts = f.Points
                (id0,id1,id2) = f.PointIndices
                if gu.isSamePoint(pts[0], pts[1], tolerance) or gu.isSamePoint(pts[0],pts[2],tolerance) or gu.isSamePoint(pts[1],pts[2],tolerance):
                    continue #otherwise we get an exception, so skip this line
                id = self.makeId(id0,id1)
                if not idMap.get(id):
                    l = Part.makeLine(pts[0],pts[1])
                    lines.append(l)
                    idMap[id] = id
                id = self.makeId(id0,id2)
                if not idMap.get(id):
                    l = Part.makeLine(pts[0],pts[2])
                    lines.append(l)
                    idMap[id] = id
                id = self.makeId(id2,id1)
                if not idMap.get(id):
                    l = Part.makeLine(pts[2],pts[1])
                    lines.append(l)
                    idMap[id] = id
            except:
                FreeCAD.Console.PrintError("Exception creating wireframe.  Typically this is caused by defective mesh object.\n")
                self.bCanceled = True
            ii += 1
            self.pb.setValue(ii)
            if ii % 10 == 0:
                QtGui.QApplication.processEvents()
            if self.bCanceled or Gui.getMainWindow().isVisible == False: #user exited with this still going, so quit
                self.bCanceled = False #for the next go around
                self.bar.removeWidget(self.pb)
                self.bar.removeWidget(self.btn)
                #QtGui.QApplication.restoreOverrideCursor()
                FreeCAD.Console.PrintMessage("MeshRemodel: WireFrame creation canceled\n")
                return
        self.bar.removeWidget(self.pb)
        self.bar.removeWidget(self.btn)
        doc.openTransaction("Create WireFrame object")
        Part.show(Part.makeCompound(lines),"MR_WireFrame")
        doc.ActiveObject.ViewObject.PointSize = point_size
        doc.ActiveObject.ViewObject.LineWidth = line_width
        doc.recompute()
        doc.commitTransaction()
        doc.recompute()
        #QtGui.QApplication.restoreOverrideCursor()
        return

    def on_clicked(self):
        self.bCanceled = True
        return

    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelectionEx()
        if len(sel) == 0:
            return False
        elif "Mesh.Feature" not in str(type(sel[0].Object)):
            return False
        else:
            self.mesh = sel[0].Object
        return True

# end create WireFrame class
####################################################################################
# Create the Mesh Cross Sections Object

class MeshRemodelCreateCrossSectionsCommandClass(object):
    """Use Mesh Design workbench Cross-Sections command"""

    def __init__(self):
        self.mesh = None

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateCrossSections.png') ,
            'MenuText': "Create cross-sections ob&ject..." ,
            'ToolTip' : "Create the cross-sections object\n\
Convenience link to the Mesh Design workbench cross-sections tool\n\
(These objects should not be directly used as wires, but rather as references for\n\
creating wires using the Mesh Remodel workbench as with the Points and WireFrame objects)\n\
"}
 
    def Activated(self):
        import MeshPartGui, FreeCADGui
        FreeCADGui.runCommand('MeshPart_CrossSections')
        return
   
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelectionEx()
        if len(sel) == 0:
            return False
        elif "Mesh.Feature" not in str(type(sel[0].Object)):
            return False
        else:
            self.mesh = sel[0].Object
        return True

# end open mesh section class
####################################################################################
# Create the Mesh Remodel Point Object

class MeshRemodelCreatePointObjectCommandClass(object):
    """Create Point Object command"""

    def __init__(self):
        self.obj = None
        self.pts = None
    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreatePointObject.png') ,
            'MenuText': "Create poin&t object" ,
            'ToolTip' : "Create a single point from the selected point"}
 
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        point_size = pg.GetFloat("PointSize",4.0)
        if len(global_picked) == 1:
            self.pts = global_picked #use preselect-picked points
        doc.openTransaction("Create point object")
        pt = doc.addObject("Part::Vertex", "MR_Point")
        if self.pts:
            pt.X = self.pts[0].x
            pt.Y = self.pts[0].y
            pt.Z = self.pts[0].z
        elif hasattr(self.obj,"Point"):
            pt.X = self.obj.Point.x
            pt.Y = self.obj.Point.y
            pt.Z = self.obj.Point.z
        elif hasattr(self.obj,"x"): #was a picked point along an edge rather than a vertex
            pt.X = self.obj.x
            pt.Y = self.obj.y
            pt.Z = self.obj.z
        
        doc.ActiveObject.ViewObject.PointSize = point_size
        doc.recompute()
        doc.commitTransaction()
        return
   
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        selobj = Gui.Selection.getSelectionEx()
        if selobj:
            sel = selobj[0].SubObjects
            if len(sel) != 1:
                return False;
            if "Vertex" in str(type(sel[0])):
                self.obj = sel[0]
                return True
            if "Edge" in str(type(sel[0])) or "Face" in str(type(sel[0])):
                if len (selobj[0].PickedPoints) == 1:
                    self.obj = selobj[0].PickedPoints[0]
                    return True
        return False

# end create point class

####################################################################################
# Create a coplanar points object from 3 selected points

class MeshRemodelCreateCoplanarPointsObjectCommandClass(object):
    """Create coplanar points object from 3 selected points"""

    def __init__(self):
        self.pts = []
        self.obj = None #original points object

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateCoplanar.png') ,
            'MenuText': "Create copla&nar points object" ,
            'ToolTip' : "\
Makes coplanar points object from 3 selected points, used to define the plane\n\
Uses internal coplanar check, (see settings -- Coplanar tolerance)\n\
(Alt+Click -- Creates empty sketch and places links to external geometry to coplanar points\n\
(Shift+Click for exploded compound, compatible with Shift+B block select)"}
 
    def Activated(self):
        if len(global_picked) == 3:
            self.pts = global_picked #use preselect-picked points
        if gu.isColinear(self.pts[0],self.pts[1],self.pts[2]):
            FreeCAD.Console.PrintError('Please select 3 non-colinear points in the plane\n')
            return
        doc = FreeCAD.ActiveDocument
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        line_width = pg.GetFloat("LineWidth",5.0)
        point_size = pg.GetFloat("PointSize",4.0)
        coplanar_tolerance = pg.GetFloat("CoplanarTolerance", .001)
        #QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        doc.openTransaction("Create coplanar")
        modifiers = QtGui.QApplication.keyboardModifiers()
        trio = [self.pts[0],self.pts[1],self.pts[2]]
        candidates = []
        if self.obj and hasattr(self.obj,"Shape"):
            candidates = self.obj.Shape.Vertexes
        coplanar = []
        bMakeSketch = False;
        if modifiers == QtCore.Qt.AltModifier or modifiers == QtCore.Qt.AltModifier.__or__(QtCore.Qt.ShiftModifier):
            bMakeSketch = True
        for v in candidates:
            if gu.isCoplanar(trio,v.Point,coplanar_tolerance):
                planar = True
            else:
                planar = False
            if planar:
                coplanar.append(Part.Point(v.Point).toShape())
        coplanar.extend([Part.Point(v).toShape() for v in trio])

        Part.show(Part.makeCompound(coplanar),"MR_Points_Coplanar")
        doc.ActiveObject.ViewObject.PointSize = point_size
        mr = doc.ActiveObject
        if bMakeSketch:
            if not "Sketcher_NewSketch" in Gui.listCommands():
                Gui.activateWorkbench("SketcherWorkbench")
                Gui.activateWorkbench("MeshRemodelWorkbench")
            Gui.runCommand("Sketcher_NewSketch")
            sketch=doc.ActiveObject
            sketch.Label = mr.Name+'_Sketch'
            sketch.MapReversed = False
            for ii in range(0,len(mr.Shape.Vertexes)):
                vname = 'Vertex'+str(ii+1)
                sketch.addExternal(mr.Name, vname)
        
            doc.recompute()

        if self.obj and hasattr(self.obj,"ViewObject"):
            self.obj.ViewObject.Visibility = False
        doc.recompute()
        doc.commitTransaction()
        if modifiers == QtCore.Qt.ShiftModifier or modifiers == QtCore.Qt.ShiftModifier.__or__(QtCore.Qt.AltModifier):
            doc.openTransaction("explode coplanar points")
            import CompoundTools.Explode
            input_obj = doc.ActiveObject
            comp = CompoundTools.Explode.explodeCompound(input_obj)
            input_obj.ViewObject.hide()
            for obj in comp[1]:
                obj.ViewObject.PointSize = point_size
            doc.recompute()
            doc.commitTransaction()
        #QtGui.QApplication.restoreOverrideCursor()
        return

    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelectionEx()
        if len(sel) == 0:
            return False
        if not hasattr(sel[0],"PickedPoints"):
            return False
        count = 0
        self.pts = []
        self.obj = sel[0].Object
        for s in sel:
            if hasattr(s,"PickedPoints"):
                p = s.PickedPoints
                for pt in s.PickedPoints:
                    self.pts.append(pt)
                    count += 1
                if len(p)==0: #might be individual part point objects
                    if len(s.Object.Shape.Vertexes)==1:
                        self.pts.append(s.Object.Shape.Vertexes[0].Point)
                        count += 1
            if count > 3:
                return False
        if count == 3:
            return True
        return False

# end create coplanar points object
####################################################################################
# Create a line from 2 selected points

class MeshRemodelCreateLineCommandClass(object):
    """Create Line from 2 selected points or from selected edge"""

    def __init__(self):
        self.pts = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateLine.png') ,
            'MenuText': "Create &line" ,
            'ToolTip' : "Create a line from 2 selected points\n(Ctrl+Click to add midpoint)\n(Ctrl+Shift+Click for only midpoint)"}
 
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        line_width = pg.GetFloat("LineWidth",5.0)
        point_size = pg.GetFloat("PointSize",4.0)
        modifiers = QtGui.QApplication.keyboardModifiers()
        #ctrl + click to include midpoint
        #ctrl + shift + click to include only the midpoint
        #QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        doc.openTransaction("Create line")
        if len(global_picked) == 2:
            self.pts = global_picked #use preselect-picked points
        line = Part.makeLine(self.pts[0],self.pts[1])
        lineName = "MR_Ref"
        if not modifiers == QtCore.Qt.ControlModifier.__or__(QtCore.Qt.ShiftModifier):
            #Part.show(line,"MR_Line")
            l = Draft.makeLine(line.Vertexes[0].Point, line.Vertexes[1].Point)
            l.Label="MR_"+l.Label
            lineName = doc.ActiveObject.Name
            doc.ActiveObject.ViewObject.LineWidth=line_width
            Gui.Selection.clearSelection()
            Gui.Selection.addSelection(doc.getObject(lineName))
            FreeCAD.Console.PrintMessage(lineName+": length = "+str(line.Length)+"\n  midpoint at "+str(gu.midpoint(line.firstVertex().Point,line.lastVertex().Point))+"\n")
        if modifiers == QtCore.Qt.ControlModifier or modifiers == QtCore.Qt.ControlModifier.__or__(QtCore.Qt.ShiftModifier):
            Part.show(Part.Point(gu.midpoint(line.firstVertex().Point,line.lastVertex().Point)).toShape(),lineName+"_Mid")
            doc.ActiveObject.ViewObject.PointSize = point_size
        doc.recompute()
        doc.commitTransaction()
        return
   
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelectionEx()
        if len(sel) == 0:
            return False
        if not hasattr(sel[0],"PickedPoints"):
            return False
        count = 0
        self.pts = []
        hasEdges = False
        for s in sel:
            if s.HasSubObjects and "Edge" in s.SubElementNames[0]:
                for sub in s.SubObjects:
                    if "Edge" in str(type(sub)):
                        self.pts.append(sub.firstVertex().Point)
                        self.pts.append(sub.lastVertex().Point)
                        count += 2
                        hasEdges = True
            if not hasEdges and hasattr(s,"PickedPoints"):
                p = s.PickedPoints
                for pt in s.PickedPoints:
                    self.pts.append(pt)
                    count += 1
                    if count > 2: #avoid parsing really long lists
                        return False
        if count == 2:
            return True
        return False

# end create line class

####################################################################################
# Create a Polygon from 3 or more selected points

class MeshRemodelCreatePolygonCommandClass(object):
    """Create Polygon from 3 or more selected points"""

    def __init__(self):
        self.pts = []
        self.edges = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreatePolygon.png') ,
            'MenuText': "Create &Polygon" ,
            'ToolTip' : "\
Create a Polygon from 3 or more selected points or 2 or more selected edges\n\
Might not always be coplanar, consider using links to external geometry in a sketch\n\
Do **not** attempt to mix selected edges and selected points, should be all edges\n\
or all points, but not a combination of the 2 object types\n\
(Makes individual lines, use Create wire to connect into a single wire object.)\n\
(Shift+Click to not close polygon) -- but selected edges never close unless connected\n\
(Alt+Click to sort selected points)\n\
"}
 
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        line_width = pg.GetFloat("LineWidth",5.0)
        point_size = pg.GetFloat("PointSize",4.0)
        #QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        if len(global_picked) > 2:
            self.pts = global_picked #use preselect-picked points
        doc.openTransaction("Create polygon")
        modifiers = QtGui.QApplication.keyboardModifiers()
        if modifiers != QtCore.Qt.ShiftModifier and modifiers != QtCore.Qt.ShiftModifier.__or__(QtCore.Qt.AltModifier):
            if len(self.pts) > 0:
                self.pts.append(self.pts[0]) #don't close polygon on shift+click

        if modifiers == QtCore.Qt.AltModifier.__or__(QtCore.Qt.ShiftModifier) or modifiers == QtCore.Qt.AltModifier:
            lineList = self.makePolygon(gu.sortPoints(self.pts))
        else:
            lineList = self.makePolygon(self.pts)

        lineObjs = []
        for line in lineList:
            #Part.show(line,"MR_Line")
            l = Draft.makeLine(line.Vertexes[0].Point, line.Vertexes[1].Point)
            l.Label="MR_"+l.Label
            doc.recompute()
            lineObjs.append(l)
        FreeCAD.Gui.Selection.clearSelection()
        for ll in lineObjs:
            ll.ViewObject.LineWidth = line_width
            FreeCAD.Gui.Selection.addSelection(ll)

        doc.commitTransaction()

        #QtGui.QApplication.restoreOverrideCursor()
        return

    def makePolygon(self,pts):
        """make list of lines out of the pts (vectors) list one line at a time, return the list
           or ignore pts if self.edges is not empty and make the list of lines out of those edges
        """
        lines=[]
        if len(self.edges) == 0:
            for ii in range(1,len(pts)):
                if pts[ii-1] != pts[ii]:
                    lines.append(Part.makeLine(pts[ii-1],pts[ii]))
        else:
            for edge in self.edges:
                lines.append(Part.makeLine(edge.firstVertex().Point, edge.lastVertex().Point))
        return lines

    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelectionEx()
        if len(sel) == 0:
            return False
        if not hasattr(sel[0],"PickedPoints"):
            return False
        count = 0
        self.pts = []
        self.edges = []
        for s in sel:
            if s.HasSubObjects and "Edge" in s.SubElementNames[0]:
                for sub in s.SubObjects:
                    if "Edge" in str(type(sub)):
                        self.edges.append(sub)
                count = len(self.edges)+1 #2 edges will work as well as 3 points
                continue
            if hasattr(s,"PickedPoints"):
                p = s.PickedPoints
                for pt in s.PickedPoints:
                    self.pts.append(pt)
                    count += 1
                if len(p)==0: #might be individual part point objects
                    if len(s.Object.Shape.Vertexes)==1:
                        self.pts.append(s.Object.Shape.Vertexes[0].Point)
                        count += 1
        if count >= 3:
            return True
        return False

# end create Polygon class
####################################################################################
# Create a BSpline from 3 or more selected points

class MeshRemodelCreateBSplineCommandClass(object):
    """Create BSpline from 3 or more selected points"""

    def __init__(self):
        self.pts = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateBSpline.png') ,
            'MenuText': "Create &BSpline" ,
            'ToolTip' : "Create a BSPline from 3 or more selected points\n(Shift+Click to not close bspline)\n(Alt+Click to sort selected points)"}
 
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        line_width = pg.GetFloat("LineWidth",5.0)
        point_size = pg.GetFloat("PointSize",4.0)
        if len(global_picked) > 2:
            self.pts = global_picked #use preselect-picked points
        #QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        doc.openTransaction("Create BSpline")
        modifiers = QtGui.QApplication.keyboardModifiers()
        is_periodic=True
        if modifiers == QtCore.Qt.ShiftModifier or modifiers == QtCore.Qt.ShiftModifier.__or__(QtCore.Qt.AltModifier):
            is_periodic=False #don't close bspline on shift+click
        if modifiers == QtCore.Qt.AltModifier or modifiers == QtCore.Qt.AltModifier.__or__(QtCore.Qt.ShiftModifier) or modifiers == QtCore.Qt.ShiftModifier:
            self.pts = gu.sortPoints(self.pts)[:-1]
        bs = Draft.makeBSpline(self.pts, is_periodic)
        bs.Label = "MR_BSpline"
        bs.ViewObject.LineWidth=line_width
        doc.recompute()
        doc.commitTransaction()
        #QtGui.QApplication.restoreOverrideCursor()
        return

    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelectionEx()
        if len(sel) == 0:
            return False
        if not hasattr(sel[0],"PickedPoints"):
            return False
        count = 0
        self.pts = []
        for s in sel:
            if hasattr(s,"PickedPoints"):
                p = s.PickedPoints
                for pt in s.PickedPoints:
                    self.pts.append(pt)
                    count += 1
                if len(p)==0: #might be individual part point objects
                    if len(s.Object.Shape.Vertexes)==1:
                        self.pts.append(s.Object.Shape.Vertexes[0].Point)
                        count += 1
        if count >= 3:
            return True
        return False

# end create BSpline class

###################################################################

# Create a Circle from first 3 selected points

class MeshRemodelCreateCircleCommandClass(object):
    """Create Circle from first 3 selected points"""

    def __init__(self):
        self.pts = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateCircle.png') ,
            'MenuText': "Create &circle" ,
            'ToolTip' : "Create a circle from first 3 selected points\n(Ctrl+Click to include Center point)\n(Ctrl+Shift+Click for only center)"}
 
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        line_width = pg.GetFloat("LineWidth",5.0)
        point_size = pg.GetFloat("PointSize",4.0)
        #QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)


        #add_varargs_method("makeCircle",&Module::makeCircle,
        #    "makeCircle(radius,[pnt,dir,angle1,angle2]) -- Make a circle with a given radius\n"
        #    "By default pnt=Vector(0,0,0), dir=Vector(0,0,1), angle1=0 and angle2=360"
        #);

        modifiers = QtGui.QApplication.keyboardModifiers()
        if len(global_picked) > 2:
            self.pts = global_picked
        poly = Part.makePolygon(self.pts)
        #Part.show(poly)
        normal = DraftGeomUtils.getNormal(poly)
        A = self.pts[0]
        B = self.pts[1]
        C = self.pts[2]

        if gu.isColinear(A,B,C):
            FreeCAD.Console.PrintError("MeshRemodel Error: Cannot make arc/circle from 3 colinear points\n")
            return

        center = gu.circumcenter(A,B,C)
        radius = gu.circumradius(A,B,C)

        doc.openTransaction("Create circle")
        circle = Part.makeCircle(radius, center, normal)
        circName="MR_Ref"
        if not modifiers == QtCore.Qt.ShiftModifier.__or__(QtCore.Qt.ControlModifier):
            #Part.show(circle,"MR_Circle")
            c = Draft.makeCircle(circle)
            c.Label="MR_"+c.Label
            circName = doc.ActiveObject.Name
            doc.ActiveObject.ViewObject.LineWidth=line_width
            FreeCAD.Console.PrintMessage(circName+": radius = "+str(radius)+"\n  center at "+str(center)+"\n")
            Gui.Selection.clearSelection()
            Gui.Selection.addSelection(doc.getObject(circName))
        if modifiers == QtCore.Qt.ControlModifier or modifiers==QtCore.Qt.ControlModifier.__or__(QtCore.Qt.ShiftModifier):
            Part.show(Part.Point(center).toShape(),circName+"_Ctr") #show the center point on ctrl click or shift+ctrl click
            doc.ActiveObject.ViewObject.PointSize = point_size

        doc.recompute()
        doc.commitTransaction()
        #QtGui.QApplication.restoreOverrideCursor()
        return



    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelectionEx()
        if len(sel) == 0:
            return False
        if not hasattr(sel[0],"PickedPoints"):
            return False
        count = 0
        self.pts = []
        for s in sel:
            if hasattr(s,"PickedPoints"):
                p = s.PickedPoints
                for pt in s.PickedPoints:
                    self.pts.append(pt)
                    count += 1
                if len(p)==0: #might be individual part point objects
                    if len(s.Object.Shape.Vertexes)==1:
                        self.pts.append(s.Object.Shape.Vertexes[0].Point)
                        count += 1
        if count >= 3:
            return True
        return False

# end create circle class

####################################################################################
# Create an Arc from first 3 selected points
class MeshRemodelCreateArcCommandClass(object):

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateArc.png') ,
            'MenuText': "Create &arc" ,
            'ToolTip' : "Create an arc from first 3 selected points\n(Ctrl+Click to include Center point)\n\
(Ctrl+Shift+Click for only center)\n\
(Alt+Click for all permutations possible with the 3 selected points -- 6 arcs -- useful where you\n\
are not getting the arc orientation you were expecting -- will need to delete unwanted extra arcs)"}

    def __init__(self):
        self.pts = []

    def Activated(self):
        doc = FreeCAD.ActiveDocument
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        line_width = pg.GetFloat("LineWidth",5.0)
        point_size = pg.GetFloat("PointSize",4.0)
        modifiers = QtGui.QApplication.keyboardModifiers()
        if len(global_picked) > 2:
            self.pts = global_picked
        poly = Part.makePolygon(self.pts)
#        Part.show(poly)
        normal = DraftGeomUtils.getNormal(poly)
        A = self.pts[0]
        B = self.pts[1]
        C = self.pts[2]

        if gu.isColinear(A,B,C):
            FreeCAD.Console.PrintError("MeshRemodel Error: Cannot make arc/circle from 3 colinear points\n")
            return
        center = gu.circumcenter(A,B,C)
        radius = gu.circumradius(A,B,C)

        doc.openTransaction("Create Arc")
        #arc = Part.ArcOfCircle(A,B,C)
        #on ctrl+shift click we only show center
        arcName="MR_Ref"
        if not modifiers == QtCore.Qt.ControlModifier.__or__(QtCore.Qt.ShiftModifier) or modifiers == QtCore.Qt.AltModifier:
            #Part.show(arc.toShape(),"MR_Arc")
            a = Draft.make_arc_3points([A,B,C])
            a.Label = "MR_"+a.Label
            if modifiers == QtCore.Qt.AltModifier:
                a = Draft.make_arc_3points([A,C,B])
                a.Label = "MR_"+a.Label
                a = Draft.make_arc_3points([B,C,A])
                a.Label = "MR_"+a.Label
                a = Draft.make_arc_3points([B,A,C])
                a.Label = "MR_"+a.Label
                a = Draft.make_arc_3points([C,B,A])
                a.Label = "MR_"+a.Label
                a = Draft.make_arc_3points([C,A,B])
                a.Label = "MR_"+a.Label
            arcName = a.Label
            doc.ActiveObject.ViewObject.LineWidth=line_width
            FreeCAD.Console.PrintMessage(arcName+": radius = "+str(radius)+"\n  center at "+str(center)+"\n")
            Gui.Selection.clearSelection()
            Gui.Selection.addSelection(a)
        if modifiers == QtCore.Qt.ControlModifier or modifiers == QtCore.Qt.ControlModifier.__or__(QtCore.Qt.ShiftModifier):
            Part.show(Part.Point(center).toShape(),arcName+"_Ctr") #show the center point
            doc.ActiveObject.ViewObject.PointSize = point_size
        doc.recompute()
        doc.commitTransaction()
        #QtGui.QApplication.restoreOverrideCursor()
        return

    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelectionEx()
        if len(sel) == 0:
            return False
        if not hasattr(sel[0],"PickedPoints"):
            return False
        count = 0
        self.pts = []
        for s in sel:
            if hasattr(s,"PickedPoints"):
                p = s.PickedPoints
                for pt in s.PickedPoints:
                    self.pts.append(pt)
                    count += 1
                if len(p)==0: #might be individual part point objects
                    if len(s.Object.Shape.Vertexes)==1:
                        self.pts.append(s.Object.Shape.Vertexes[0].Point)
                        count += 1
        if count >= 3:
            return True
        return False

# end create arc class

####################################################################################

# Make a sketch from selected objects
class MeshRemodelCreateSketchCommandClass(object):
    """Create sketch from selected objects"""

    def __init__(self):
        self.objs = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateSketch.png') ,
            'MenuText': "Create s&ketch" ,
            'ToolTip' : "\
Create a new empty sketch, optionally attaching to selected objects, e.g. 3 points to define a plane.\n\
(Ctrl+Click to make a sketch out of selected objects, e.g. circles, polygons, etc.)\n\
(Alt+Click to make a separate sketch from each selected object, and then merge them together.)\n\
(Shift+Click to make a new points object from all picked points, attach a new sketch to first 3 points, add all\n\
picked points as links to external geometry -- useful to produce a coplanar wire)\n\
"}
 
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        modifiers = QtGui.QApplication.keyboardModifiers()
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        prec = pg.GetInt("SketchRadiusPrecision", 1)

        if modifiers == QtCore.Qt.NoModifier:
            if not "Sketcher_NewSketch" in Gui.listCommands():
                Gui.activateWorkbench("SketcherWorkbench")
                Gui.activateWorkbench("MeshRemodelWorkbench")
            Gui.runCommand("Sketcher_NewSketch")
            return
        if modifiers == QtCore.Qt.AltModifier:
            #alternative method: on alt+click make separate sketch from each object, then merge them together
            sketches=[]
            for obj in self.objs:
                sketches.append(Draft.makeSketch(obj,autoconstraints=True,radiusPrecision=prec))
            doc.recompute()
            FreeCADGui.Selection.clearSelection()
            for sk in sketches:
                if sk:
                    FreeCADGui.Selection.addSelection(sk)
            if len(sketches) >= 2:
                if not "Sketcher_NewSketch" in Gui.listCommands():
                    Gui.activateWorkbench("SketcherWorkbench")
                    Gui.activateWorkbench("MeshRemodelWorkbench")
                FreeCADGui.runCommand("Sketcher_MergeSketches")

            sketch = doc.ActiveObject
            doc.recompute()
            for sk in sketches:
                if sk:
                    doc.removeObject(sk.Name)
        elif modifiers == QtCore.Qt.ControlModifier:
            #on ctrl+click make single sketch out of selected objects
            sketch = Draft.makeSketch(self.objs,autoconstraints=True,radiusPrecision=prec)
            doc.recompute()
        elif modifiers == QtCore.Qt.ShiftModifier:
            #on shift+click map sketch to first 3 picked points as a plane, add all picked points as links to external geometry
            sel = FreeCADGui.Selection.getSelectionEx()
            picked = []
            if global_picked:
                picked = global_picked
            else:
                for s in sel:
                    picked.extend(list(s.PickedPoints))
            if len(picked) < 3:
                FreeCAD.Console.PrintMessage("MeshRemodel: Not enough picked points, need at least 3.\n")
                return
            #make a compound of the picked points
            part_pts = []
            for m in picked:
                p = Part.Point(m)
                part_pts.append(p.toShape())
            Part.show(Part.makeCompound(part_pts),"MR_Picked_Points")
            sk_pts = doc.ActiveObject
            doc.recompute()
            FreeCADGui.Selection.clearSelection()
            FreeCADGui.Selection.addSelection(sk_pts, ["Vertex1", "Vertex2", "Vertex3"])
            if not "Sketcher_NewSketch" in Gui.listCommands():
                Gui.activateWorkbench("SketcherWorkbench")
                Gui.activateWorkbench("MeshRemodelWorkbench")
            Gui.runCommand("Sketcher_NewSketch")
            sketch=doc.ActiveObject
            sketch.Label = 'MR_Picked_Sketch'
            sketch.MapReversed = False
            for ii in range(0,len(sk_pts.Shape.Vertexes)):
                vname = 'Vertex'+str(ii+1)
                sketch.addExternal(sk_pts.Name, vname)
            doc.recompute()

        for o in self.objs:
            if hasattr(o,"ViewObject"):
                o.ViewObject.Visibility=False

        doc.recompute()
        return

    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelectionEx()
        if len(sel) == 0:
            return False
        count = 0
        self.objs = []
        for s in sel:
            if hasattr(s,"Object"):
                self.objs.append(s.Object)
                count += 1
        if count >= 1:
            return True
        return False

# end create sketch class
####################################################################################

# Make a wire from selected objects
class MeshRemodelCreateWireCommandClass(object):
    """Create wire from selected objects"""

    def __init__(self):
        self.objs = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateWire.png') ,
            'MenuText': "Create &wire" ,
            'ToolTip' : "Create a wire from selected objects\n(All selected objects should be connected.)\n\
(Runs draft upgrade)\n\
Ctrl+Click to downgrade to edges\n\
Tip: You can also use this to upgrade a wire to a face, which can be converted to a sketch to avoid some coplanar issues\n"}
 
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        modifiers = QtGui.QApplication.keyboardModifiers()
        if (modifiers == QtCore.Qt.ControlModifier):
            doc.openTransaction("Downgrade to edges")
            Draft.downgrade(self.objs)
            doc.recompute()
        else:
            selbackup = FreeCAD.Gui.Selection.getSelection()
            doc.openTransaction("Create wire (upgrade)")
            Draft.upgrade(self.objs)
            doc.recompute()
            for o in self.objs:
                if hasattr(o,"ViewObject"):
                    o.ViewObject.Visibility=False
            FreeCAD.Gui.Selection.clearSelection()
            for sel in selbackup:
                FreeCAD.Gui.Selection.addSelection(sel)
        doc.commitTransaction()

        #QtGui.QApplication.restoreOverrideCursor()
        return
   
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelectionEx()
        if len(sel) == 0:
            return False
        count = 0
        self.objs = []
        for s in sel:
            if hasattr(s,"Object"):
                self.objs.append(s.Object)
                count += 1
        if count >= 1:
            return True
        return False

# end create wire class

####################################################################################
class MeshRemodelSelectionObserver():
    """Selection Observer to select preselected vertices"""

    def __init__(self):
        self.mode = ["Vertex"] #can also be ["Edge","Vertex"] or ["Edge"]
        pass

    def setMode(self,mode):
        self.mode = mode

    def getMode(self):
        return self.mode

    def isEdgeMode(self):
        return "Edge" in self.mode

    def setEdgeMode(self):
        if not "Edge" in self.mode:
            self.mode.append(["Edge"])

    def setVertexMode(self):
        if not "Vertex" in self.mode:
            self.mode.append(["Vertex"])

    def isVertexMode(self):
        return "Vertex" in self.mode

    def setPreselection(self,doc,obj,sub):                # Preselection object
        modifiers = QtGui.QApplication.keyboardModifiers()
        if not modifiers == QtCore.Qt.ControlModifier:
            return
        if self.isVertexMode() and "Vertex" in str(sub):
            Gui.Selection.addSelection(doc,obj,str(sub))
            idx = int(sub[6:])
            thisobj = FreeCAD.ActiveDocument.getObject(obj)
            p = thisobj.Shape.Vertexes[idx-1].Point
            if not gu.hasPoint(p,global_picked,.0001):
                global_picked.append(p)

        elif self.isEdgeMode() and "Edge" in str(sub):
            Gui.Selection.addSelection(doc,obj,str(sub))

    def addSelection(self,doc,obj,sub,pnt):               # Selection object
        pass
        #App.Console.PrintMessage("addSelection"+ "\n")
        #App.Console.PrintMessage(str(doc)+ "\n")          # Name of the document
        #App.Console.PrintMessage(str(obj)+ "\n")          # Name of the object
        #App.Console.PrintMessage(str(sub)+ "\n")          # The part of the object name
        #App.Console.PrintMessage(str(pnt)+ "\n")          # Coordinates of the object
        #App.Console.PrintMessage("______"+ "\n")

    def removeSelection(self,doc,obj,sub):                # Delete the selected object
        #FreeCAD.Console.PrintMessage("removeSelection"+ "\n")
        if self.isVertexMode() and "Vertex" in str(sub):
            idx = int(sub[6:])
            thisobj = FreeCAD.ActiveDocument.getObject(obj)
            p = thisobj.Shape.Vertexes[idx-1].Point
            if gu.hasPoint(p, global_picked, .0001):
                for ii in range(0,len(global_picked)):
                    try:
                        if gu.isSamePoint(global_picked[ii],p,.0001):
                            global_picked.remove(global_picked[ii])
                    except:
                        pass
        pass

    def setSelection(self,doc):                           # Selection in ComboView
        #App.Console.PrintMessage("setSelection"+ "\n")
        pass

    def clearSelection(self,doc):                         # If click on the screen, clear the selection
        #FreeCAD.Console.PrintMessage("clearSelection"+ "\n")  # If click on another object, clear the previous object
        global_picked.clear()
        pass


#end MeshRemodelSelectionObserver class

#####################################################################

# add a selection observer
class MeshRemodelAddSelectionObserverCommandClass(object):
    """Add a selection observer to enable automatic selection of edges or vertices on preselection"""

    def __init__(self):
        self.sel = None
        self.btn = None
        self.bar = None

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'AddSelectionObserver.png') ,
            'MenuText': "Add Selection &Observer" ,
            'ToolTip' : "Allows to select vertices by Ctrl + preselection\n\
This is intended for use with BSpline selection of points, but\n\
can be used with all MeshRemodel tools that use selected points.  DOES NOT work if you \n\
mix points selected in this mode with normal selection mode in the same operation\n"}

    def on_clicked(self):
        self.bar.removeWidget(self.btn)
        FreeCADGui.Selection.removeObserver(self.sel)
        self.bar = None
        self.btn = None
        self.sel = None
        global_picked.clear()

    def Activated(self):
        doc = FreeCAD.ActiveDocument
        if not self.sel:
            self.sel = MeshRemodelSelectionObserver()
            FreeCADGui.Selection.addObserver(self.sel)
            self.btn = QtGui.QPushButton("Cancel Preselect Mode")
            self.btn.setToolTip("Cancel Preselect selection mode and remove observer")
            self.btn.clicked.connect(self.on_clicked)
            self.bar = Gui.getMainWindow().statusBar()
            self.bar.addWidget(self.btn)
            self.btn.show()
            global_picked.clear()
        else:
            FreeCADGui.Selection.removeObserver(self.sel)
            self.sel = None
            self.bar.removeWidget(self.btn)
            FreeCADGui.Selection.removeObserver(self.sel)
            self.bar = None
            self.btn = None
            global_picked.clear()
        return
   
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        return True

# end MeshRemodelAddSelectionObserverCommandClass(object)
####################################################################################


# Merge selected sketches
class MeshRemodelMergeSketchesCommandClass(object):
    """Merge selected sketches"""

    def __init__(self):
        self.objs = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'MergeSketches.png') ,
            'MenuText': "&Merge sketches" ,
            'ToolTip' : "Merge selected sketches with sketcher merge sketches tool"}
 
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        #doc.openTransaction("Merge sketches")  #not needed since the command does this
        if not "Sketcher_NewSketch" in Gui.listCommands():
            Gui.activateWorkbench("SketcherWorkbench")
            Gui.activateWorkbench("MeshRemodelWorkbench")
        Gui.runCommand("Sketcher_MergeSketches")
        doc.recompute()
        for o in self.objs:
            if hasattr(o,"ViewObject"):
                o.ViewObject.Visibility=False
        #doc.commitTransaction()
        #QtGui.QApplication.restoreOverrideCursor()
        Gui.activateWorkbench("MeshRemodelWorkbench")
        return
   
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelectionEx()
        if len(sel) == 0:
            return False
        count = 0
        self.objs = []
        for s in sel:
            if hasattr(s,"Object"):
                if "Sketch" in s.Object.Name:
                    self.objs.append(s.Object)
                    count += 1
        if count >= 2:
            return True
        return False

# end merge sketches

####################################################################################

# validate sketch
class MeshRemodelValidateSketchCommandClass(object):
    """Validate sketch"""

    def __init__(self):
        self.objs = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'ValidateSketch.png') ,
            'MenuText': "&Validate sketch" ,
            'ToolTip' : "Validate selected sketch with sketcher validate sketch tool"}
 
    def Activated(self):
        if not "Sketcher_NewSketch" in Gui.listCommands():
            Gui.activateWorkbench("SketcherWorkbench")
            Gui.activateWorkbench("MeshRemodelWorkbench")
        Gui.runCommand("Sketcher_ValidateSketch")
        return
   
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelectionEx()
        if len(sel) == 0:
            return False
        count = 0
        self.objs = []
        for s in sel:
            if hasattr(s,"Object"):
                if "Sketch" in s.Object.Name:
                    self.objs.append(s.Object)
                    count += 1
        if count == 1:
            return True
        return False

# end validate sketch

##################################################################################################
# initialize

def initialize():
    if FreeCAD.GuiUp:
        Gui.addCommand("MeshRemodelCreatePointsObject", MeshRemodelCreatePointsObjectCommandClass())
        Gui.addCommand("MeshRemodelCreateWireFrameObject",MeshRemodelCreateWireFrameObjectCommandClass())
        Gui.addCommand("MeshRemodelCreateCrossSectionsObject",MeshRemodelCreateCrossSectionsCommandClass())
        Gui.addCommand("MeshRemodelAddSelectionObserver",MeshRemodelAddSelectionObserverCommandClass())
        Gui.addCommand("MeshRemodelCreatePointObject", MeshRemodelCreatePointObjectCommandClass())
        Gui.addCommand("MeshRemodelCreateCoplanarPointsObject", MeshRemodelCreateCoplanarPointsObjectCommandClass())
        Gui.addCommand("MeshRemodelCreateLine", MeshRemodelCreateLineCommandClass())
        Gui.addCommand("MeshRemodelCreatePolygon", MeshRemodelCreatePolygonCommandClass())
        Gui.addCommand("MeshRemodelCreateBSpline", MeshRemodelCreateBSplineCommandClass())
        Gui.addCommand("MeshRemodelCreateCircle", MeshRemodelCreateCircleCommandClass())
        Gui.addCommand("MeshRemodelCreateArc", MeshRemodelCreateArcCommandClass())
        Gui.addCommand("MeshRemodelCreateWire", MeshRemodelCreateWireCommandClass())
        Gui.addCommand("MeshRemodelCreateSketch", MeshRemodelCreateSketchCommandClass())
        Gui.addCommand("MeshRemodelMergeSketches", MeshRemodelMergeSketchesCommandClass())
        Gui.addCommand("MeshRemodelValidateSketch", MeshRemodelValidateSketchCommandClass())
        Gui.addCommand("MeshRemodelSettings", MeshRemodelSettingsCommandClass())


initialize()
