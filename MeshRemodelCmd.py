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
__date__    = "2024.08.14"
__version__ = "1.10.1"

import FreeCAD, FreeCADGui, Part, os, math
from PySide import QtCore, QtGui
import Draft, DraftGeomUtils, DraftVecUtils, Mesh, MeshPart
import time
import numpy as np
import shiboken2 as shiboken


if FreeCAD.GuiUp:
    from FreeCAD import Gui

__dir__ = os.path.dirname(__file__)
iconPath = os.path.join( __dir__, 'Resources', 'icons' )
keepToolbar = False
windowFlags = QtCore.Qt.WindowTitleHint | QtCore.Qt.WindowCloseButtonHint #no ? in title bar
global_picked = [] #picked points list for use with selection by preselection observer
FC_VERSION = float(FreeCAD.Version()[0]) + float(FreeCAD.Version()[1]) #e.g. 0.20, 0.18, 1.??
epsilon = Part.Precision.confusion() #1e-07

def fixTip(tip):
    if FC_VERSION >= 0.20:
        return tip.replace("\n","<br/>")
    else:
        return tip




######################################################################################
# geometry utilities

class MeshRemodelGeomUtils(object):
    """Geometry Utilities"""

    #progress bar on status bar with cancel button
    class MRProgress:
        def __init__(self):
            self.pb = None
            self.btn = None
            self.bar = None
            self.bCanceled = False
            self.value = 0
            self.total = 0
            self.mw = FreeCADGui.getMainWindow()
            self.lastUpdate = time.time()

        def makeProgressBar(self,total=0,buttonText = "Cancel",tooltip = "Cancel current operation",updateInterval = .5):
            """total is max value for progress bar, mod = number of updates you want"""
            self.btn = QtGui.QPushButton(buttonText)
            self.btn.setToolTip(tooltip)
            self.btn.clicked.connect(self.on_clicked)
            self.pb = QtGui.QProgressBar()
            self.bar = self.mw.statusBar()
            self.bar.addWidget(self.pb)
            self.bar.addWidget(self.btn)
            self.btn.show()
            self.pb.show()
            self.pb.reset()
            self.value = 0
            self.pb.setMinimum(0)
            self.updateInterval = updateInterval
            self.pb.setMaximum(total);
            self.total = total
            self.pb.setFormat("%v/%m")
            self.bAlive = True #hasn't been killed yet
            self.bCanceled = False

        def on_clicked(self):
            self.bCanceled = True
            self.killProgressBar()

        def isCanceled(self):
            self.value += 1
            timeNow = time.time()
            if timeNow - self.lastUpdate >= self.updateInterval:
                self.lastUpdate = timeNow
                self.pb.setValue(self.value)
                FreeCADGui.updateGui()
            if self.mw.isHidden():
                self.bCanceled = True
                self.killProgressBar()
            return self.bCanceled

        def killProgressBar(self):
            if self.bAlive: #check if it has already been removed before removing
                self.bar.removeWidget(self.pb)
                self.bar.removeWidget(self.btn)
                self.bAlive = False
            self.pb.hide()
            self.btn.hide()
            self.value = 0
            self.total = 0
            
    #### end progress bar class
    

    def getFloatFromUser(self, title, msg, value, min=-float("inf"), max=float("inf"), flags=None, step=0.1):
        """Open a dialog and get a floating point value from the user"""
        val,ok = QtGui.QInputDialog.getDouble(FreeCADGui.getMainWindow(), title, 
                msg, value, min, max, 8,FreeCADGui.getMainWindow().windowFlags(), step)
        return val if ok else None

    def getFacetsFromFacetIndices(self, facet_indices, mesh):
        """getFacetsFromFacetIndices(self, facet_indices, mesh)"""
        return [mesh.Facets[idx] for idx in facet_indices]

    def getFacetIndicesFromPointIndices(self, point_indices, mesh):
        """getFacetIndicesFromPointINdices(self, point_indices, mesh)
        return the facet indices given the point indices and the mesh"""
        return [facet.Index for facet in mesh.Facets \
                if self.facetIsInList(facet, point_indices)]
    
    def getPointIndicesInMesh(self, pts, mesh):
        """gets the point indices of the points in the mesh"""
        return [self.findPointInMesh(mesh, pt) for pt in pts]

    def facetIsInList(self, facet, point_indices):
        """all 3 of the facet's points must be in the point_indices"""
        for idx in facet.PointIndices:
            if not idx in point_indices:
                return False
        return True

    def findPointInMesh(self, mesh, pt):
        """find the point in the mesh, return the index of the point in
        the points list"""
        topology = mesh.Topology
        # topology is a tuple ([vectorlist],[facetlist])
        # facet list is a list of tuples
        # each tuple of the form (pt1_idx, pt2_idx, pt3_idx)
        indices = []
        points = mesh.Points
        for point in points:
            if gu.isSamePoint(pt, point.Vector, 1e-5):
                if not point.Index in indices:
                    indices.append(point.Index)
        if len(indices) == 1:
            return indices[0]
        return []
        
    def checkComponents(self,mesh):
        """check for multiple components in a mesh and warn user"""
        if mesh.countComponents() > 1:
            FreeCAD.Console.PrintWarning("\
MeshRemodel: selected mesh has multiple components.  Consider using \
Split by components operation in Mesh workbench to separate these into \
component objects before attempting modifications.\n")

    def getBaseAndNormal(self, trio):
        """return base,normal of plane defined by trio, list of Vector"""
        trio = np.array([np.array(v.Point) for v in trio[:3]])
        normal = np.cross(trio[1] - trio[0], trio[2] - trio[0])
        if not any(normal): #all 0's
            return (None,None)
        divisor = np.linalg.norm(normal)
        normal /= np.linalg.norm(normal) #normalize
        norm = FreeCAD.Vector(normal[0], normal[1], normal[2])
        return trio[0], norm

    def isCoplanar(self, trio, pt, tol=1e-5):
        """ check if pt is one the same plane as the one defined by
        trio, a trio of points.  pt is considered to be on the plane if
        its distance to the plane <= tol"""
        def distance_to_plane(point, plane_point, plane_normal):
            point = np.array(point)
            plane_point = np.array(plane_point)
            plane_normal = np.array(plane_normal)
            vector_to_point = point - plane_point
            distance = np.abs(np.dot(vector_to_point, plane_normal)) / np.linalg.norm(plane_normal)
            return distance
 
        trio = np.array([np.array(v.Point) for v in trio])
        pt = np.array(pt)
        normal = np.cross(trio[1] - trio[0], trio[2] - trio[0])
        normal /= np.linalg.norm(normal) #normalize
        return distance_to_plane(pt, trio[0], normal) <= tol

    def hasPoint(self,pt,lis,tol):
        """hasPoint(pt,lis,tol)"""
        for l in lis:
            if self.isSamePoint(pt,l,tol):
                return True
        return False
        
    def hasLine(self, line, lis, tol=1e-6):
        """hasLine(self, line, lis, tol=1e-6)"""
        for l in lis:
            if self.isSameLine(line, l, tol):
                return True
        return False

    def isSamePoint(self,A,B,tol):
        """isSamePoint(A,B,tol)"""
        dis = self.dist(A,B)
        if dis < tol:
            return True
        return False
        
    def isSameLine(self, A, B, tol=1e-6):
        """isSameLine(self, A, B, tol)"""
        if abs(A.Length - B.Length) > tol:
            return False
        if self.isSamePoint(A.Vertex1.Point, B.Vertex1.Point, tol):
            if self.isSamePoint(A.Vertex2.Point, B.Vertex2.Point, tol):
                return True
        if self.isSamePoint(A.Vertex1.Point, B.Vertex2.Point, tol):
            if self.isSamePoint(A.Vertex2.Point, B.Vertex1.Point, tol):
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
        
    def planeFromPoints(self, pts):
        """returns (base,normal) of pts by using pts[0], pts[mid], and pts[-1]
           to find the plane defined by those 3 points"""
        if len(pts) < 3:
            raise Exception("planeFromPoints requires 3 or more points")
        points = np.array([[v.x, v.y, v.z] for v in pts])
        
        # Define the plane using the three points: first, mid, last
        mid = int(len(pts)/2)
        #print(f"points used: {points[0], points[mid], points[-1]}")
        plane_origin = points[0]
        plane_normal = np.cross(np.array(points[mid]) - np.array(points[0]), np.array(points[-1]) - np.array(points[0]))
        plane_normal /= np.linalg.norm(plane_normal)
        base = FreeCAD.Vector(plane_origin[0], plane_origin[1], plane_origin[2])
        normal = FreeCAD.Vector(plane_normal[0], plane_normal[1], plane_normal[2])
        return (base,normal)

    def projectVectorToPlane(self, point, base, normal):
        """uses FreeCAD.Vectors, which are converted to numpy and returned as FreeCAD.Vector"""
        pt = np.array(point)
        b = np.array(base)
        n = np.array(normal)
        #print(f"point:{point},base:{base},normal:{normal}")
        projected = self.projectPointToPlane(pt, b, n)
        #print(f"projected: {projected}")
        return FreeCAD.Vector(projected[0], projected[1], projected[2])

    def projectPointToPlane(self, point, plane_origin, plane_normal):
        """project np point to np plane, return np array"""
        # Calculate vector from the plane origin to the point
        v = np.array(point) - np.array(plane_origin)
        
        # Project the vector onto the plane
        projected_v = v - np.dot(v, plane_normal) * plane_normal
        
        # Calculate the new position of the point on the plane
        projected_point = np.array(plane_origin) + projected_v
        
        return projected_point.tolist()    

    def projectPointsToPlane(self, pts):
        """find the plane defined by the first 3 points in the points list, and then
           project all remaining points to that plane using numpy"""

        points = np.array([[v.x, v.y, v.z] for v in pts])
        
        # Define the plane using the three points: first, mid, last
        mid = int(len(pts)/2)
        plane_origin = points[0]
        plane_normal = np.cross(np.array(points[mid]) - np.array(points[0]), np.array(points[-1]) - np.array(points[0]))
        plane_normal /= np.linalg.norm(plane_normal)
        
        # Project all points onto the plane
        projected_points = [self.projectPointToPlane(point, plane_origin, plane_normal) for point in points]
        FCPoints = [FreeCAD.Vector(x,y,z) for x,y,z in projected_points]
        return FCPoints
        
    def flattenPoints(self, pts, align_plane):
        """ project points to align_plane.
            pts is list of Part.Vertex objects
            align_plane is part::plane object or 
            can be any object with a face (face1 will be used)
            returns: new list of Part.Vertex objects on the plane"""
        plane = align_plane.Shape.Faces[0]
        normal = plane.normalAt(0,0)
        base = align_plane.Shape.Vertexes[0].Point
        verts = []
        flatPts = [] #eliminate duplicate points
        pb = gu.MRProgress()
        pb.makeProgressBar(len(pts),buttonText = "Cancel Phase1/2",tooltip="Cancel projecting points to plane")
        for p in pts:
            fpt = p.Point.projectToPlane(base,normal)
            if not self.hasPoint(fpt,flatPts,epsilon):
                flatPts.append(fpt)
            if pb.isCanceled():
                FreeCAD.Console.PrintWarning("Phase 1 / 2 canceled, object may be incomplete.\n")
                break
        pb.killProgressBar()

        pb.makeProgressBar(len(flatPts),"Cancel Phase2/2")
        for p in flatPts:
            verts.append(Part.Vertex(p))
            if pb.isCanceled():
                FreeCAD.Console.PrintWarning("Phase 2 / 2 canceled, object may be incomplete.\n")
                break
        pb.killProgressBar()
        return verts

    def nearestPoint(self, pt, pts, exclude):
        """ nearestPoint(pt, pts, exclude)
            pt is a vector, pts a list of vectors
            exclude is a list of vectors to exclude from process
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
        return {'Pixmap'  : os.path.join( iconPath , 'Settings.svg') , # the name of an icon file available in the resources
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
        coplanar_tol = pg.GetFloat("CoplanarTolerance",.01)
        wireframe_tol = pg.GetFloat("WireFrameTolerance",.01)
        items=[("","*")[keep]+"Keep the toolbar active",
            ("","*")[not keep]+"Do not keep the toolbar active",
            "Change point size ("+str(point_size)+")",
            "Change line width ("+str(line_width)+")",
            "Change sketch radius precision ("+str(prec)+")",
            "Change coplanar tolerance ("+str(coplanar_tol)+")",
            "Change wireframe tolerance("+str(wireframe_tol)+")",
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
            new_coplanar_tol, ok = QtGui.QInputDialog.getDouble(window,"Coplanar tolerance", "Enter coplanar tolerance\n(Used when creating coplanar points.  Increase if some points are missing.)", coplanar_tol,.0000001,1,8)
            if ok:
                pg.SetFloat("CoplanarTolerance", new_coplanar_tol)
        elif ok and item==items[6]:
            new_wireframe_tol, ok = QtGui.QInputDialog.getDouble(window,"Wireframe tolerance", "Enter wireframe tolerance\n(Used when creating wireframes to check if 2 points are the same.)", wireframe_tol,.0000001,1,8)
            if ok:
                pg.SetFloat("WireFrameTolerance", new_wireframe_tol)
        return

    def IsActive(self):
        return True

#end settings class


####################################################################################
# Create the Mesh Remodel Points Object

class PointsObject:
    """Part::FeaturePython object to serve as proxy for mesh object so we can have
    selectable points to work with"""
    def __init__(self, obj, meshObj=None, className="PointsObject"):
        obj.Proxy = self
        obj.addProperty("App::PropertyLink","Source",className,"Source mesh document object").Source = meshObj
        obj.addProperty("App::PropertyString", "OriginalDisplayMode",className,\
                "Used to restore mesh to original display mode upon deleting points object").OriginalDisplayMode=\
                meshObj.ViewObject.DisplayMode if meshObj else "Shaded"
        obj.addProperty("App::PropertyBool", "OriginalSelectable", className, \
                "Used to restore mesh to original Selectable status after deleting points object")\
                .OriginalSelectable = meshObj.ViewObject.Selectable if meshObj else True
        obj.setEditorMode("OriginalDisplayMode", 2) #hidden
        obj.setEditorMode("OriginalSelectable", 2)
        obj.addProperty("App::PropertyBool","FlatLines",className,\
                        "Whether to show mesh object in Flat Lines mode").FlatLines = True
        obj.addProperty("App::PropertyBool","OneSideLighting",className,\
                        "Whether to set mesh object Lighting to 'One side'").OneSideLighting = True
        obj.addProperty("App::PropertyFloat","PointSize", className,\
                    "Size of the points in pixels, changing also changes value in settings"\
                    ).PointSize=self.PointSize
        obj.addProperty("App::PropertyInteger","Transparency",className,\
                "Transparency of Source object").Transparency = 25
        obj.addProperty("App::PropertyBool", "Selectable", className, \
                "Toggle to false to make the mesh object non-selectable").Selectable = False

    @property
    def PointSize(self):
        """manage parameter"""
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        point_size = pg.GetFloat("PointSize",4.0)
        return point_size
    
    @PointSize.setter
    def PointSize(self, psize):
        """manager parameter"""
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        pg.SetFloat("PointSize", psize)
        
    def onChanged(self, fp, prop):
        if prop == "OneSideLighting" and fp.Source != None:
            if fp.OneSideLighting:
                fp.Source.ViewObject.Lighting = "One side"
            else:
                fp.Source.ViewObject.Lighting = "Two side"
        elif prop == "FlatLines" and fp.Source != None:
            if fp.FlatLines:
                fp.Source.ViewObject.DisplayMode = "Flat Lines"
            else:
                fp.Source.ViewObject.DisplayMode = "Shaded"
        elif prop == "PointSize":
            fp.ViewObject.PointSize = fp.PointSize
            self.PointSize = fp.PointSize
        elif prop == "Selectable":
            if hasattr(fp,"Source"):
                fp.Source.ViewObject.Selectable = fp.Selectable
            
    def execute(self, fp):
        if not fp.Source:
            fp.Shape = Part.Shape()
            return
        doc = fp.Document
        pts=[]
        if hasattr(fp.Source,"Mesh"):
            meshpts =fp.Source.Mesh.Points
            for m in meshpts:
                p = Part.Point(m.Vector)
                pts.append(p.toShape())
        #print(f"pts = {pts}")        
        fp.Shape = Part.makeCompound(pts)
        fp.ViewObject.PointSize = fp.PointSize
        fp.Source.ViewObject.Transparency = fp.Transparency
        
    def refresh(self, fp):
        """when the mesh changes we don't always get the usual touched notification,
        such as when a point is moved or a facet is added.  This is good because if
        a point is accidentally removed when removing facets the point structure will
        still be there in the points object proxy, at least until the user manually
        refreshes via this function
        """
        if not fp.Source:
            return
        fp.Source.touch()
        fp.Document.openTransaction("Recompute points object")
        fp.Document.recompute()
        fp.Document.commitTransaction()
        
    def evaluate(self, fp):
        """evaluate the Source object in mesh workbench"""
        curWB = Gui.activeWorkbench().name()
        Gui.activateWorkbench("MeshWorkbench")
        Gui.activateWorkbench(curWB)
        Gui.Selection.clearSelection()
        Gui.Selection.addSelection(fp.Source)
        Gui.runCommand("Mesh_Evaluation", 0)
        
    def harmonizeNormals(self, fp):
        """harmonize normals of source object"""
        copy = fp.Source.Mesh.copy()
        copy.rebuildNeighbourHood()
        copy.fixIndices()
        copy.harmonizeNormals()
        fp.Document.openTransaction("Harmonize normals")
        fp.Source.Mesh = copy
        fp.Document.commitTransaction()
        
    def flipNormals(self, fp):
        """flip the normals of the source mesh object"""
        copy = fp.Source.Mesh.copy()
        copy.flipNormals()
        fp.Document.openTransaction("Flip normals")
        fp.Source.Mesh = copy
        fp.Document.commitTransaction()
        

class PointsObjectVP:
    """view provider for PointsObject object"""
    def __init__(self, vobj):
        vobj.Proxy = self
    
    def getIcon(self):
        return os.path.join( iconPath , 'CreatePointsObject.svg')
        
    def attach(self, vobj):
        self.Object = vobj.Object
        
    def claimChildren(self):
        return [self.Object.Source]
        
    def setupContextMenu(self, vobj, menu):
        refresh_action = menu.addAction(f"Recompute {vobj.Object.Label}")
        refresh_action.triggered.connect(lambda: vobj.Object.Proxy.refresh(vobj.Object))
        evaluate_action = menu.addAction(f"Evaluate {vobj.Object.Source.Label} for defects")
        evaluate_action.triggered.connect(lambda: vobj.Object.Proxy.evaluate(vobj.Object))
        harmonize_action = menu.addAction(f"Harmonize Normals")
        harmonize_action.triggered.connect(lambda: vobj.Object.Proxy.harmonizeNormals(vobj.Object))
        flip_action = menu.addAction(f"Flip Normals")
        flip_action.triggered.connect(lambda: vobj.Object.Proxy.flipNormals(vobj.Object))

    def onDelete(self, vobj, subelements):
        vobj.Object.Source.ViewObject.DisplayMode = vobj.Object.OriginalDisplayMode
        vobj.Object.Source.ViewObject.Selectable = vobj.Object.OriginalSelectable
        return True
        
    def setEdit(self, vobj, modNum):
        pass

    def __getstate__(self):
        '''When saving the document this object gets stored using Python's json module.\
                Since we have some un-serializable parts here -- the Coin stuff -- we must define this method\
                to return a tuple of all serializable objects or None.'''
        return {"name": self.Object.Name}
 
    def __setstate__(self,state):
        '''When restoring the serialized object from document we have the chance to set some internals here.\
                Since no data were serialized nothing needs to be done here.'''
        self.Object = FreeCAD.ActiveDocument.getObject(state["name"])
        return None


class MeshRemodelCreatePointsObjectCommandClass(object):
    """Create Points Object command"""

    def __init__(self):
        self.mesh = None

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreatePointsObject.svg') ,
            'MenuText': "Create points &object" ,
            'ToolTip' : 
"""Create a parametric PointsObject linked to the selected mesh object.  The points
object serves as a proxy for the mesh object.  We can easily select the points of
the points object whereas selecting directly the points of a mesh object is not 
supported in FreeCAD.  See also the WireFrame object, which serves a similar role,
but that also has selectable edges that can be used in conjunction with the add or 
remove facet tool.

The advantage of the points object over the wireframe object is it is much faster
to recompute when working with complex meshes with many points and facets.  Also,
the way the WireFrame selection of edges works is the 2 vertices of the edge are
added to the points list kept in memory, but the order the vertices are added is
undefined, meaning you have less control over the order of selection for the vertices
when using the wireframe edge selection for adding multiple facets in one go.  For
that operation it is highly recommended to select the vertices.

Tip: Select the vertices in counter-clockwise fashion to get the facet oriented
with the normal side outward when adding new facets to a mesh object.

"""}
 
    def Activated(self):
        doc = self.mesh.Document
        if hasattr(self.mesh,"Mesh"): #might be a points workbench object
            gu.checkComponents(self.mesh.Mesh)
        elif hasattr(self.mesh,"Points") and hasattr(self.mesh.Points,"Points"):
            #for a points cloud we just create the old non-parametric object
            meshpts = self.mesh.Points.Points
            pts = []
            for m in meshpts:
                p = Part.Point(m)
                pts.append(p.toShape())
            doc.openTransaction("Create points object")
            ptsobj = doc.addObject("Part::Feature","PointsObject")
            ptsobj.Shape = Part.Compound(pts)
            pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
            point_size = pg.GetFloat("PointSize",4.0)
            ptsobj.ViewObject.PointSize = point_size
            self.mesh.ViewObject.Visibility = False
            doc.commitTransaction()
            #create new empty mesh
            doc.openTransaction("Make empty mesh")
            meshobj = doc.addObject("Mesh::Feature", self.mesh.Label)
            meshobj.Mesh = Mesh.Mesh()
            meshobj.ViewObject.DisplayMode = "Flat Lines"
            meshobj.ViewObject.Lighting = "One side"
            doc.commitTransaction()
            return

        doc.openTransaction("Create points object")
        fp = doc.addObject("Part::FeaturePython","PointsObject")
        PointsObject(fp, self.mesh)
        PointsObjectVP(fp.ViewObject)
        doc.commitTransaction()
        doc.recompute()
        return

    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelection()
        if len(sel) == 0:
            return False
        if not bool(sel[0].isDerivedFrom("Mesh::Feature") or sel[0].isDerivedFrom("Points::Feature")):
            return False
        self.mesh = sel[0]
        return True

# end create points class
####################################################################################
# Create the Mesh Remodel WireFrame Object

class WireFrameObject(PointsObject):
    def __init__(self, obj, meshObj, className):
        super(WireFrameObject, self).__init__(obj, meshObj, className)
        obj.FlatLines = False
        
    def execute(self, fp):
        if not fp.Source:
            return
        mesh = fp.Source.Mesh.copy()
        lines = []
        for facet in mesh.Facets:
            pt1 = FreeCAD.Vector(facet.Points[0])
            pt2 = FreeCAD.Vector(facet.Points[1])
            pt3 = FreeCAD.Vector(facet.Points[2])
            try:
                line1 = Part.LineSegment(pt1, pt2).toShape()
                line2 = Part.LineSegment(pt2, pt3).toShape()
                line3 = Part.LineSegment(pt3, pt1).toShape()
                lines.extend([line1, line2, line3])
            except:
                pass
            # this works, but takes too long, and if we multiFuse
            # then that resolves the self-intersections, anyway
            # if not gu.hasLine(line1, lines, 1e-6):
            #     lines.append(line1)
            # if not gu.hasLine(line2, lines, 1e-6):
            #     lines.append(line2)
            # if not gu.hasLine(line3, lines, 1e-6):
            #     lines.append(line3)
        if len(lines) > 1 and len(lines) < 10000:    
            fp.Shape = lines[0].multiFuse(lines[1:])
        elif len(lines) > 1 and len(lines) > 10000:
            FreeCAD.Console.PrintWarning(\
"""Compounding WireFrame rather than fusing since there are more than 10,000 edges.
""")
            fp.Shape = Part.Compound(lines)
        else:
            fp.Shape = Part.Shape()
        
class WireFrameObjectVP(PointsObjectVP):
    def __init__(self, vobj):
        super(WireFrameObjectVP, self).__init__(vobj)
            
    def getIcon(self):
        return os.path.join( iconPath , 'CreateWireFrameObject.svg')
    
   
class MeshRemodelCreateWireFrameObjectCommandClass(object):
    """Create WireFrame Object command"""

    def __init__(self):
        self.mesh = None
        self.pb = None
        self.btn = None
        self.bar = None
        self.bCanceled = False

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateWireFrameObject.svg') ,
            'MenuText': "Create Wire&Frame object" ,
            'ToolTip' : \
"""Create a parametric WireFrame object.  The wireframe will be bound to the mesh
object and can be used as a proxy for the mesh object's edges and vertices, which
can be selected for use with the various mesh editing tools in MeshRemodel:

Add or remove facet(s)
Move point
Remove point

Note that the wire frame object does not automatically recompute when changes are
made to the underlying mesh object.  You can do a manual recompute via the wire frame
object's context menu.  Recomputing can take a long time depending on the complexity
of the mesh and it can also be useful to have the points and edges available after
removing them from the mesh object in order to add them back as a means of curing
defects.
"""}
    def Activated(self):
        gu.checkComponents(self.mesh.Mesh)
        doc = self.mesh.Document
        doc.openTransaction("Create WireFrame object")
        fp = doc.addObject("Part::FeaturePython","WireFrameObject")
        WireFrameObject(fp, self.mesh,"WireFrameObject")
        WireFrameObjectVP(fp.ViewObject)
        doc.commitTransaction()
        doc.recompute()
        return
        
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelection()
        if len(sel) == 0:
            return False
        if not sel[0].isDerivedFrom("Mesh::Feature"):
            return False
        self.mesh = sel[0]
        return True
# end create WireFrame class
################################################################################

#MeshBoundaryWires
class MeshRemodelMeshBoundaryWiresCommandClass(object):
    def __init__(self):
        self.mesh = None

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'MeshBoundaryWires.svg'),
            'MenuText': "Makes boundary wire objects from meshes with holes in them" ,
            'ToolTip' : \
"""Makes boundary wire objects from meshes with holes in them.  The wires produced
are not parametric and might or might not be planar.  A mesh without any holes,
meaning no missing facets, will not produce any wires.

This tool can be used to detect holes that might not be easily visible upon brief
inspection.  Can also aid in diagnosing otherwise difficult to find problem areas,
such as 2 points in very close proximity that are causing self-intersections.

If you have trimmed a mesh with a plane in Mesh workbench, then this can help to 
fill the hole that was created by that process.

If a planar wire can be produced, then a mesh face will also be created from it 
and added as a new document object.  You can merge the face with the mesh in the
Mesh workbench.  After merging use the analyze, evaluate, and repair tool to remove
the duplicated points and to check for any additional defects that might exist.
"""} 

    def Activated(self):
        copy = self.mesh.Mesh.copy()
        wires = MeshPart.wireFromMesh(copy)
        FreeCAD.Console.PrintMessage(f"MeshRemodel: {len(wires)} wires created\n")
        doc = self.mesh.Document
        if not wires:
            return
        doc.openTransaction("Mesh boundary wires")
        for wire in wires:
            obj = doc.addObject("Part::Feature", f"{self.mesh.Label}_BoundaryWire")
            obj.Shape = wire
            plane = wire.findPlane()
            if plane:
                face = Part.makeFace(wire,"Part::FaceMakerBullseye")
                mface = MeshPart.meshFromShape(face, 0.25, 50)
                mobj = doc.addObject("Mesh::Feature",f"{self.mesh.Label}_MeshFace")
                mobj.Mesh = mface
        doc.commitTransaction()

    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        self.mesh = None
        sel = Gui.Selection.getSelection()
        self.mesh = sel[0] if len(sel) == 1 else None
        if self.mesh and self.mesh.isDerivedFrom("Mesh::Feature"):
            return True
        return False

################################################################################

#create a simple copy of a mesh object
class MeshRemodelMeshSimpleCopyCommandClass(object):
    """creates a simple copy of a mesh object"""
    
    def __init__(self):
        self.mesh = None

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'MeshSimpleCopy.svg'),
            'MenuText': "Create a simple copy of a mesh object or a new empty mesh" ,
            'ToolTip' : \
"""Create a simple copy of a mesh object.  If no mesh is selected, then a new
empty mesh object is created.
"""}

    def Activated(self):
        doc = self.mesh.Document if self.mesh else FreeCAD.ActiveDocument
        if self.mesh:
            doc.openTransaction("Mesh simple copy")
            duplicate = doc.addObject("Mesh::Feature",self.mesh.Label)
            duplicate.Mesh = self.mesh.Mesh.copy()
            doc.commitTransaction()
        else:
            doc.openTransaction("Make new empty mesh")
            duplicate = doc.addObject("Mesh::Feature","Mesh")
            duplicate.Mesh = Mesh.Mesh()
            doc.commitTransaction() 

    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        self.mesh = None
        sel = Gui.Selection.getSelection()
        self.mesh = sel[0] if len(sel) == 1 else None
        if self.mesh and not self.mesh.isDerivedFrom("Mesh::Feature"):
            return False
        return True    

################################################################################

#offset mesh

class OffsetMesh:
    """offset a mesh by moving its points along their normals"""
    def __init__(self, obj, meshobj=None):
        obj.Proxy = self
        obj.addProperty("App::PropertyLink","Source","OffsetMesh","Source mesh object").Source=meshobj
        obj.addProperty("App::PropertyFloat","Offset","OffsetMesh", \
                "Amount of offset, can be negative value").Offset=2.0
        obj.addProperty("App::PropertyFloat","ZMin","OffsetMesh", \
                "Not really sure what this does yet, it's an option offered by the API").ZMin = 0.0
        obj.addProperty("App::PropertyFloat","ZMax","OffsetMesh",\
                "Not really sure what this does yet, it's an option offered by the API").ZMax = 0.0
        obj.addProperty("App::PropertyBool","MergeSource","OffsetMesh",\
                "If true, merge the original source mesh with the offset one.").MergeSource = False
        obj.addProperty("App::PropertyBool","FlipNormals","OffsetMesh",\
                "Flip normals if true, done before the offset").FlipNormals = False
        obj.addProperty("App::PropertyBool","HarmonizeNormals","OffsetMesh",\
                "Harmonize normals if true, done before the offset").HarmonizeNormals = False
        
    def execute(self, fp):
        if not fp.Source:
            return
        copy = fp.Source.Mesh.copy()
        if fp.FlipNormals:
            copy.flipNormals()
        if fp.HarmonizeNormals:
            copy.harmonizeNormals()
        if fp.Offset != 0.0 and bool(fp.ZMin or fp.ZMax):
            copy.offsetSpecial(fp.Offset, fp.ZMin, fp.ZMax)
        elif fp.Offset != 0.0:
            copy.offset(fp.Offset)
        if fp.MergeSource:
            copy.addMesh(fp.Source.Mesh.copy())
        fp.Mesh = copy

class OffsetMeshVP:        
    def __init__(self, obj):
        obj.Proxy = self
    
    def attach(self, obj):
        self.Object = obj.Object
        
    def updateData(self, fp, prop):
        '''If a property of the handled feature has changed we have the chance to handle this here'''
        # fp is the handled feature, prop is the name of the property that has changed
        pass

    def setDisplayMode(self,mode):
        '''Map the display mode defined in attach with those defined in getDisplayModes.\
                Since they have the same names nothing needs to be done. This method is optional'''
        return mode
 
    def claimChildren(self):
        return [self.Object.Source]
        
    def onDelete(self, vobj, subelements):
        self.Object.Source.Visibility = True
        return True
 
    def onChanged(self, vp, prop):
        '''Here we can do something when a single property got changed'''
        #FreeCAD.Console.PrintMessage("Change property: " + str(prop) + "\n")
        
    def getIcon(self):
        return os.path.join( iconPath , 'OffsetMesh.svg')
    
    def __getstate__(self):
        '''When saving the document this object gets stored using Python's json module.\
                Since we have some un-serializable parts here -- the Coin stuff -- we must define this method\
                to return a tuple of all serializable objects or None.'''
        return {"name": self.Object.Name}
 
    def __setstate__(self,state):
        '''When restoring the serialized object from document we have the chance to set some internals here.\
                Since no data were serialized nothing needs to be done here.'''
        self.Object = FreeCAD.ActiveDocument.getObject(state["name"])
        return None


class MeshRemodelOffsetMeshCommandClass(object):
    """offset a mesh"""
    def __init__(self):
        self.mesh = None

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'OffsetMesh.svg'),
            'MenuText': "Offset a mesh" ,
            'ToolTip' : \
"""Offset a mesh
Offsetting a mesh involves simply moving its points out on their normals (or in 
if the offset is a negative value).  It doesn't always work as expected, but can
at times cure some issues with a mesh, and sometimes it works.

Offset: the amount of offset (in mm)
ZMin: Not sure what this does.  It is an option provided by the API.
ZMax: Not sure what this does.  It is an option provided by the API.
"""}    


    def Activated(self):
        doc = self.mesh.Document
        doc.openTransaction("Offset Mesh")
        fp = doc.addObject("Mesh::FeaturePython","OffsetMesh")
        OffsetMesh(fp, self.mesh)
        OffsetMeshVP(fp.ViewObject)
        self.mesh.ViewObject.Visibility = False
        doc.commitTransaction()
        doc.recompute()


    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelection()
        if len(sel) != 1:
            return False
        self.mesh = sel[0] if sel[0].isDerivedFrom("Mesh::Feature") else None
        if not self.mesh:
            return False
        return True



################################################################################

#duplicate selected facets

class MeshRemodelDuplicateSelectedFacetsCommandClass(object):
    """duplicate selected facets"""
    def __init__(self):
        self.mesh1 = None
        self.mesh2 = None
        self.pts = []
        
    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'DuplicateSelectedFacets.svg'),
            'MenuText': "Duplicated selected facets" ,
            'ToolTip' : \
"""Duplicate selected facets.
Usage: Select 1st the source mesh, then select 2nd the destination mesh.  Finally,
select the facets (by edge or vertex) to be duplicated into the destination mesh.
You can create an empty mesh object to put the duplicates into with the Create 
simple copy tool, which creates an empty mesh when nothing is selected.  Facet 
selections are done on the associated PointsObject or WireFrameObject, representing 
the source mesh.  Box selection (Shift+E) is supported for WireFrameObject types
using edge selection.  (Press Ctrl while using the mouse to do the box selection 
after pressing Shift+E to get into box selection mode.)
"""}

    def Activated(self):
        doc = self.mesh1.Document
        point_indices = gu.getPointIndicesInMesh(self.pts, self.mesh1.Mesh)
        #print (f"point_indices: {point_indices}")
        facet_indices = gu.getFacetIndicesFromPointIndices(point_indices, self.mesh1.Mesh)
        #print(f"facet_indices: {facet_indices}")
        facets = gu.getFacetsFromFacetIndices(facet_indices, self.mesh1.Mesh)
        #print(f"facets = {facets}")
        copy = self.mesh2.Mesh.copy()
        copy.addFacets(facets)
        print(f"Source: {self.mesh1.Label} -> Destination: {self.mesh2.Label} \
{len(facets)} facets duplicated.")
        doc.openTransaction("Duplicate facets")
        self.mesh2.Mesh = copy
        self.mesh2.ViewObject.DisplayMode = "Flat Lines"
        self.mesh2.ViewObject.Lighting = "One side"
        doc.commitTransaction()
       
    # def facetIsInList(self, facet, point_indices):
    #     """all 3 of the facet's points must be in the point_indices"""
    #     for idx in facet.PointIndices:
    #         if not idx in point_indices:
    #             return False
    #     return True

    def getMesh(self, s):
        """get the mesh object from the selection object, can be a direct selection
        or a Points object or WireFrame object, return None if not a mesh"""
        if s.Object.isDerivedFrom("Part::FeaturePython") and \
                bool("PointsObject" in s.Object.Name or "WireFrameObject"\
                in s.Object.Name):
            return s.Object.Source
        if s.Object.isDerivedFrom("Mesh::Feature"):
            return s.Object
        return None

    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getCompleteSelection()

        count = 0
        self.mesh1 = None
        self.mesh2 = None
        self.pts = []
        for idx,s in enumerate(sel):
            #usage: select from mesh first, then select destination mesh, then select
            #points/edges to duplicate from sourc to destination
            if idx == 0:
                self.mesh1 = self.getMesh(s)
            if idx == 1:
                self.mesh2 = self.getMesh(s)
            if idx == 2 and not bool(self.mesh1 and self.mesh2) or self.mesh1 == self.mesh2:
                return False
            if s.SubElementNames and s.SubElementNames[0]:
                sub = s.Object.getSubObject(s.SubElementNames[0])
                if "Vertex" in s.SubElementNames[0]:
                    self.pts.append(sub.Point)
                    count += 1
                if "Edge" in s.SubElementNames[0] or "Face" in s.SubElementNames[0]:
                    for v in sub.Vertexes:
                        if not gu.hasPoint(v.Point, self.pts, 1e-5):
                            self.pts.append(v.Point)
                            count += 1

        if count >= 3:
            return True
        return False

################################################################################
#add or remove a point to or from a mesh

class MeshRemodelRemovePointCommandClass(object):
    """Remove a point from a mesh object"""
    def __init__(self):
        self.mesh = None
        self.pt = None
        self.edge = None

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'RemovePoint.svg') ,
            'MenuText': "Remove a point from (or add a point to) a mesh object" ,
            'ToolTip' : \
"""Remove a point from (or add a point to) a mesh object.
Select the mesh object and a vertex, then run the command
to remove the selected point.  This also removes all facets connected
to that point.

Alt+Click = add a new point to the mesh at the picked point on the selected
edge of a facet.  Edge must be selected from a WireFrame object.  New point will 
go at the picked point location (location on edge where the mouse pointer was when
the edge was selected).  There must be at least one facet containing that edge.
"""}        
    
    def Activated(self):
        modifiers = QtGui.QApplication.keyboardModifiers()
        if modifiers & QtCore.Qt.AltModifier:
            if self.edge:
                pts = [v.Point for v in self.edge.Vertexes]
                pt_indices = gu.getPointIndicesInMesh(pts, self.mesh.Mesh)
                facet_indices = []
                for facet in self.mesh.Mesh.Facets:
                    count = 0
                    for ptIdx in pt_indices:
                        if ptIdx in facet.PointIndices:
                            count += 1
                    if count == 2:
                        facet_indices.append(facet.Index)
                    if len(facet_indices) == 2:
                        break
                if not bool(len(facet_indices) == 1 or len(facet_indices) == 2):
                    FreeCAD.Console.PrintError(\
f"""MeshRemodel Error: Facet not found in {self.mesh.Label}. Select an edge from
 the wireframe object.  The facet(s) will be in the middle of that edge where the
 new point is added.""")
                    return
                facets = gu.getFacetsFromFacetIndices(facet_indices, self.mesh.Mesh)
                copy = self.mesh.Mesh.copy()
                for facet in facets:
                    mid = self.mid
                    pt3_index = [idx for idx in facet.PointIndices if not idx in pt_indices][0]
                    pt3 = copy.Points[pt3_index].Vector
                    copy.splitFacet(facet.Index, mid, pt3)
                self.mesh.Document.openTransaction("Add point")
                self.mesh.Mesh = copy
                self.mesh.touch()
                self.mesh.Document.recompute()
                self.mesh.Document.commitTransaction()
        else:
            if self.pt:
                new_mesh = self.removeFacets(self.mesh.Mesh.copy(), self.pt)
                self.mesh.Document.openTransaction("Remove point")
                self.mesh.Mesh = new_mesh
                self.mesh.Document.commitTransaction()
            else:
                FreeCAD.Console.PrintError("MeshRemodel: No point selected to remove.\n")
            

    def removeFacets(self, mesh, pt):
        """remove the point, and all the facets containing the point
        and return the new mesh"""
        copy = mesh.copy()
        topology = copy.Topology
        # topology is a tuple ([vectorlist],[facetlist])
        # facet list is a list of tuples
        # each tuple of the form (pt1_idx, pt2_idx, pt3_idx)
        idx = None #will be the index of the point in the points
        points = mesh.Points
        for point in points:
            if gu.isSamePoint(pt, point.Vector, 1e-5):
                idx = point.Index
        facets = [tuple(sorted(f)) for f in topology[1]]
        #print(f"facets = {facets}")
        indices = []
        for f_idx,facet_tuple in enumerate(facets):
            if idx in facet_tuple:
                indices.append(f_idx)
        if indices:
            copy.removeFacets(indices)
            return copy
        else:
            FreeCAD.Console.PrintError("MeshRemodel: Unable to remove facets\n")
            return mesh
        
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getCompleteSelection()
        if len(sel) > 2:
            return False
        count = 0
        edgecount = 0
        meshcount = 0
        self.edge = None
        self.pt = None
        self.mid = None
        for s in sel:
            if hasattr(s.Object,"Mesh") or s.Object.isDerivedFrom("Mesh::Feature") :
                meshcount = 1
                self.mesh = s.Object
            if s.Object.isDerivedFrom("Part::FeaturePython") and \
                    bool("PointsObject" in s.Object.Name or "WireFrameObject" in s.Object.Name):
                self.mesh = s.Object.Source if s.Object.Source else self.mesh
                if s.Object.Source:
                    meshcount = 1
            if s.SubElementNames and "Vertex" in s.SubElementNames[0]:
                sub = s.Object.getSubObject(s.SubElementNames[0])
                self.pt = sub.Point
                count += 1
            if s.SubElementNames and "Edge" in s.SubElementNames[0]:
                sub = s.Object.getSubObject(s.SubElementNames[0])
                self.edge = sub
                edgecount += 1
                self.mid = s.PickedPoints[0]
                
        if count == 1 and meshcount == 1 and gu.findPointInMesh(self.mesh.Mesh, self.pt):
            return True
        elif count == 0 and meshcount == 1 and edgecount == 1 and self.edge \
                and len(self.edge.Vertexes) == 2 \
                and gu.findPointInMesh(self.mesh.Mesh, self.edge.Vertex1.Point) != []\
                and gu.findPointInMesh(self.mesh.Mesh, self.edge.Vertex2.Point) != []:
            return True
        return False
        
################################################################################
#move a point in a mesh object

class VectorEditor(QtGui.QDialog):
    def __init__(self, cmd, pt, normal):
        super(VectorEditor, self).__init__()
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowStaysOnTopHint)
        self.ptBase = pt #before moving in normal direction
        self.pt = pt
        self.normalDir = normal.normalize()
        self.finished = False
        self.cmd = cmd
        self.Vertex = None
        self.blockingSignals = False
        self.setWindowTitle("MeshRemodel: Move mesh point")
        self.setGeometry(100, 100, 300, 200)
        self.layout = QtGui.QVBoxLayout(self)

        self.label = QtGui.QLabel()
        self.layout.addWidget(self.label)
        self.label.setText(\
"""Normal is the distance to move the vector in its normal direction.  This overrides
and resets the x,y,z values back to their defaults, and then the normal movement is
applied to that default position, which is the point's original position when the dialog
is first opened.  The normal direction of a mesh point is outward from the mesh.

A red vertex on the screen is there to show the destination.  The Ghost object will be
removed when the dialog is closed.

Note: this operation changes the mesh object topology and can lead to potential issues.

Use Undo to undo the changes if desired.
        
""")

        #create function also adds to layout
        self.spinner_x = self.create_spinner("X:", self.pt.x)
        self.spinner_y = self.create_spinner("Y:", self.pt.y)
        self.spinner_z = self.create_spinner("Z:", self.pt.z)
        self.spinner_normal = self.create_spinner("Normal:", 0.0)
        
        buttonbox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel, self)
        buttonbox.accepted.connect(self.accept)
        buttonbox.rejected.connect(self.reject)
        self.layout.addWidget(buttonbox)
        
        self.doc = self.cmd.mesh.Document
        self.Vertex = self.doc.addObject("Part::Vertex","Ghost")
        self.Vertex.ViewObject.PointColor = (255,0,0)
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        point_size = pg.GetFloat("PointSize",4.0)
        self.Vertex.ViewObject.PointSize = point_size
        self.updateData(self.spinner_normal)
        self.updateVertex()

    def removeVertex(self):
        """removes the Ghost object"""
        if self.Vertex:
            self.doc.removeObject(self.Vertex.Name)
            self.Vertex = None  

    def closeEvent(self, event):
        self.finished = True
        self.removeVertex()
        event.accept()
        
    def reject(self):
        self.finished = True
        self.pt = self.ptBase
        self.removeVertex()
        super().reject()
        
    def accept(self):
        self.finished = True
        self.removeVertex()
        super().accept()
        
    def create_spinner(self, label_text, val):
        def trigger(value,spinner): return lambda value, spinner=spinner: self.spinner_value_changed(value,spinner)
        
        spinner = QtGui.QDoubleSpinBox(self)
        spinner.setObjectName(f"{label_text[:-1]}") #no colon
        spinner.setValue(val)
        spinner.setRange(-1000,1000)
        spinner.setSingleStep(0.1)

        label = QtGui.QLabel(label_text)
        spinner.valueChanged.connect(trigger(val, spinner))
        layout = QtGui.QVBoxLayout()
        layout.addWidget(label)
        layout.addWidget(spinner)

        widget = QtGui.QWidget()
        widget.setLayout(layout)
        self.layout.addWidget(widget)
        return spinner

    def spinner_value_changed(self, val, spinner):
        self.updateData(spinner) #updates the value of self.pt based on spinboxes
        self.updateVertex()
        
    def updateData(self, spinner):
        """update self.pt based on spinner boxes"""
        if self.blockingSignals:
            return
        
        if "Normal" in spinner.objectName():
            self.pt = self.ptBase + self.normalDir * self.Normal
            self.blockingSignals = True
            self.X = self.pt.x
            self.Y = self.pt.y
            self.Z = self.pt.z
            self.blockingSignals = False
        self.pt = FreeCAD.Vector(self.X, self.Y, self.Z)

            
    def updateVertex(self):
        self.Vertex.X = self.X
        self.Vertex.Y = self.Y
        self.Vertex.Z = self.Z
        
    @property
    def X(self):
        return self.spinner_x.value()

    @X.setter
    def X(self,val):
        self.spinner_x.setValue(val)

    @property
    def Y(self):
        return self.spinner_y.value()

    @Y.setter
    def Y(self,val):
        self.spinner_y.setValue(val)
        
    @property
    def Z(self):
        return self.spinner_z.value()

    @Z.setter
    def Z(self,val):
        self.spinner_z.setValue(val)        

    @property
    def Normal(self):
        return self.spinner_normal.value()

    @Normal.setter
    def Normal(self,val):
        self.spinner_normal.setValue(val)  


class MeshRemodelMovePointCommandClass(object):
    """Move a point in a mesh object"""
    def __init__(self):
        self.mesh = None
        self.pts = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'MovePoint.svg') ,
            'MenuText': "Move a point in a mesh object" ,
            'ToolTip' : \
"""
This will relocate a point in a mesh object, moving all connected
facets with it.  Usage:

Select 1 point in the points object, then click the toolbar icon.  In
a dialog you will be able to set the new coordinates for the point.

Select 2 points.  The first point is the point to move, the 2nd point
is the destination for the first point.  Then all duplicated points will
be removed, so this, in effect, would remove the first selected point from
the mesh, but keep the facets that are connected to it by moving them to
the new destination point, assuming the 2nd point is also a point in the
mesh object.  If it's not then the point is simply moved to that location.

Tip: Moving a point often results in folds on surface errors with the mesh
analyze and repair tool.  These errors can often be resolved by converting
the mesh to a shape in Part workbench, and applying a tiny chamfer or fillet
to the affected edges, and then converting back to a mesh.
"""}

    def Activated(self):
        doc = self.mesh.Document
        if len(global_picked) >= 2:
            self.pts = global_picked 
        mesh = self.mesh.Mesh.copy()
        if len(self.pts) == 2:
            destination = self.pts[1]
        else:
            destination = self.getDestination()
        if gu.isSamePoint(self.pts[0], destination, 1e-5):
            return
        idx = gu.findPointInMesh(mesh, self.pts[0])
        if idx == []:
            return
        doc.openTransaction("Move point")
        mesh.setPoint(idx, destination)
        mesh.removeDuplicatedPoints()
        self.checkForDuplicateFacets(mesh)
        self.mesh.Mesh = mesh
        doc.commitTransaction()
        doc.recompute()
        
    def checkForDuplicateFacets(self, mesh):
        num_points = len(mesh.Points)
        num_facets = len(mesh.Facets)
        calculated = num_points * 2 - 4
        if calculated < num_facets:
            mesh.removeFoldsOnSurface()
            FreeCAD.Console.PrintWarning("MeshRemodel: removed folds on surface\n")
        
    def getDestination(self):
        """get the destination for the point from the user"""
        #for now, just return the origin
        idx = gu.findPointInMesh(self.mesh.Mesh, self.pts[0])
        if idx == None:
            return None #calling function returns anyway if idx can't be found
        normal = self.mesh.Mesh.getPointNormals()[idx]
        dlg = VectorEditor(self, self.pts[0], normal)
        dlg.show()
        while not dlg.finished:
            FreeCADGui.updateGui()
        pt = dlg.pt
        dlg.deleteLater()
        return pt
        
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getCompleteSelection()
        if len(sel) == 0 or len(sel) > 3:
            return False
        count = 0
        meshcount = 0
        self.pts = []
        for s in sel:
            if s.SubElementNames and s.SubElementNames[0]:
                sub = s.Object.getSubObject(s.SubElementNames[0])
                if "Vertex" in s.SubElementNames[0]:
                    self.pts.append(sub.Point)
                    count += 1
                if "Edge" in s.SubElementNames[0]:
                    return False
            if s.Object.isDerivedFrom("Part::FeaturePython") and \
                    bool("PointsObject" in s.Object.Name or "WireFrameObject" in s.Object.Name):
                meshcount = 1
                self.mesh = s.Object if hasattr(s.Object,"Mesh") else s.Object.Source
        if bool(count == 1 or count == 2) and meshcount == 1:
            if gu.findPointInMesh(self.mesh.Mesh, self.pts[0]) != []:
                return True
        return False


################################################################################
#add a facet to a mesh

class MeshRemodelAddOrRemoveFacetCommandClass(object):
    """Add a facet(s) to a mesh object"""
    def __init__(self):
        self.mesh = None
        self.pts = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'AddOrRemoveFacet.svg') ,
            'MenuText': "Add or remove facet(s) to or from a mesh object" ,
            'ToolTip' : \
"""Add or remove facet(s) to or from a mesh object.  Select the mesh, select the
points or edges to define the facet(s), and then execute the command.  You can select an 
edge of a WireFrame object in lieu of selecting the end points of that edge or 2 
adjacent edges instead of 3 points or all 3 edges that define the facet.

Click = insert facet
Ctrl+Click = insert facet with flipped normal
Alt+Click = remove facet(s)

Tip: Do not add a flipped normal directly over another facet as this will result
in a defective mesh.  Remove the flipped facet first, with Alt+Click, then add
it back with Ctrl+Click.

Tip: If you have a mesh with one or more missing facets, create a points
object, then select the 3 points on the points object that represent the
missing facet, along with the mesh object to add the facet to the mesh.

If you notice the normal is flipped, remove it with Alt+Click (or Undo),
and then add it back with Ctrl+Click.

If you have multiple facets to make you can select first the anchor point,
which will be common to all the facets, then the remaining points that
will make up the rest of the facets.  These will be added in a loop as:

anchor, pt2, pt3
anchor, pt3, pt4
anchor, pt4, pt5...

Care should be taken to ensure the facets do not intersect previous facets.  Make
the anchor point the that has clear line of sight to the most nearby points.

Note: This changes the mesh object, but you can use Undo to undo the
change if it is not successful or work on a duplicate of the mesh instead
of working on the original.
"""}

    def Activated(self):
        self.mesh.ViewObject.Lighting = "One side"
        doc = self.mesh.Document
        if len(global_picked) > 2:
            self.pts = global_picked 
        modifiers = QtGui.QApplication.keyboardModifiers()
        mesh = self.mesh.Mesh.copy()
        if len(self.pts) == 3:
            if modifiers & QtCore.Qt.AltModifier:
                #remove facet
                removed = self.removeFacet(mesh, self.pts)
                doc.openTransaction("Remove facet")
                self.mesh.Mesh = removed 
                doc.commitTransaction()
                doc.recompute()
                return
            
            if modifiers & QtCore.Qt.ControlModifier:
                mesh.addFacet(self.pts[2],self.pts[1],self.pts[0])
                print(f"MeshRemodel: facet added: {self.pts[2], self.pts[1], self.pts[0]}")
            else:
                mesh.addFacet(self.pts[0],self.pts[1],self.pts[2])
                print(f"MeshRemodel: facet added: {self.pts[0], self.pts[1], self.pts[2]}")
            doc.openTransaction("Add facet")
            self.mesh.Mesh = mesh
            doc.commitTransaction()
        elif QtCore.Qt.AltModifier & modifiers:
            #remove multiple facets
            point_indices = gu.getPointIndicesInMesh(self.pts, mesh)
            facet_indices = gu.getFacetIndicesFromPointIndices(point_indices, mesh)
            mesh.removeFacets(facet_indices)
            print(f"Removing {len(facet_indices)} facets")
            doc.openTransaction("Remove facets")
            self.mesh.Mesh = mesh
            doc.commitTransaction()
            return
        # when more than 3 points are selected we make facets of all of them
        # and add all the facets.  First selected point becomes the anchor
        # point, then we go 1,2,3; 1,3,4; 1;4,5 and so on
        else:
            anchor = self.pts[0]
            for cur,next in zip(self.pts[1:-1], self.pts[2:]):
                if modifiers & QtCore.Qt.ControlModifier:
                    mesh.addFacet(next, cur, anchor)
                    print(f"MeshRemodel: facet added: {next, cur, anchor}")
                else:
                    mesh.addFacet(anchor, cur, next)
                    print(f"MeshRemodel: facet added: {anchor, cur, next}")
            doc.openTransaction("Add facets")
            self.mesh.Mesh = mesh
            doc.commitTransaction()
        doc.recompute()
        
    def removeFacet(self, mesh, pts):
        """remove the facet containing the 3 points and return the new mesh"""
        copy = mesh.copy()
        topology = copy.Topology
        # topology is a tuple ([vectorlist],[facetlist])
        # facet list is a list of tuples
        # each tuple of the form (pt1_idx, pt2_idx, pt3_idx)
        indices = []
        points = mesh.Points
        for point in points:
            for pt in pts:
                if gu.isSamePoint(pt, point.Vector, 1e-4):
                   if not point.Index in indices:
                       indices.append(point.Index)
        #print(f"indices = {indices}")
        facet_tuple = tuple(sorted(indices))
        #print(f"facet_tuple = {facet_tuple}")
        facets = [tuple(sorted(f)) for f in topology[1]]
        #print(f"facets = {facets}")
        if facet_tuple in facets:
            facet_idx = facets.index(facet_tuple)
            print(f"MeshRemodel: removing facet: {copy.Facets[facet_idx]}")
            copy.removeFacets([facet_idx])

            return copy
        else:
            FreeCAD.Console.PrintError("MeshRemodel: Unable to remove facet\n")
            return mesh
        
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getCompleteSelection()
        if len(sel) < 2:
            return False
        count = 0
        meshcount = 0
        self.pts = []
        for s in sel:
            if s.SubElementNames and s.SubElementNames[0]:
                sub = s.Object.getSubObject(s.SubElementNames[0])
                if "Vertex" in s.SubElementNames[0]:
                    self.pts.append(sub.Point)
                    count += 1
                if "Edge" in s.SubElementNames[0] or "Face" in s.SubElementNames[0]:
                    for v in sub.Vertexes:
                        if not gu.hasPoint(v.Point, self.pts, 1e-5):
                            self.pts.append(v.Point)
                            count += 1
            if s.Object.isDerivedFrom("Part::FeaturePython") and \
                    bool("PointsObject" in s.Object.Name or "WireFrameObject"\
                    in s.Object.Name) or s.Object.isDerivedFrom("Mesh::Feature"):
                meshcount = 1
                self.mesh = s.Object if hasattr(s.Object,"Mesh") else s.Object.Source
        if bool(count >= 3) and meshcount == 1:
            return True
        return False

####################################################################################
#Expand a mesh object with a plane

class MeshExpander:
    def __init__(self, obj, mesh=None, base=FreeCAD.Vector(0,0,0), normal=FreeCAD.Vector(0,0,1)):
        obj.Proxy = self
        obj.addProperty("App::PropertyVector","Base","MeshExpander","Base position of the expansion plane").Base=base
        obj.addProperty("App::PropertyVector", "Normal", "MeshExpander", "The normal vector of the expansion plane").Normal=normal
        obj.addProperty("App::PropertyLink","SourceMesh","MeshExpander","The source mesh object").SourceMesh=mesh
        obj.addProperty("App::PropertyFloat","Expansion","MeshExpander","Distance of the expansion in the normal direction of the plane").Expansion=5
        obj.addProperty("App::PropertyFloat","ExpansionRev","MeshExpander","Distance of expansion in the reverse direction of the plane, typically a negative number is needed for this.").ExpansionRev=0
        obj.addProperty("App::PropertyBool","FlipNormals","MeshExpander","Whether to flip the normals of the generated center").FlipNormals = True
        obj.addProperty("App::PropertyBool","ShowCrossSection","MeshExpander","Whether to show the cross-section(s) at the plane, creates new document object(s) on every recompute").ShowCrossSection = False

        
    def execute(self, fp):
        doc = fp.Document
        source = fp.SourceMesh.Mesh.copy()
        # leave source unchanged and always work on a copy of it
        trimmed = source.copy() #top
        trimmed.trimByPlane(fp.Base, fp.Normal.negative())
        trimmedRev = source.copy() #bottom
        trimmedRev.trimByPlane(fp.Base, fp.Normal)
        center = Mesh.Mesh() #expanded center

        if fp.Expansion != 0 or fp.ExpansionRev != 0:
            #pb = gu.MRProgress()
            pts1 = trimmed.Points
            pts2 = trimmedRev.Points
            f1 = trimmed.Facets
            f2 = trimmedRev.Facets
            lines1 = []# lines from facets on the plane
            for f in f1:
                indices = self.facetIsOnPlane(f, fp)
                if len(indices) == 2: #facet has 2 points on the plane
                    t1 = pts1[indices[0]]
                    t2 = pts1[indices[1]]
                    lines1.append(Part.LineSegment(t1.Vector, t2.Vector).toShape())
            #now that we have cross-section made we can separate the top and bottom
            vec = fp.Normal * fp.Expansion
            trimmed.translate(vec.x, vec.y, vec.z)
            rev = fp.Normal * fp.ExpansionRev
            trimmedRev.translate(rev.x, rev.y, rev.z)
            
            # lines need to be sorted to ensure vertex1, vertex2, etc. are in proper order
            # sortEdges() returns a list of lists, with all the lines joined where they belong

            sorted_lines = Part.sortEdges(lines1, 1e-4)
            wires = []
            for srt in sorted_lines:
                fusion = srt[0]
                for line in srt[1:]:
                    fusion = fusion.fuse(line)
                if fp.ShowCrossSection:
                    wires.append(fusion)                    
                    section = doc.addObject("Part::Feature",f"{fp.SourceMesh.Name}_Section")
                    section.Shape = fusion
                mapped = []
                for v in fusion.Vertexes:
                    mapped.append([self.nearestMeshPoint(v.Point, trimmed),
                                self.nearestMeshPoint(v.Point, trimmedRev)])
                #print(f"mapped: {mapped}")
                for cur,next in zip(mapped[:-1],mapped[1:]):
                    p1 = trimmed.Points[cur[0]].Vector
                    p2 = trimmedRev.Points[cur[1]].Vector
                    p3 = trimmed.Points[next[0]].Vector
                    p4 = trimmedRev.Points[next[1]].Vector
                    #print(f"p1,p2,p3,p4: {p1,p2,p3,p4}")
                    center.addFacet(p3, p2, p1)
                    center.addFacet(p2, p3, p4)
                p1 = trimmed.Points[mapped[-1][0]].Vector #last edge
                p2 = trimmedRev.Points[mapped[-1][1]].Vector
                p3 = trimmed.Points[mapped[0][0]].Vector #first edge
                p4 = trimmedRev.Points[mapped[0][1]].Vector
                #print(f"p1,p2,p3,p4: {p1,p2,p3,p4}")
                center.addFacet(p3, p2, p1)
                center.addFacet(p2, p3, p4)
            center.removeDuplicatedPoints()
            if fp.FlipNormals:
                center.flipNormals()

            merged = center.copy()
            merged.addMesh(trimmed)
            merged.addMesh(trimmedRev)
            merged.removeDuplicatedPoints()
            merged.mergeFacets()
            fp.Mesh = merged

                
    def nearestMeshPoint(self, pt, mesh):
        """find the nearest point in the mesh to pt, return that index"""
        least = 10**20
        index = None
        for idx,mp in enumerate(mesh.Points):
            dist = abs(mp.Vector.distanceToPoint(pt))
            if dist < least:
                least = dist
                index = idx
        return index

    def facetIsOnPlane(self, facet, fp):
        """check if 2 points in this facet are on the plane, if so returns a list of
        the point indices"""
        planar = []
        for idx,pt in enumerate(facet.Points):
            p = FreeCAD.Vector(pt[0], pt[1], pt[2])
            dist = abs(p.distanceToPlane(fp.Base, fp.Normal))
            if dist < 1e-5:
                planar.append(facet.PointIndices[idx])
        return planar
        

class MeshExpanderVP:
    def __init__(self, obj):
        obj.Proxy = self
    
    def attach(self, obj):
        self.Object = obj.Object
        
    def updateData(self, fp, prop):
        '''If a property of the handled feature has changed we have the chance to handle this here'''
        # fp is the handled feature, prop is the name of the property that has changed
        pass

    def setDisplayMode(self,mode):
        '''Map the display mode defined in attach with those defined in getDisplayModes.\
                Since they have the same names nothing needs to be done. This method is optional'''
        return mode
 
    def claimChildren(self):
        return [self.Object.SourceMesh]
        
    def onDelete(self, vobj, subelements):
        self.Object.SourceMesh.Visibility = True
        return True
 
    def onChanged(self, vp, prop):
        '''Here we can do something when a single property got changed'''
        #FreeCAD.Console.PrintMessage("Change property: " + str(prop) + "\n")
        
    def getIcon(self):
        return os.path.join( iconPath , 'ExpandMesh.svg')
    
    def __getstate__(self):
        '''When saving the document this object gets stored using Python's json module.\
                Since we have some un-serializable parts here -- the Coin stuff -- we must define this method\
                to return a tuple of all serializable objects or None.'''
        return {"name": self.Object.Name}
 
    def __setstate__(self,state):
        '''When restoring the serialized object from document we have the chance to set some internals here.\
                Since no data were serialized nothing needs to be done here.'''
        self.Object = FreeCAD.ActiveDocument.getObject(state["name"])
        return None


            
class MeshRemodelExpandedMeshCommandClass(object):
    """Expand a mesh object with a plane"""

    def __init__(self):
        self.plane = None
        self.mesh = None

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'ExpandMesh.svg') ,
            'MenuText': "Expand a mesh with a plane" ,
            'ToolTip' : \
"""Expands a mesh with a plane.  This function creates a mesh python object
that can expand itself via its Expansion and ExpansionRev properties.  The
expanded area will be identical to the cross-section at the intersection of
the plane.

Usage:
Create a plane and position it where you want to put the expansion.
Select the plane and the mesh, click the toolbar icon.  Edit the Expansion
and ExpansionRev properties to adjust the level of expansion.
"""}    

    def Activated(self):
        gu.checkComponents(self.mesh.Mesh)
        base = self.plane.Shape.Surface.Position
        normal = self.plane.Shape.Surface.Axis
        doc = self.plane.Document
        fp = doc.addObject("Mesh::FeaturePython","MeshExpander")
        MeshExpander(fp, self.mesh, base, normal)
        MeshExpanderVP(fp.ViewObject)
        fp.ViewObject.Lighting = "One side"
        #fp.ViewObject.Proxy = 0
        self.mesh.Visibility = False
        self.plane.Visibility = False
        doc.recompute()
        
        
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelection()
        if len(sel) != 2:
            return False
        self.mesh = None
        self.plane = None
        for s in sel:
            self.mesh = s if s.isDerivedFrom("Mesh::Feature") else self.mesh
            self.plane = s if s.isDerivedFrom("Part::Plane") or s.isDerivedFrom("PartDesign::Plane")\
                          or bool(s.isDerivedFrom("Part::Feature") and len(s.Shape.Faces) == 1\
                           and s.Shape.Face1.Surface.TypeId == "Part::GeomPlane") else self.plane
        if not bool(self.mesh and self.plane):
            return False
        return True


####################################################################################
# Create the Mesh Cross Sections Object

class MeshRemodelCreateCrossSectionsCommandClass(object):
    """Use Mesh Design workbench Cross-Sections command"""

    def __init__(self):
        self.mesh = None

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateCrossSections.svg') ,
            'MenuText': "Create cross-sections ob&ject..." ,
            'ToolTip' : "Create the cross-sections object \
Convenience link to the Mesh Design workbench cross-sections tool \
(These objects should not be directly used as wires, but rather as references for \
creating wires using the Mesh Remodel workbench as with the Points and WireFrame objects) \
"}
 
    def Activated(self):
        gu.checkComponents(self.mesh.Mesh)
        import MeshPartGui, FreeCADGui
        FreeCADGui.runCommand('MeshPart_CrossSections')
        return
   
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelection()
        if len(sel) == 0:
            return False
        if not sel[0].isDerivedFrom("Mesh::Feature"):
            return False
        self.mesh = sel[0]
        return True

# end open mesh section class

####################################################################################
# Convenience links to some oft-used Part Solid commands: Extrude, Sweep, and Loft

class MeshRemodelPartSolidCommandClass(object):
    """Convenience links to Part Solid commands Extrude, Loft, Sweep"""

    def __init__(self):
        pass

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'PartSolid.svg') ,
            'MenuText': "Part Sol&id" ,
            'ToolTip' : fixTip("Perform a Part Solid command:\n\n\
No Modifier = Extrude\n\
Ctrl + Click = Sweep\n\
Shift + Click = Loft\n\
Alt + Click = Revolution \n\
")}
 
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        selobj = FreeCADGui.Selection.getSelectionEx()
        modifiers = QtGui.QApplication.keyboardModifiers()
        if not modifiers == (QtCore.Qt.ControlModifier & QtCore.Qt.ShiftModifier & QtCore.Qt.AltModifier): #one or more modifiers
            if modifiers == QtCore.Qt.ControlModifier: #sweep
                sections = [obj.Object for obj in selobj if obj.HasSubObjects == False]
                if len(selobj) < 1+len(sections) or not sections: #user did not select an edge to use for the spline
                    FreeCAD.Console.PrintMessage("Running Sweep command using the Gui\n")
                    FreeCAD.Console.PrintMessage("Note: you can bypass the Gui if you select one or more profiles **in the tree view**\n")
                    FreeCAD.Console.PrintMessage("and one or more edges (of a single object) in the 3D view to use as a spine,\n")
                    FreeCAD.Console.PrintMessage("and then run this command.\n")
                    FreeCADGui.runCommand("Part_Sweep",0)
                    return
                if len(selobj) > len(sections) + 1: #need a subshapebinder
                    window = QtGui.QApplication.activeWindow()
                    items=["Make detached (nonparametric) SubShapeBinder","Make synchronized (parametric) SubShapeBinder","Cancel"]
                    item,ok = QtGui.QInputDialog.getItem(window,'Mesh Remodel v'+__version__,
'Part::Sweep requires Spine edges all be from the same object.\n\
If you like a SubShapeBinder can be created for this purpose.\n\
You will need to do the sweep again, this time selecting the \n\
desired edges of the SubShapeBinder to use for the Sweep Spine.\n\n\
But it is generally better to make a wire from the edges and use\n\
the edges of the wire as this ensures they are connected.\n\
\n',items,0,False,windowFlags)
                    if ok and item != items[-1]:
                        doc.openTransaction("Create SubShapeBinder")
                        binder = doc.addObject("PartDesign::SubShapeBinder","Binder")
                        support = []
                        for ff in range(len(sections),len(selobj)):
                            subnames = [elem for elem in selobj[ff].SubElementNames if "Edge" in elem]
                            selobj[ff].Object.ViewObject.Visibility = False
                            for sub in subnames:
                                support.append((selobj[ff].Object,sub))
                        binder.Support = support
                        binder.Fuse = True
                        if item == items[0]: #nonparametric
                            binder.BindMode = "Detached"
                        elif item == items[1]:
                            binder.BindMode = "Synchronized"
                            binder.ClaimChildren = True
                        doc.commitTransaction()
                        doc.recompute()
                        FreeCADGui.Selection.clearSelection()
                        for sec in sections:
                            FreeCADGui.Selection.addSelection(doc.Name,sec.Name)
                    return
                doc.openTransaction("Part Sweep")
                f = doc.addObject("Part::Sweep","Sweep")
                f.Sections = sections
                for sec in f.Sections:
                    sec.ViewObject.Visibility = False
                subnames = [elem for elem in selobj[len(f.Sections)].SubElementNames]
                f.Spine = (selobj[len(f.Sections)].Object, subnames)
                selobj[len(f.Sections)].Object.ViewObject.Visibility = False
                f.Solid = self.checkClosed(selobj[0].Object)
                doc.commitTransaction()
            elif modifiers == QtCore.Qt.ShiftModifier: #loft
                if len(selobj)<2:
                    FreeCAD.Console.PrintMessage("Running Loft command using the Gui\n")
                    FreeCAD.Console.PrintMessage("Note: you can bypass the Gui if you select two or more profiles,\n")
                    FreeCAD.Console.PrintMessage("and then run this command.\n")
                    FreeCADGui.runCommand("Part_Loft",0)
                    return
                sections = [obj.Object for obj in selobj]
                doc.openTransaction("Part Loft")
                f = doc.addObject("Part::Loft","Loft")
                f.Sections = sections
                for sec in f.Sections:
                    sec.ViewObject.Visibility = False
                f.Solid = self.checkClosed(selobj[0].Object)
                doc.commitTransaction()
            elif modifiers == QtCore.Qt.AltModifier: #revolution
                doc.openTransaction("Part Revolve")
                f = doc.addObject("Part::Revolution","Revolve")
                f.Source = selobj[0].Object #profile
                f.Source.ViewObject.Visibility = False
                if len(selobj) >= 2:
                    subnames = [elem for elem in selobj[1].SubElementNames]
                    f.AxisLink = (selobj[1].Object, subnames)
                f.Solid = self.checkClosed(selobj[0].Object)
                doc.commitTransaction()

        else: #extrude
            doc.openTransaction("Part Extrude")
            f = doc.addObject("Part::Extrusion","Extrude")
            f.Solid = self.checkClosed(selobj[0].Object)
            f.Base = selobj[0].Object #profile
            f.Base.ViewObject.Visibility = False
            if len(selobj) >= 2:
                f.DirMode = "Edge"
                subnames = [elem for elem in selobj[1].SubElementNames]
                f.DirLink = (selobj[1].Object, subnames)
            else:
                f.DirMode = self.checkNormal(f.Base)
                f.LengthFwd = 10
            doc.commitTransaction()
        doc.recompute()
        return

    def checkClosed(self,obj):
        import DraftGeomUtils as dgu
        bClosed = False
        if hasattr(obj.Shape,"OuterWire"):
            bClosed = obj.Shape.OuterWire.isClosed()
            if not bClosed:
                FreeCAD.Console.PrintWarning("OuterWire of "+obj.Label+" is not closed.  Setting Solid property to False.\n")
        elif obj.Shape.Wires:
            bClosed = obj.Shape.Wires[0].isClosed()
        else:
            return False #must be point or line
        bPlanar = dgu.isPlanar(obj.Shape)
        if not bPlanar:
            FreeCAD.Console.PrintWarning("Shape is non planar.  Setting Solid to False.\n")
        return bClosed and bPlanar

    def checkNormal(self,obj):
        if hasattr(obj.Shape,"OuterWire"):
            return "Normal"
        return "Custom"

    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        selobj = Gui.Selection.getSelectionEx()
        if len(selobj) < 1:
             return False
        elif hasattr(selobj[0].Object,"Shape"):
                return True
        return False

# end part solid


####################################################################################
# Create the Mesh Remodel Point Object

class MeshRemodelCreatePointObjectCommandClass(object):
    """Create Point Object command"""

    def __init__(self):
        self.obj = None
        self.pts = None
    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreatePointObject.svg') ,
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

class CoplanarPoints:
    def __init__(self,obj):
        obj.addExtension("Part::AttachExtensionPython") #make attachable, why not?
        obj.addProperty("App::PropertyVectorList","Points","CoplanarPoints","Not readonly, but they will be updated with each recompute")
        obj.addProperty("App::PropertyLinkSubList","Trio","CoplanarPoints","3 points that define the plane")
        obj.addProperty("App::PropertyFloatConstraint","Tolerance","CoplanarPoints","Bigger means more points can be included, Zero = include all points").Tolerance = (0.3,0.0,float("inf"),0.1)
        obj.addProperty("App::PropertyLink","BasePointsObject","CoplanarPoints","The base points object from which these coplanar points are selected")
        obj.addProperty("App::PropertyFloat","PointSize","CoplanarPoints","Point size taken from settings")
        obj.addProperty("App::PropertyBool","ExplodeCompound","Triggers","Whether to explode compound of points shapes").ExplodeCompound = False
        obj.addProperty("App::PropertyBool","MakeSketch","Triggers","Whether to make a sketch and add points to it as external geometry links").MakeSketch = False
        obj.addProperty("App::PropertyBool","FlattenToPlane","CoplanarPoints","Flatten all points to the same plane defined by Trio").FlattenToPlane=True
        obj.addProperty("App::PropertyString","Version","CoplanarPoints","Version of MeshRemodel used to create this object").Version = __version__
        obj.setEditorMode("Version",1) #readonly
        obj.Proxy = self
        self.inhibitRecomputes = False

    def makeSketch(self,fp):
        doc = FreeCAD.ActiveDocument
        sketch=doc.addObject("Sketcher::SketchObject","Sketch")
        trio = []
        #for vertName in fp.Trio[0][1]:
        #    trio.append(fp.Trio[0][0].Shape.Vertexes[int(vertName[6:])-1].Point)
        sketch.Support = fp.Trio
        sketch.MapMode = "ThreePointsPlane"
        for ii in range(0,len(fp.Shape.Vertexes)):
            vname = 'Vertex'+str(ii+1)
            sketch.addExternal(fp.Name, vname)

    def explodeCompound(self,fp):
        doc = FreeCAD.ActiveDocument
        doc.openTransaction("explode coplanar points")
        import CompoundTools.Explode
        input_obj = fp
        comp = CompoundTools.Explode.explodeCompound(input_obj)
        input_obj.ViewObject.hide()
        for obj in comp[1]:
            obj.ViewObject.PointSize = fp.PointSize
        doc.recompute()
        doc.commitTransaction()

    def execute(self,fp):
        #FreeCAD.Console.PrintMessage("execute called.\n")
        if self.inhibitRecomputes:
            self.inhibitRecomputes = False;
            return
        doc = FreeCAD.ActiveDocument
        #QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        candidates = []
        if fp.BasePointsObject:
            candidates = fp.BasePointsObject.Shape.Vertexes
        coplanar = []
        trio = []
        for vertName in fp.Trio[0][1]:
            trio.append(fp.Trio[0][0].Shape.Vertexes[int(vertName[6:])-1])
        base,normal = gu.getBaseAndNormal(trio)
        if base is None and normal is None:
            FreeCAD.Console.PrintError("MeshRemodel: Unable to find base/normal\n")
            return
        if fp.Tolerance == 0:
            tolerance = float("inf")
        else:
            tolerance = fp.Tolerance
        for v in candidates:
            if gu.isCoplanar(trio,v.Point,tolerance):
                coplanar.append(Part.Vertex(v.Point))
        coplanar = [Part.Vertex(v) for v in trio] + coplanar
        if hasattr(fp, "FlattenToPlane") and fp.FlattenToPlane:
            base,normal = gu.getBaseAndNormal(coplanar[:3])
            if base is None and normal is None:
                FreeCAD.Console.PrintError("MeshRemodel: Unable to find base/normal\n")
                return
            coplanar = [Part.Vertex(gu.projectVectorToPlane(v.Point, base, normal)) for v in coplanar]
        self.inhibitRecomputes = True
        fp.Points = [v.Point for v in coplanar]
        fp.Shape = Part.makeCompound(coplanar)
        fp.ViewObject.PointSize = fp.PointSize

    def onChanged(self,fp,prop):
        #FreeCAD.Console.PrintMessage(prop+" changed\n")
        self.inhibitRecomputes = False
        if prop == "PointSize":
            self.inhibitRecomputes = True
            fp.ViewObject.PointSize = fp.PointSize
        elif prop == "ExplodeCompound" and fp.ExplodeCompound == True:
            self.inhibitRecomputes = True
            self.explodeCompound(fp)
            fp.ExplodeCompound = False
        elif prop == "MakeSketch" and fp.MakeSketch == True:
            self.inhibitRecomputes = True
            self.makeSketch(fp)
            fp.MakeSketch = False
        elif prop == "Points":
            self.inhibitRecomputes = True


class CoplanarPointsVP:
    """View Provider for Coplanar Points FP object"""
    def __init__(self, obj):
        '''Set this object to the proxy object of the actual view provider'''
        obj.Proxy = self
 
    def attach(self, obj):
        '''Setup the scene sub-graph of the view provider, this method is mandatory'''
        self.Object = obj.Object
 
    def updateData(self, fp, prop):
        '''If a property of the handled feature has changed we have the chance to handle this here'''
        # fp is the handled feature, prop is the name of the property that has changed
        pass
 
    def getDisplayModes(self,obj):
        '''Return a list of display modes.'''
        modes=[]
        modes.append("Flat Lines")
        return modes
 
    def getDefaultDisplayMode(self):
        '''Return the name of the default display mode. It must be defined in getDisplayModes.'''
        return "Flat Lines"
 
    def setDisplayMode(self,mode):
        '''Map the display mode defined in attach with those defined in getDisplayModes.\
                Since they have the same names nothing needs to be done. This method is optional'''
        return mode
 
    def onChanged(self, vp, prop):
        '''Here we can do something when a single property got changed'''
        #FreeCAD.Console.PrintMessage("Change property: " + str(prop) + "\n")

 
    def getIcon(self):
        '''Return the icon in XPM format which will appear in the tree view. This method is\
                optional and if not defined a default icon is shown.'''
        return """
/* XPM */
static char *_631571683001[] = {
/* columns rows colors chars-per-pixel */
"64 64 5 1 ",
"  c black",
". c gray68",
"X c #B9B9B9B9B9B9",
"o c gray78",
"O c none",
/* pixels */
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOO    OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOO .. OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOO .. OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOO    OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOO    OOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOO .. OOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOO .. OOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOO    OOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOO",
"O                                                             OO",
"O                                                             OO",
"O                                                             OO",
"O   OOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOO   OO",
"O   OOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOO   OO",
"O   OOOOOO    OOOOOOOOOO    OOOOOOO    OOOOOOOOOOO    OOOOO   OO",
"O   OOOOOO .. OOOOOOOOOO .. OOOOOOO .. OOOOOOOOOOO .. OOOOO   OO",
"O   OOOOOO .. OOOOOOOOOO .. OOOOOOO X. OOOOOOOOOOO .. OOOOO   OO",
"O   OOOOOO    OOOOOOOOOO    OOOOOOO    OOOOOOOOOOO    OOOOO   OO",
"O   OOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOO   OO",
"O   OOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOO   OO",
"O                                                             OO",
"O                                                             OO",
"O                                                             OO",
"OOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOO    OOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOO .. OOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOO .. OOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOO    OOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO    OOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO .. OOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO .. OOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO    OOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"
};


"""
 
    def __getstate__(self):
        '''When saving the document this object gets stored using Python's json module.\
                Since we have some un-serializable parts here -- the Coin stuff -- we must define this method\
                to return a tuple of all serializable objects or None.'''
        return {"name": self.Object.Name}
 
    def __setstate__(self,state):
        '''When restoring the serialized object from document we have the chance to set some internals here.\
                Since no data were serialized nothing needs to be done here.'''
        self.Object = FreeCAD.ActiveDocument.getObject(state["name"])
        return None


class MeshRemodelCreateCoplanarPointsObjectCommandClass(object):
    """Create coplanar points object from 3 selected points"""

    def __init__(self):
        self.pts = []
        self.obj = None #original points object
        self.vertexNames = [] #for attaching a part::plane

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateCoplanar.svg') ,
            'MenuText': "Create copla&nar points object" ,
            'ToolTip' : """\
Filters the vertices in the selected object based on the first 3 vertices
selected from that object.  Only those vertices that are on the plane (within
the given tolerance set in the Tolerance property) are put into the coplanar
points object.

See also MeshRemodel settings to set the default Tolerance.  Tolerance is the
shortest distance from the point to the plane.  (Formerly it was something else.)

If FlattenToPlane = True, then the points are forced to the same plane even
if the original vertices are not exactly on the plane.

Tip: If you are using this to remodel a mesh, then leave FlattenToPlane True,
but if you are adding facets to a mesh, set the FlattenToPlane property to 
False because you want them referencing the exact positions of the points in
the original mesh object.
"""}

    def Activated(self):
        if len(global_picked) == 3:
            self.pts = global_picked #use preselect-picked points
        if gu.isColinear(self.pts[0],self.pts[1],self.pts[2]):
            FreeCAD.Console.PrintError('Please select 3 non-colinear points in the plane\n')
            return
        doc = FreeCAD.ActiveDocument
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        point_size = pg.GetFloat("PointSize",4.0)
        coplanar_tolerance = pg.GetFloat("CoplanarTolerance", .00001)
        doc.openTransaction("Create coplanar")
        cp = doc.addObject("Part::FeaturePython","MR_Coplanar_Points")
        CoplanarPoints(cp)
        CoplanarPointsVP(cp.ViewObject)
        cp.BasePointsObject = self.obj
        cp.PointSize = point_size
        cp.Tolerance = coplanar_tolerance
        cp.Trio = (self.obj,self.vertexNames)
        cp.BasePointsObject.ViewObject.Visibility = False
        doc.commitTransaction()
        FreeCADGui.Selection.clearSelection()
        FreeCADGui.Selection.addSelection(doc.Name,cp.Name)
        doc.recompute()
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
            self.vertexNames = sel[0].SubElementNames
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
        return {'Pixmap'  : os.path.join( iconPath , 'CreateLine.svg') ,
            'MenuText': "Create &line" ,
            'ToolTip' : fixTip("Create a line from 2 selected points\n\
(Ctrl+Click to add midpoint)\n\
(Ctrl+Shift+Click for only midpoint)")}
 
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
        return {'Pixmap'  : os.path.join( iconPath , 'CreatePolygon.svg') ,
            'MenuText': "Create &Polygon" ,
            'ToolTip' : fixTip("\
Create a Polygon from 3 or more selected points or 2 or more selected edges\n\
Might not always be coplanar, consider using links to external geometry in a sketch\n\
Do **not** attempt to mix selected edges and selected points, should be all edges\n\
or all points, but not a combination of the 2 object types\n\
(Makes individual lines, use Create wire to connect into a single wire object.)\n\
(Shift+Click to not close polygon) -- but selected edges never close unless connected\n\
(Alt+Click to sort selected points)\n\
")}
 
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
        return {'Pixmap'  : os.path.join( iconPath , 'CreateBSpline.svg') ,
            'MenuText': "Create &BSpline" ,
            'ToolTip' : fixTip("Create a BSPline from 3 or more selected points\n\
(Ctrl+Click to not flatten bspline)\n\
(Shift+Click to not close bspline)\n\
(Alt+Click to sort selected points)")}
 
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
        is_periodic = not modifiers & QtCore.Qt.ShiftModifier
        do_sort = modifiers & QtCore.Qt.AltModifier
        do_flatten = not modifiers & QtCore.Qt.ControlModifier
        if do_flatten:
            self.pts = gu.projectPointsToPlane(self.pts)
        if do_sort:
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

####################################################################################
# Flatten an existing Draft BSpline object to a plane defined by first, last, and a
# point from the middle of the points list

class MeshRemodelFlattenDraftBSplineCommandClass(object):
    """Flatten an existing Draft BSpline"""

    def __init__(self):
        self.spline = None

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'FlattenDraftBSpline.svg') ,
            'MenuText': "Flatt&en BSpline" ,
            'ToolTip' : fixTip("Flatten an existing Draft BSpline\n")}

    def Activated(self):
        doc = FreeCAD.ActiveDocument
        pg = FreeCAD.ParamGet("User parameter:Plugins/MeshRemodel")
        line_width = pg.GetFloat("LineWidth",5.0)
        point_size = pg.GetFloat("PointSize",4.0)
        doc.openTransaction("Flatten Draft BSpline")
        projected = gu.projectPointsToPlane(self.spline.Points)
        new_spline = Part.BSplineCurve()
        new_spline = Draft.makeBSpline(projected, closed=self.spline.Closed)
        new_spline.Label = f"{self.spline.Label}_flattened"
        new_spline.ViewObject.LineWidth = line_width
        self.spline.ViewObject.Visibility = False
        doc.recompute()
        doc.commitTransaction()
        return

    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelection()
        if len(sel) != 1:
            return False
        self.spline = sel[0]
        if not hasattr(self.spline,"Points"):
            return False
        if not "BSpline" in self.spline.Name:
            return False
        return True

# end create BSpline class
################################################################################

# rotate an object in place
class MeshRemodelRotateObjectCommandClass(object):
    """Rotate object about a given axis"""
    
    def __init__(self):
        self.obj = None
        self.normal = None
        self.center = None
        self.angle = None
        
    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'RotateObject.svg') ,
            'MenuText': "Rotate object" ,
            'ToolTip' : \
"""Rotate an object about a selected subobject

Click = rotate 1 degree
Ctrl+Click = rotate 0.1 degree
Shift+Click = rotate 15 degrees
Alt+Click = rotate in opposite direction
Ctrl+Shift+Click = enter custom angle in dialog

Supported selections:

Planar face or circle edge: axis = normal
Line segment edge: axis =  line direction
2 vertices: axis = direction from one to the other
3 vertices: axis = 3 point plane, center = 2nd point, angle = 1st to 3rd point
2 faces/circles: axis,angle = axis,angle needed to bring them into coplanarity, 
if they are already planar, then Ctrl+Shift+Click, enter 180 degrees to flip.

If you are working in Part Design and wish to rotate a Part Design Feature, such
as a Pad or Pocket, then usually the underlying sketch is what should be rotated,
but if it is attached, then the attachment will override the rotation on recompute.
To move the entire body, set the Body's View tab DisplayModeBody property to Tip
temporarily, and then back to Through after completing the rotation.
"""}
            
    def Activated(self):
        #FreeCAD.Console.PrintMessage(f"obj: {self.obj.Label}, normal: {self.normal}\n")
        angle = self.angle if self.angle != None else 1
        modifiers = QtGui.QApplication.keyboardModifiers()
        if modifiers & QtCore.Qt.ControlModifier:
            angle *= 0.1
        if modifiers & QtCore.Qt.ShiftModifier:
            angle *= 15.0
        if modifiers & QtCore.Qt.AltModifier:
            angle *= -1.0
        if modifiers & QtCore.Qt.ControlModifier and modifiers & QtCore.Qt.ShiftModifier:
            angle = self.angle if self.angle != None else 1
            angle = gu.getFloatFromUser("Rotate", "Enter value for angle:", value=angle, step=0.1)
            if angle == None:
                return #user cancel
        doc = self.obj.Document
        doc.openTransaction("Rotate object")
        plm = self.obj.Placement
        plm.rotate(self.center, self.normal, angle, True)
        self.obj.Placement = plm
        
        if self.obj.isDerivedFrom("PartDesign::Feature"):
            if self.obj._Body:
                if self.obj._Body.ViewObject.DisplayModeBody != "Tip":
                    FreeCAD.Console.PrintWarning(\
f"""MeshRemodel Rotate Object: {self.obj.Label} is a Part Design feature and its
View tab DisplayModeBody is not in Tip mode.  Did you want to rotate the body or
the feature?  Feature placement is usually defined by sketch placement or other
attachment.
""")

        doc.commitTransaction()
        
    def IsActive(self):
        self.obj = None
        self.normal = None
        self.center = None
        self.angle = None
        sel = FreeCADGui.Selection.getCompleteSelection()
        if len(sel) == 0 or len(sel) > 3:
            return False
        if len(sel) == 2:
            pts = []
            if bool("Edge" in sel[0].SubElementNames[0] or "Face" in sel[0].SubElementNames[0]) \
                    and bool("Edge" in sel[1].SubElementNames[0] or "Face" in sel[1].SubElementNames[0]):
                sub1 = sel[0].Object.getSubObject(sel[0].SubElementNames[0])
                sub2 = sel[1].Object.getSubObject(sel[1].SubElementNames[0])
                if hasattr(sub1,"Surface") or \
                        bool(hasattr(sub1,"Curve") and sub1.Curve.TypeId == "Part::GeomCircle"):
                    self.center = sub1.Curve.Center
                    axis = sub1.Curve.Axis if hasattr(sub1,"Curve") else sub1.Surface.Axis
                    if hasattr(sub2,"Surface") or \
                            bool(hasattr(sub2,"Curve") and sub2.Curve.TypeId == "Part::GeomCircle"):
                        axis2 = sub2.Curve.Axis if hasattr(sub2,"Curve") else sub2.Surface.Axis
                        self.normal = axis.cross(axis2)
                        try:
                            self.angle = math.degrees(math.acos(axis.dot(axis2)))
                        except:
                            #same normal, so let's tweak slightly to allow for flipping 180 degrees
                            axis.x = axis.x + 1e-5
                            axis.y = axis.y - 1e-5
                            axis.z = axis.z + 1e-5
                            self.angle = math.degrees(math.acos(axis.dot(axis2)))
                        self.obj = sel[0].Object
                        return True
                if not gu.hasPoint(sub1.Vertex1.Point, pts, 1e-6):
                    pts.append(sub1.Vertex1.Point)
                if not gu.hasPoint(sub1.Vertex2.Point, pts, 1e-6):
                    pts.append(sub1.Vertex2.Point)
                if not gu.hasPoint(sub2.Vertex1.Point, pts, 1e-6):
                    pts.append(sub2.Vertex1.Point)
                if not gu.hasPoint(sub2.Vertex2.Point, pts, 1e-6):
                    pts.append(sub2.Vertex2.Point)
                if len(pts) == 3:
                    base,normal = gu.getBaseAndNormal([Part.Vertex(p) for p in pts])
                    if normal is None:
                        if hasattr(sub1,"Curve") and hasattr(sub1.Curve,"XAxis"):
                            normal = sub1.Curve.XAxis
                        else:
                            return False
                    self.normal = normal
                    base = sub1.Curve.intersect(sub2.Curve)[0]
                    self.center = FreeCAD.Vector(base.X, base.Y, base.Z)
                    self.obj = sel[0].Object
                    self.angle = math.degrees(sub1.Curve.Direction.getAngle(sub2.Curve.Direction))
                    return True
            elif "Vertex" in sel[0].SubElementNames[0] and "Vertex" in sel[1].SubElementNames[0]:
                sub1 = sel[0].Object.getSubObject(sel[0].SubElementNames[0])
                sub2 = sel[1].Object.getSubObject(sel[1].SubElementNames[0])
                self.normal = sub1.Point - sub2.Point
                self.center = sub1.Point
                self.obj = sel[0].Object
                return True
        elif len(sel) == 3:
            if "Vertex" in sel[0].SubElementNames[0] and "Vertex" in sel[1].SubElementNames[0]\
                    and "Vertex" in sel[2].SubElementNames[0]:
                verts = [s.Object.getSubObject(s.SubElementNames[0]) for s in sel]
                base,normal = gu.getBaseAndNormal(verts)
                if normal is None:
                    return False
                self.normal = normal
                self.center = verts[1].Point
                dir1 = verts[0].Point - verts[1].Point
                dir2 = verts[2].Point - verts[1].Point
                self.angle = math.degrees(dir2.getAngle(dir1))
                self.obj = sel[0].Object
                return True
            elif bool("Edge" in sel[0].SubElementNames[0] or "Face" in sel[0].SubElementNames[0]) \
                    and "Vertex" in sel[1].SubElementNames[0] and "Vertex" in sel[2].SubElementNames[0]:
                sub1 = sel[0].Object.getSubObject(sel[0].SubElementNames[0])
                if hasattr(sub1,"Curve") and sub1.Curve.TypeId == "Part::GeomCircle":
                    self.normal = sub1.Curve.Axis
                    self.center = sub1.Curve.Center
                elif hasattr(sub1,"Curve") and sub1.Curve.TypeId == "Part::GeomLine":
                    self.normal = sub1.Curve.Direction
                    self.center = sub1.Vertex1.Point
                elif hasattr(sub1,"Surface") and sub1.Surface.TypeId == "Part::GeomPlane":
                    self.normal = sub1.Surface.Axis
                    self.center = sub1.CenterOfGravity
                sub2 = sel[1].Object.getSubObject(sel[1].SubElementNames[0])
                sub3 = sel[2].Object.getSubObject(sel[2].SubElementNames[0])
                self.obj = sel[0].Object
                dir1 = sub1.CenterOfGravity - sub2.Point
                dir2 = sub1.CenterOfGravity - sub3.Point
                self.angle = math.degrees(dir2.getAngle(dir1))
                return True
        elif len(sel) == 1:        
            self.obj = sel[0].Object
            if not hasattr(self.obj, "Shape"):
                return False
            sub = sel[0].SubObjects[0]
            if hasattr(sub,"Curve"):
                if sub.Curve.TypeId == "Part::GeomLine":
                    self.normal = sub.Curve.Direction
                    self.center = sub.CenterOfGravity
                    return True
                if sub.Curve.TypeId == "Part::GeomCircle":
                    self.normal = sub.Curve.Axis
                    self.center = sub.Curve.Center
                    return True
                
            if hasattr(sub,"Surface"):
                if not sub.Surface.isPlanar():
                    cog = sub.CenterOfGravity
                    pr = sub.ParameterRange
                    u = (pr[0] + pr[1]) / 2.0
                    v = (pr[2] + pr[3]) / 2.0
                    self.normal = sub.valueAt(u,v) - cog
                    self.center = sub.CenterOfGravity
                    return True
            shp = sub if sub else self.obj.Shape
            try:
                self.normal = shp.normalAt(0,0)
                self.center = shp.CenterOfGravity
                return True
            except:
                try:
                    self.normal = shp.normalAt(0)
                    self.center = shp.CenterOfGravity
                    return True
                except:
                    return False
                return False

################################################################################

# Move an object in the axial direction
class MeshRemodelMoveAxialCommandClass(object):
    """Move an object in the normal direction, if a normal can be found"""
    
    def __init__(self):
        self.obj = None
        self.normal = None
        self.distance = None
        
    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'MoveAxial.svg') ,
            'MenuText': "Move A&xial" ,
            'ToolTip' :\
"""
Move an object in the axial (normal) direction of selected subobject(s)
(Edits object's Placement property)
If subobject is a planar face, direction is normal to that face
If subobject is a non planar face, direction is normal to the face center
If subobject is a line segment, direction is the direction of the line

2 subobjects support:

If 2 vertices are selected the direction is the direction from the first to the 
second and the distance will be the distance between the 2 points.  The first
selected vertex belongs to the object that will be moved if the vertices are from
separate objects.

2 faces/circles: direction/distance = center to center.

Click = move 1 mm
Ctrl+Click = move 0.1 mm
Shift+Click = move 10 mm
Alt+Click = move in opposite direction
Ctrl+Shift+Click = enter custom distance in dialog

If you are working in Part Design and wish to move a Part Design Feature, such
as a Pad or Pocket, then usually the underlying sketch is what should be rotated,
but if it is attached, then the attachment will override the move on recompute.
To move the entire body, set the Body's View tab DisplayModeBody property to Tip
temporarily, and then back to Through after completing the move.
"""}
            
    def Activated(self):
        #FreeCAD.Console.PrintMessage(f"obj: {self.obj.Label}, normal: {self.normal}\n")
        distance = 1 if self.distance == None else self.distance
        modifiers = QtGui.QApplication.keyboardModifiers()
        if modifiers & QtCore.Qt.ControlModifier:
            distance *= 0.1
        if modifiers & QtCore.Qt.ShiftModifier:
            distance *= 10.0
        if modifiers & QtCore.Qt.AltModifier:
            distance *= -1.0
        if modifiers & QtCore.Qt.ControlModifier and modifiers & QtCore.Qt.ShiftModifier:
            distance = self.distance if self.distance != None else 1
            distance = gu.getFloatFromUser("Move Axial", "Enter value for distance:", value=distance, step=0.1)
            if distance == None:
                return #user cancel
        doc = self.obj.Document
        doc.openTransaction("MR_Move")
        self.obj.Placement = FreeCAD.Placement(self.obj.Placement.Base 
                        + self.normal.normalize() * distance, self.obj.Placement.Rotation)
        if self.obj.isDerivedFrom("PartDesign::Feature"):
            if self.obj._Body:
                if self.obj._Body.ViewObject.DisplayModeBody != "Tip":
                    FreeCAD.Console.PrintWarning(\
f"""MeshRemodel Move Axial: {self.obj.Label} is a Part Design feature and its
View tab DisplayModeBody is not in Tip mode.  Did you want to move the body or
the feature?  Feature placement is usually defined by sketch placement or other
attachment.
""")

        doc.commitTransaction()
        
    def IsActive(self):
        sel = FreeCADGui.Selection.getCompleteSelection()
        self.distance = None
        if len(sel) > 2 or len(sel) == 0:
            return False
        self.obj = sel[0].Object
        if not hasattr(self.obj, "Shape"):
            return False
        if len(sel) == 2:
            subs = []
            centers = []
            for s in sel:
                subs.append(s.Object.getSubObject(s.SubElementNames[0]))
                centers.append(subs[-1].CenterOfGravity)
            if len(subs) == 2:
                self.normal = centers[1] - centers[0]
                self.distance = centers[1].distanceToPoint(centers[0])
                self.obj = sel[0].Object
                return True
                
        sub = sel[0].SubObjects[0]
        if hasattr(sub,"Curve"):
            if sub.Curve.TypeId == "Part::GeomLine":
                self.normal = sub.Curve.Direction
                return True
            if sub.Curve.TypeId == "Part::GeomCircle":
                self.normal = sub.Curve.Axis
                return True
                
        if hasattr(sub,"Point"):
            self.normal = sub.CenterOfGravity - sel[0].Object.Shape.CenterOfGravity
            return True
        
        if hasattr(sub,"Surface"):
            if not sub.Surface.isPlanar():
                cog = sub.CenterOfGravity
                pr = sub.ParameterRange
                u = (pr[0] + pr[1]) / 2.0
                v = (pr[2] + pr[3]) / 2.0
                self.normal = sub.valueAt(u,v) - cog
                return True
        shp = sub if sub else self.obj.Shape
        try:
            self.normal = shp.normalAt(0,0)
            return True
        except:
            try:
                self.normal = shp.normalAt(0)
                return True
            except:
                return False
            return False
        return False

# end Move Axial command class
################################################################################

class BlockSelectorDialog(QtGui.QDialog):
    def __init__(self, blocks, parent=None):
        def trigger(idx): 
            return lambda: self.on_block_selected(idx)
        super(BlockSelectorDialog, self).__init__(parent)
        self.blocks = blocks
        self.selected_index = None

        self.setWindowTitle('Go back to a previous selection')
        self.setWindowFlags(QtCore.Qt.Dialog | QtCore.Qt.WindowCloseButtonHint)
        self.setModal(True)

        layout = QtGui.QVBoxLayout()

        # Create a scroll area
        scroll_area = QtGui.QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_content = QtGui.QWidget()
        scroll_layout = QtGui.QVBoxLayout(scroll_content)
        scroll_layout.addWidget(QtGui.QLabel("Tooltip for each button shows the selection objects\nDeleted objects will not appear."))

        # Create buttons for each block
        for i, block in enumerate(blocks):
            button = QtGui.QPushButton(f'Block {i + 1}')
            tooltip_text = "\n".join(block)
            button.setToolTip(tooltip_text)
            button.clicked.connect(trigger(i+1))
            scroll_layout.addWidget(button)

        scroll_content.setLayout(scroll_layout)
        scroll_area.setWidget(scroll_content)

        layout.addWidget(scroll_area)
        self.setLayout(layout)

    def on_block_selected(self, index):
        self.selected_index = index
        self.accept()  # Close the dialog and set the result to Accepted


# go back selection
class MeshRemodelGoBackSelectionCommandClass(object):
    """Provide a means of returning to previous selection states"""
    def __init__(self):
        self.num_lines_to_remove = 0
        self.pc = None
    
    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'GoBackSelection.svg') ,
            'MenuText': "Go back selection" ,
            'ToolTip' :\
"""
Return to the previous selection state.

The function makes use of the commented Gui.Selection lines from the python 
console.  Each time a selection changes it gets recorded as a comment in the
python console.  This function parses that history and uses that to return to
the previous selection state.

In the dialog you will find a column of buttons labeled Block 1, Block 2, etc.
Each button has a tooltip that will display the selection commands that will be
invoked when the button is clicked.  Only the selections made previously in
documents with the same name as the current active document will be available.

"""}

    def Activated(self):
        mw = Gui.getMainWindow()
        self.pc = mw.findChild(QtGui.QPlainTextEdit,"Python console")
        txt = self.pc.toPlainText()
        txtLines = txt.split('\n')
        lines = [line[6:] for line in txtLines \
                if line.startswith(">>> # Gui.Selection")]
        blocks = self.make_blocks(lines)
        blocks = self.filter_blocks(blocks, FreeCAD.ActiveDocument.Name)
        if not blocks:
            FreeCAD.Console.PrintError("MeshRemodel: Not enough selection history\n")
            return
        which = self.ask_user_which_block(blocks)
        if which == -1:
            return

        if len(blocks) >= which:
            cmd = "\n".join(blocks[which-1])
            self.num_lines_to_remove = len(blocks[which -1])*2+2
            Gui.Selection.clearSelection()
            Gui.doCommand(cmd)
            timer = QtCore.QTimer()
            timer.singleShot(250, self.remove_lines)
        else:
            FreeCAD.Console.PrintError("MeshRemodel: Not enough selection history\n")
        
    def remove_lines(self):
        """remove the new selection lines we just added when we selected the items in the block"""
        n = self.num_lines_to_remove
        text_to_keep = self.pc.toPlainText().splitlines()[:-n]
        self.pc.clear()
        self.pc.setPlainText('\n'.join(text_to_keep))
        cursor = self.pc.textCursor()
        cursor.movePosition(cursor.End)
        self.pc.setTextCursor(cursor)

    def ask_user_which_block(self, blocks):
        """ask user which block to go back to in a dialog, return index into blocks"""

        dialog = BlockSelectorDialog(blocks)
        result = dialog.exec_()
        
        if result == QtGui.QDialog.Accepted:
            print(f"returning {dialog.selected_index}")
            return dialog.selected_index
        else:
            return -1
            
    def filter_blocks(self, blocks, doc_name):
        """we only want the blocks that reference the current active document 
        and that reference existing objects, not deleted ones
        """
        filtered = []
        for block in blocks:
            if any(doc_name in s for s in block):
                filtered.append(block)
                
        filtered2 = []       
        for block in filtered:
            filtered_block = []
            for s in block:
                if doc_name in s:
                    # Extract the object name, assuming it is the argument following the document name
                    parts = s.split(',')
                    if len(parts) > 1:
                        obj_name = parts[1].strip().strip("'")
                        if self.ifExists(doc_name, obj_name):
                            filtered_block.append(s)
                else:
                    filtered_block.append(s)
    
            if filtered_block:
                filtered2.append(filtered_block)
        
        return filtered2

    def ifExists(self, doc_name, obj_name):
        """check if object still exists"""
        doc = FreeCAD.getDocument(doc_name)
        obj = doc.getObject(obj_name)
        if obj:
            return True
        else:
            return False

    def make_blocks(self, lines):
        """convert the lines into blocks of selections in between the clear selection calls"""
        result = []
        current_block = []

        bSkipping = False #skip MeshRemodelGoBackSelection blocks
        for item in reversed(lines):
            if item.startswith("Gui.Selection.clearSelection("):
                if current_block:
                    rev = reversed(current_block)
                    cb = [b for b in rev]
                    result.append(cb)
                    current_block = []
            else:
                if not item in current_block:
                    current_block.append(item)

        if current_block:
            if not result or bool(current_block != result[-1]):
                result.append(current_block)
        # remove duplicates
        result = [lst for i, lst in enumerate(result) if lst not in result[:i]]
        return result

    def IsActive(self):
        if FreeCAD.ActiveDocument:
            return True
        return False
################################################################################
#SubObjectLoft

class ParametricLine:
    def __init__(self, obj, sub1=None, sub2=None):
        obj.Proxy = self
        obj.addProperty("App::PropertyLinkSub", "Link1", "ParametricLine").Link1 = sub1
        obj.addProperty("App::PropertyLinkSub", "Link2", "ParametricLine").Link2 = sub2
        obj.addProperty("App::PropertyDistance","Link1Extend","ParametricLine","Extend line before Link1")
        obj.addProperty("App::PropertyDistance","Link2Extend","ParametricLine","Extend line after Link2")
        obj.addProperty("App::PropertyString","Info","ParametricLine", "Status information about line")
        
    def execute(self, fp):
        p1 = getattr(fp.Link1[0].Shape, fp.Link1[1][0]).Point if fp.Link1 else None
        p2 = getattr(fp.Link2[0].Shape, fp.Link2[1][0]).Point if fp.Link2 else None
        if not bool(p1 and p2):
            fp.Shape = Part.Shape()
            return
        dir = (p2 - p1).normalize()
        p1 = p1 + (dir * fp.Link1Extend)
        p2 = p2 + (dir * fp.Link2Extend)
        
        fp.Info = f"Line Length: {p1.distanceToPoint(p2)}\nLineDirection: {dir}"
        
        if p1 == p2:
            fp.Shape = Part.Vertex(p1)
            return
            
        line = Part.LineSegment(p1,p2)  
        fp.Shape = line.toShape()
        
class ParametricLineVP:
    def __init__(self, obj):
        obj.Proxy = self
        
    def attach(self, obj):
        self.Object = obj
        
    def __getstate__(self):
        '''When saving the document this object gets stored using Python's json module.\
                Since we have some un-serializable parts here -- the Coin stuff -- we must define this method\
                to return a tuple of all serializable objects or None.'''
        return None
    
    def __setstate__(self,state):
        '''When restoring the serialized object from document we have the chance to set some internals here.\
                Since no data were serialized nothing needs to be done here.'''
        return None
    
    def getDisplayModes(self,obj):
        '''Return a list of display modes.'''
        modes=[]
        modes.append("Flat Lines")
        return modes

    def getIcon(self):
        
        icon = """        
/* XPM */
static char *dummy[]={
"64 64 5 1",
"b c black",
"a c #55ff00",
"c c #55ff7f",
"# c #bcbcbc",
". c None",
"................................................................",
"................................................................",
"................................................................",
"..........#####.................................................",
".........#######................................................",
"........###aaaa##...............................................",
".......###aaaaa###..............................................",
".......##aaaaaaa##..............................................",
".......##aaaaaaa##..............................................",
".......##aaaaaaa###.............................................",
".......##aaaaaaa####............................................",
"........##aaaaa##b###...........................................",
".........#######bbb###..........................................",
"..........######bbbb##..........................................",
"..............###bbbb##.........................................",
"...............##bbbbb##........................................",
"................##bbbb###.......................................",
".................##bbbb##.......................................",
".................###bbbb##......................................",
"..................##bbbbb##.....................................",
"...................##bbbb###....................................",
"....................##bbbb###...................................",
"....................###bbbb##...................................",
".....................###bbbb##..................................",
"......................##bbbbb##.................................",
".......................##bbbb###................................",
"........................##bbbb###...............................",
"........................###bbbb##...............................",
".........................###bbbb##..............................",
"..........................##bbbbb##.............................",
"...........................##bbbb###............................",
"............................##bbbb##............................",
"............................###bbbb##...........................",
".............................##bbbbb##..........................",
"..............................##bbbb###.........................",
"...............................##bbbb###........................",
"...............................###bbbb##........................",
"................................###bbbb##.......................",
".................................##bbbbb##......................",
"..................................##bbbb###.....................",
"...................................##bbbb###....................",
"...................................###bbbb##....................",
"....................................###bbbb##...................",
".....................................##bbbbb##..................",
"......................................##bbbb###.................",
".......................................##bbbb##.................",
".......................................###bbbb##................",
"........................................##bbbbb##...............",
".........................................##bbbb###..............",
"..........................................##bbbb#####...........",
"..........................................###bb#######..........",
"...........................................######cccc##.........",
"............................................####ccccc###........",
".............................................##ccccccc##........",
".............................................##ccccccc##........",
".............................................##ccccccc##........",
".............................................##ccccccc##........",
"..............................................##ccccc##.........",
"...............................................#######..........",
"................................................#####...........",
"................................................................",
"................................................................",
"................................................................",
"................................................................"};
"""
        return icon

class SubObjectLoft:
    def __init__(self, obj, sub1=None, sub2=None):
        obj.Proxy = self
        obj.addProperty("App::PropertyLinkSub", "Link1", "SubObjectLoft").Link1 = sub1
        obj.addProperty("App::PropertyLinkSub", "Link2", "SubObjectLoft").Link2 = sub2
        obj.addProperty("App::PropertyFloat", "VertexFaceScale", "SubObjectLoft","""
When lofting a face to a vertex a tiny copy of the face is created to take the place
of the vertex in the loft.  This property sets the scale of that face.""" \
).VertexFaceScale = 1e-5
        obj.addProperty("App::PropertyRotation","VertexFaceRotation","SubObjectLoft","""
Rotation offset of the scaled copy of the face made when lofting to a vertex.  This
property is relative to the face's local coordinate system, z=1 always being the normal.""")

    
    def execute(self, fp):
        sub1 = getattr(fp.Link1[0].Shape, fp.Link1[1][0]) if fp.Link1 else None
        sub2 = getattr(fp.Link2[0].Shape, fp.Link2[1][0]) if fp.Link2 else None
        if not bool(sub1 and sub2):
            fp.Shape = Part.Shape()
            return
        #subs have already been filtered to be either 2 faces or 1 face and 1 vertex
        if hasattr(sub1,"Surface") and hasattr(sub2,"Surface"): #2 faces
            fp.setEditorMode("VertexFaceRotation", 1) #hidden
            fp.setEditorMode("VertexFaceScale", 1)
            inners = []
            outer = None
            for wire1,wire2 in zip(sub1.Wires, sub2.Wires):
                if wire1.isSame(sub1.OuterWire):
                    outer = Part.makeLoft([wire1, wire2], True)
                else:
                    inners.append(Part.makeLoft([wire1, wire2], True))
            if outer:
                for inner in inners:
                    outer = outer.cut(inner)
                fp.Shape = outer
        elif hasattr(sub1, "Point") or hasattr(sub2, "Point"):
            fp.setEditorMode("VertexFaceRotation", 0) #shown
            fp.setEditorMode("VertexFaceScale", 0)
            if hasattr(sub1, "Point"):
                face = sub2
                v = sub1
            else:
                face = sub1
                v = sub2
                
            copy = face.copy()
            copy.scale(fp.VertexFaceScale, copy.CenterOfGravity)
            dir = (face.CenterOfGravity - v.CenterOfGravity).normalize()
            copy.Placement.Rotation = FreeCAD.Rotation(copy.normalAt(0,0), dir).multiply(fp.VertexFaceRotation)
            copy.translate(v.CenterOfGravity - copy.CenterOfGravity)

            inners = []
            outer = None
            for wire1, wire2 in zip(face.Wires, copy.Wires):
                if wire1.isSame(face.OuterWire):
                    outer = Part.makeLoft([wire1, wire2], True)
                else:
                    inners.append(Part.makeLoft([wire1, wire2], True))
            for inner in inners:
                outer = outer.cut(inner)
            fp.Shape = outer
        else:
            FreeCAD.Console.PrintError("SubObjectLoft error: Links must be either 2 faces or 1 face and 1 vertex\n")
            fp.Shape = Part.Shape()
        
    def show(self, shape, label):
        """debug helper, can show a shape without throwing recursive recompute error"""
        obj = FreeCAD.ActiveDocument.addObject("Part::Feature",label)
        obj.Shape = shape

#end SubObjectLoftVP class
        
        
class SubObjectLoftVP:
    def __init__(self, obj):
        obj.Proxy = self
        
    def attach(self, obj):
        self.Object = obj
        
    def __getstate__(self):
        return None
    
    def __setstate__(self,state):
        return None
    
    def getDisplayModes(self,obj):
        '''Return a list of display modes.'''
        modes=[]
        modes.append("Flat Lines")
        return modes

    def getIcon(self):
        cmd = MeshRemodelSubObjectLoftCommandClass()
        return cmd.GetResources()["Pixmap"]
#end SubObjectLoftVP class

class MeshRemodelSubObjectLoftCommandClass(object):
    def __init__(self):
        self.subs = [] #list of tuples [(obj1,(subnames)),(obj2,(subnames))]
        self.face_count = 0 #counts of selected faces, edges, vertices
        self.edge_count = 0
        self.vertex_count = 0
    
    def GetResources(self):
        return {'Pixmap': os.path.join(iconPath, "SubObjectLoft.svg"),
                "MenuText": "SubObject Loft",
                "ToolTip":\
"""Loft between selected subobjects.

Selection (result)
2 faces (SubObjectLoft -> Solid)
1 face, 1 vertex (SubObjectLoft -> Solid)
2 vertices (ParametricLine -> Line)
2 edges (Part::Loft -> Face)
1 edge, 1 vertex (Part::Loft -> Face)

Not supported:
1 edge, 1 face

Alt+Click = reverses references to reverse object for all but SubObjectLofts
"""
        }
        
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        doc.openTransaction("MeshRemodel:SubObjectLoft")

        if self.face_count >= 1:
            loft = doc.addObject("Part::FeaturePython","SubObjectLoft")
            SubObjectLoft(loft, self.subs[0], self.subs[1])
            SubObjectLoftVP(loft.ViewObject)

        elif self.vertex_count == 2:
            loft = doc.addObject("Part::FeaturePython","ParametricLine")
            if not QtCore.Qt.AltModifier & QtGui.QApplication.keyboardModifiers():
                ParametricLine(loft, self.subs[0], self.subs[1])
            else:
                ParametricLine(loft, self.subs[1], self.subs[0])
            ParametricLineVP(loft.ViewObject)
        else:
            loft = doc.addObject("Part::Loft","SubObjectLoft")
            #Part::Loft cannot use subobjects, so we use subshapebinders
            binder1 = doc.addObject("PartDesign::SubShapeBinder","Binder")
            binder1.Support = [self.subs[0]]
            binder2 = doc.addObject("PartDesign::SubShapeBinder","Binder")
            binder2.Support = [self.subs[1]]
            if not QtCore.Qt.AltModifier & QtGui.QApplication.keyboardModifiers():
                loft.Sections = [binder1, binder2]
            else:
                loft.Sections = [binder2, binder1]
            binder1.ViewObject.Visibility = False
            binder2.ViewObject.Visibility = False
            loft.Solid = False
        doc.commitTransaction()
        doc.recompute()
        
    def IsActive(self):
        self.subs = []
        self.face_count = 0
        self.edge_count = 0
        self.vertex_count = 0
        if not FreeCAD.ActiveDocument:
            return False
        sel = FreeCADGui.Selection.getCompleteSelection()
        for s in sel:
            if not s.Object.isDerivedFrom("Part::Feature"):
                return False
            self.subs.append((s.Object, s.SubElementNames))
            if "Face" in s.SubElementNames[0]:
                self.face_count += 1
            elif "Edge" in s.SubElementNames[0]:
                self.edge_count += 1
            elif "Vertex" in s.SubElementNames[0]:
                self.vertex_count += 1
        if len(self.subs) == 2:
            if self.face_count == 1 and self.edge_count == 1:
                return False #cannot loft between face and edge
            return True
        return False
        

# end subobject loft command class
###################################################################


class GridSurface:
    def __init__(self, obj):
        obj.Proxy = self

        grp = "Dimensions"
        obj.addProperty("App::PropertyLength","YInterval",grp,"Length in Y direction, triggers grid rebuild").YInterval=10
        obj.addProperty("App::PropertyLength","XInterval",grp,"Length in X direction, triggers grid rebuild").XInterval=10

        obj.addProperty("App::PropertyIntegerConstraint","CountRows",grp,"Whenu using Custom grid size this works.").CountRows = (3,1,10000,1)
        obj.addProperty("App::PropertyIntegerConstraint","CountColumns",grp,"When using Custom grid size this works.").CountColumns = (3,1,10000,1)

        grp = "Vertices"
        obj.addProperty("App::PropertyLinkList","Vertices",grp,"Part::Vertex objects making up grid")
        obj.addProperty("App::PropertyInteger","PointSize",grp,"Point size applied to Part::Vertexes").PointSize = 8
        obj.addProperty("App::PropertyBool","GridVisibility",grp,"Visibility toggler for Part::Vertexes").GridVisibility = True

        grp = "GridSurface"
        obj.addProperty("App::PropertyBool","ColumnMode",grp,"If True, the grid is treated as column-oriented instead of row-oriented").ColumnMode = False
        obj.addProperty("App::PropertyBool","Reversed",grp,"Flips the normal of the face").Reversed = False
        obj.addProperty("App::PropertyEnumeration","Method",grp,"Method for making the bspline surface, LoftedRows makes bsplines out each row and lofts them").Method = ["Interpolate Points","Lofted Edges","Lofted Edges Ruled"]
        obj.Method = "Lofted Edges"
        obj.addProperty("App::PropertyEnumeration","RowWireType",grp,"Only for loft methods, the type of wire to make from each row of points").RowWireType = \
        ["Open BSpline", "Closed BSpline", "Open Polyline", "Closed Polyline"]
        obj.addProperty("App::PropertyEnumeration","Output",grp, "Surface, Edges, Gordon Template").Output = ["Surface","Edges", "Gordon Template"]
        obj.Output="Surface"
        obj.addProperty("App::PropertyLinkSubList", "PreWireEdges",grp,"When using Loft Row method the wire formed from these edges is lofted to row 0 as part of the loft.")
        obj.addProperty("App::PropertyLinkSubList","PostWireEdges",grp,"When using Loft Row method the wire formed form these edges is lofted to the last row as part of the loft.")
        obj.addProperty("App::PropertyString","Version",grp,"Version of GridSurface used to create this object.").Version = __version__

        self.done = True
        self.blockSignals = False


    def onChanged(self, fp, prop):
        if not hasattr(self, "done"):
            return
        rebuilders = ["YInterval","XInterval"]
        if prop in rebuilders:
            fp.ViewObject.Proxy.onDelete(None, None) #delete all children to trigger rebuild of grid
        elif prop == "CountRows":
            num_rows = self.countRows(fp)
            while fp.CountRows > num_rows:
                self.addRow(fp)
                num_rows = self.countRows(fp)
            while fp.CountRows < num_rows:
                self.removeRow(fp)
                num_rows = self.countRows(fp)
        elif prop == "CountColumns":
            num_cols = self.countColumns(fp)
            while fp.CountColumns > num_cols:
                self.addColumn(fp)
                num_cols = self.countColumns(fp)
            while fp.CountColumns < num_cols:
                self.removeColumn(fp)
                num_cols = self.countColumns(fp)


    def makeSplines(self, fp, vrows):
        """if fp.Output is "Edges" or "Gordon Template", then the shape becomes the rows of edges instead of a surface
        """

        splines = []
        left = [] #left and right edges of gordon template
        right = []
        for row in vrows:
            spline = Part.BSplineCurve()
            if "BSpline" in fp.RowWireType:
                spline.interpolate(row, PeriodicFlag = "Closed BSpline" == fp.RowWireType)
            elif "Open Polyline" == fp.RowWireType:
                spline = Part.makePolygon(row)
            elif "Closed Polyline" in fp.RowWireType:
                spline = Part.makePolygon(row + [row[0]])
            splines.append(spline.toShape() if hasattr(spline, "toShape") else spline)
            if fp.Output == "Gordon Template":
                left.append(row[0])
                right.append(row[-1])
        if fp.Output == "Gordon Template":
            left_spline = Part.BSplineCurve()
            right_spline = Part.BSplineCurve()
            if "BSpline" in fp.RowWireType:
                left_spline.interpolate(left)
                right_spline.interpolate(right)
                splines = [left_spline.toShape(), right_spline.toShape()] + splines
            elif "Polyline" in fp.RowWireType:
                left_spline = Part.makePolygon(left)
                right_spline = Part.makePolygon(right)
                splines = [left_spline, right_spline] + splines



        fp.Shape = Part.makeCompound(splines)
        return

    def makeWire(self, edges):
        """make the post and pre wire edges"""
        if not edges:
            return None
        subobjects = []
        for edge in edges:
            feat = edge[0]
            subs = edge[1]
            for sub in subs:
                if "Edge" in sub:
                    subobjects.append(feat.getSubObject(sub))
                else:
                    FreeCAD.Console.PrintWarning(f"GridSurface: {sub} is not supported subobject type, skipping...\n")
        if not subobjects:
            return None
        return Part.Wire(subobjects)

        return None



    def execute(self,fp):
        if self.blockSignals:
            self.blockSignals = False
            return

        vertices = fp.Vertices
        if not vertices:
            self.buildGrid(fp)
            vertices = fp.Vertices

        vert_dict = self.buildDictionary(fp)
        rows = fp.CountRows
        cols = fp.CountColumns
        if fp.ColumnMode:
            rows,cols = (cols,rows)
        vrows = []
        for rr in range(rows):
            this_row = []
            for cc in range(cols):
                if not fp.ColumnMode:
                    vert = self.fetch(fp,vert_dict,rr,cc)
                else:
                    vert = self.fetch(fp,vert_dict,cc,rr)
                if not vert: #allow user to delete individual vertices
                    continue
                try:
                    this_row.append(vert.Shape.CenterOfGravity)
                except:
                    this_row.append(vert.Placement.Base)
            if this_row:
                vrows.append(this_row)
        if fp.Output == "Edges" or fp.Output == "Gordon Template":
            self.makeSplines(fp, vrows)
            return

        surf = Part.BSplineSurface()
        if fp.Method == "Interpolate Points":
            if len(vrows) > 1:
                surf.interpolate(vrows)
            else:
                spline = Part.BSplineCurve()
                spline.interpolate(vrows[0])
                fp.Shape = spline.toShape()
                return
        elif "Lofted Edges" in fp.Method:

            preSpline = self.makeWire(fp.PreWireEdges) #returns None if no edges are defined
            postSpline = self.makeWire(fp.PostWireEdges)
            splines = []
            ruled = True if "Ruled" in fp.Method else False
            for row in vrows:
                if "BSpline" in fp.RowWireType:
                    spline = Part.BSplineCurve()
                    periodicFlag = "Closed" in fp.RowWireType
                    spline.interpolate(row, PeriodicFlag = periodicFlag)
                    splines.append(spline)
                else: # polyines
                    if "Closed" in fp.RowWireType:
                        row.append(row[0])
                    poly = Part.makePolygon(row)
                    splines.append(poly)
            if preSpline:
                splines = [preSpline] + splines
            if postSpline:
                splines = splines + [postSpline]
            if len(splines) > 1:
                loft = Part.makeLoft(splines, ruled = ruled)
                surf = loft
            else: #just make a bspline edge since there is only 1 row
                fp.Shape = splines[0].toShape() if splines else Part.Shape()
                return
        fp.Shape = surf.toShape() if hasattr(surf, "toShape") else surf
        if fp.Reversed:
            fp.Shape = fp.Shape.copy().reversed()

    def buildDictionary(self, fp):
        """build and return a dictionary object of all the fp.Vertices objects"""
        vert_dict = {}
        verts = [v for v in fp.Vertices if v is not None]
        for vert in verts:
            vert.ViewObject.PointSize = fp.PointSize
            vert.ViewObject.Visibility = fp.GridVisibility
            row,col = (vert.Row, vert.Column)
            if row == -1 and col == -1:
                continue
            vert_dict[vert.Name] = (vert,row,col)
        return vert_dict


    def fetch(self,fp,dictionary,row,col):
        """fetch the Part::Vertex from the dictionary"""
        if len(dictionary.items()) == 0:
            return None
        for k,v in dictionary.items():
            if v[1] == row and v[2] == col:
                return v[0]
        return None


    def countRows(self,fp):
        """count the number of rows in the current grid"""
        verts_dict = self.buildDictionary(fp)
        high_row = 0
        for k,v in verts_dict.items():
            vert,row,col = v
            if row > high_row:
                high_row = row
        return high_row + 1

    def countColumns(self,fp):
        """count the number of columns in the longest row"""
        verts_dict = self.buildDictionary(fp)
        high_col = 0
        for k,v in verts_dict.items():
            vert, row, col = v
            if col > high_col:
                high_col = col
        return high_col + 1


    def buildGrid(self, fp):
        """constructs the grid of vertices, only called from execute() when fp.Vertices is empty"""
        if self.blockSignals:
            self.blockSignals = False
            return
        doc = fp.Document
        num_x = fp.CountColumns
        num_y = fp.CountRows
        used = []
        for xx in range(num_x):
            for yy in range(num_y):
                if not (xx,yy) in used:
                    used.append((xx,yy))
                    x = xx * fp.XInterval
                    y = yy * fp.YInterval
                    self.addVertex(fp, xx, yy, x, y)


    def addColumn(self,fp):
        """add a new column vertex to each row"""
        num_rows = self.countRows(fp)
        num_cols = self.countColumns(fp)
        x = (num_cols) * fp.XInterval.Value
        for rr in range(num_rows):
            y = rr * fp.YInterval
            vert = self.addVertex(fp, num_cols, rr, x, y)


    def removeColumn(self,fp):
        """remove the last column"""
        num_cols = self.countColumns(fp)
        to_remove = [v for v in fp.Vertices if v.Column == num_cols - 1]
        for t in to_remove:
           # self.blockSignals = True
            fp.Document.removeObject(t.Name)

    def addVertex(self, fp, col, row, x, y):

        vert = fp.Document.addObject("Part::Vertex",f"{fp.Name}_Vertex")
        vert.addProperty("App::PropertyInteger","Row",f"{fp.Label}",f"Row in the {fp.Label} object").Row = row
        vert.addProperty("App::PropertyInteger","Column",f"{fp.Label}",f"Column in the {fp.Label} object").Column = col
        vert.setEditorMode("Row",1)
        vert.setEditorMode("Column",1)
        vert.Placement.Base.x = x
        vert.Placement.Base.y = y
        vert.recompute()
        fp.Vertices += [vert]
        return vert

    def addRow(self, fp):
        """add a new row to end of the current rows"""
        num_rows = self.countRows(fp)
        y = (num_rows) * fp.YInterval.Value
        verts_added = []
        for cc in range(fp.CountColumns):
            x = cc * fp.XInterval
            vert = self.addVertex(fp, cc, num_rows, x, y)


    def removeRow(self, fp):
        """remove the last row of vertices"""
        num_rows = self.countRows(fp)
        to_remove = [v for v in fp.Vertices if v.Row == num_rows -1 ]
        for t in to_remove:
            fp.Document.removeObject(t.Name)



class GridSurfaceVP:
    def __init__(self, vobj):
        vobj.Proxy = self
        self.name = vobj.Object.Name
        self.outputMode = "Surface"

    @property
    def FP(self):
        return FreeCAD.ActiveDocument.getObject(self.name)

    def setupContextMenu(self, vobj, menu):
        def selectRowTrigger(row): return lambda: self.selectRow(row)
        def selectColTrigger(col): return lambda: self.selectCol(col)
        def outputMode(mode): return lambda: self.switchOutput(mode)
        if not shiboken.isValid(menu):
            FreeCAD.Console.PrintMessage(f"{self.FP.Name}: cannot create context menu.  You may try restarting FreeCAD.\n")
            return
        fp = vobj.Object
        num_rows = fp.CountRows
        num_cols = fp.CountColumns
        selectRowMenu = menu.addMenu("Select row")
        selectColMenu = menu.addMenu("Select column")

        for row in range(num_rows):
            action = selectRowMenu.addAction(f"Row {row}")
            action.triggered.connect(selectRowTrigger(row))

        for col in range(num_cols):
            action = selectColMenu.addAction(f"Column {col}")
            action.triggered.connect(selectColTrigger(col))

        if fp.PostWireEdges or fp.PreWireEdges:
            swapAction = menu.addAction("Swap Pre and Post Wire Edges")
            swapAction.triggered.connect(self.swapPrePostWireEdges)

        outputMenu = menu.addMenu("Output Mode")

        if fp.Output != "Surface":
            outputAction = outputMenu.addAction("Surface output")
            outputAction.triggered.connect(outputMode("Surface"))

        if fp.Output != "Edges":
            outputAction = outputMenu.addAction("Edges output")
            outputAction.triggered.connect(outputMode("Edges"))

        if fp.Output != "Gordon Template":
            outputAction = outputMenu.addAction("Gordon Template output")
            outputAction.triggered.connect(outputMode("Gordon Template"))

        if fp.Output == "Gordon Template":
            selectAction = menu.addAction("Select for Gordon Surface")
            selectAction.triggered.connect(self.selectForGordonTemplate)

        visibilityAction = menu.addAction("Toggle visibility of vertices")
        visibilityAction.triggered.connect(self.toggleVisibility)

        mode = "Row Mode" if fp.ColumnMode else "Column Mode"
        columnModeAction = menu.addAction(f"Switch to {mode}")
        columnModeAction.triggered.connect(self.toggleColumnMode)

        menu.addSeparator()

    def toggleColumnMode(self, checked):
        fp = self.FP
        fp.ColumnMode = not fp.ColumnMode
        fp.recompute()

    def toggleVisibility(self, checked):
        """Toggle visibility of Part::Vertex objects"""
        fp = self.FP
        for vert in fp.Vertices:
            vert.Visibility = not vert.Visibility

    def selectForGordonTemplate(self, checked):
        """select the edges in proper order to make a Gordon surface object in Curves workbench"""
        fp = self.FP
        sel = FreeCADGui.Selection
        sel.clearSelection()
        for idx,edge in enumerate(fp.Shape.Edges):
            sel.addSelection(fp.Document.Name, fp.Name, f"Edge{idx + 1}")

    def switchOutput(self, mode):
        """switch the output type to Surface, Edges, or Gordon Template"""
        fp = self.FP
        fp.Output = mode
        fp.recompute()

    def swapPrePostWireEdges(self, checked):
        """swaps the pre and post wire edge properties"""
        fp = self.FP
        fp.PreWireEdges,fp.PostWireEdges = (fp.PostWireEdges,fp.PreWireEdges)
        fp.recompute()

    def selectRow(self, row):
        fp = self.FP
        FreeCADGui.Selection.clearSelection()
        for vert in fp.Vertices:
            if vert.Row == row:
                FreeCADGui.Selection.addSelection(fp.Document.Name, vert.Name)

    def selectCol(self, col):
        fp = self.FP
        FreeCADGui.Selection.clearSelection()
        for vert in fp.Vertices:
            if vert.Column == col:
                FreeCADGui.Selection.addSelection(fp.Document.Name, vert.Name)

    def claimChildren(self):
        fp = self.FP
        return fp.Vertices

    def getIcon(self):
        cmd = MeshRemodelCreateGridSurfaceCommandClass()
        return cmd.GetResources()["Pixmap"]

    def onDelete(self, vobj, obj):
        fp = self.FP
        verts = fp.Vertices
        for vert in verts:
            if vert:
                fp.Document.removeObject(vert.Name)
        return True

    def __getstate__(self):
        '''When saving the document this object gets stored using Python's json module.\
            If we have any unserializable stuff return them here or None'''
        return None

    def __setstate__(self,state):
        '''When restoring the serialized object from document we have the chance to set some internals here.\
                Since no data were serialized nothing needs to be done here.'''
        return None

class MeshRemodelCreateGridSurfaceCommandClass(object):
    def __init__(self):
        pass
    
    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateGridSurface.svg') ,
            'MenuText': "Create GridSurfacee" ,
            'ToolTip' : fixTip("""
Create a GridSurface object.

This object consists of a grid of Part::Vertex objects, each of which
can be manipulated individually.  From these vertices the GridSurface
object can produce a surface, edges from the rows or columns, or a Gordon 
Surface template object.  The Gordon Surface Template can be the 
source object for a Gordan Surface object created in Curves workbench.

In Surface mode you can try different methods, including Lofting the edges
and interpolating the points directly into a surface.  The loft can be ruled
or (default) unruled.  Edges can be open or closed bsplines or polylines.

In Loft mode you can also select edges of other objects to form part of the 
surface by adding them before (Pre Wire Edges) and after (Post Wire Edges)
the edges in the GridSurface object.

""")}
    
    def Activated(self):

        doc = FreeCAD.ActiveDocument if FreeCAD.ActiveDocument else FreeCAD.newDocument()
        fp = doc.addObject("Part::FeaturePython", "GridSurface")
        GridSurface(fp)
        GridSurfaceVP(fp.ViewObject)
        doc.recompute()

    def isActive(self):
        return True

# end grid surface command class
###################################################################
# Create a wire from selected subobjects
class MeshRemodelCreateWireCommandClass(object):
    """Create a wire from selected objects or subobjects"""
    
    def __init__(self):
        self.shapes = [] #can be full shapes or subobjects
        
    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateWire.svg') ,
            'MenuText': "Create &Wire" ,
            'ToolTip' : fixTip("\
Create a flattened wire from selected subobjects.  Nonplanar objects\n\
are flattened to the plane defined by the first arc or circle found, \n\
otherwise by some key points: first point of first subobject, last point \n\
of last object, and the middle point of the key points, which are the starting\n\
and ending points of the various edges.\n\
\n\
Needs at least 2 edges.  Wire should be closed or nearly closed.  Works on \n\
SubShapeBinders, too, by selecting entire object.  Subobjects need not be \n\
from a single object.\n\
\n\
Alt+Click = Use camera orientation for plane normal (circles become bsplines)\n")}      
            
    def Activated(self):
        def rv(vector):
            """used for printing display only"""
            return FreeCAD.Vector(round(vector.x,6), round(vector.y,6), round(vector.z,6))
        useCameraOrientation = False
        modifiers = QtGui.QApplication.keyboardModifiers()
        if modifiers & QtCore.Qt.AltModifier:
            useCameraOrientation = True
        plm = FreeCAD.Placement()
        if not useCameraOrientation:
            #look for a circle/arc and use that normal if found
            base,normal = None,None
            for shape in self.shapes:
                if hasattr(shape,"Curve") and shape.Curve.TypeId == "Part::GeomCircle":
                    base = shape.Curve.Location
                    normal = shape.Curve.Axis
                    break
            if not base:
                keypts = self.getKeyPoints()
                if keypts[0] == keypts[-1]:
                    keypts.pop()
                base,normal = gu.planeFromPoints(keypts)
        else:
            #get from camera orientation
            rot = FreeCADGui.ActiveDocument.ActiveView.getCameraOrientation()
            base = self.shapes[0].Vertex1.Point
            normal = rot.multVec(FreeCAD.Vector(0,0,1))

        new_shapes = []
        firstStartPoint = None
        previousEndPoint = None
        wasReversed = False
        self.shapes.append(self.shapes[0])
        for idx,(shape,next) in enumerate(zip(self.shapes[:-1],self.shapes[1:])):
            if not hasattr(shape,"Curve"):
                FreeCAD.Console.PrintWarning("MeshRemodel Create Flat Wire: Faces not supported at this time.\n")
                continue
            firstEdge = bool(idx == 0)
            lastEdge = bool(idx == len(self.shapes)-2)
            if shape.Curve.TypeId == "Part::GeomLine":
                startPoint = shape.valueAt(shape.ParameterRange[0])
                endPoint = shape.valueAt(shape.ParameterRange[1])
                nextStart = next.valueAt(next.ParameterRange[0])
                nextEnd = next.valueAt(next.ParameterRange[1])
                startPoint = gu.projectVectorToPlane(startPoint, base, normal)
                endPoint = gu.projectVectorToPlane(endPoint, base, normal)
                nextStart = gu.projectVectorToPlane(nextStart, base, normal)
                nextEnd = gu.projectVectorToPlane(nextEnd, base, normal)
                if wasReversed:
                    startPoint,endPoint = endPoint,startPoint
                #print(f"startPoint:{rv(startPoint)}, endPoint:{rv(endPoint)}, nextStart:{rv(nextStart)}, nextEnd:{rv(nextEnd)}")
                needsReversing = self.needsReversing(startPoint, endPoint, nextStart, nextEnd)
                #print(f"needsReversing: {needsReversing}, idx: {idx}")
                if needsReversing:
                    startPoint,endPoint = endPoint,startPoint
                    wasReversed = True
                else:
                    wasReversed = False
                previousEndPoint = endPoint
                if firstEdge:
                    firstStartPoint = startPoint
                elif lastEdge:
                    endPoint = firstStartPoint
                #print(f"making line segment from {rv(startPoint)} to {rv(endPoint)}")
                new_shape = Part.LineSegment(startPoint, endPoint).toShape()
                new_shapes.append(new_shape)
                #print(f"appending shape")
            elif shape.Curve.TypeId == "Part::GeomBSplineCurve" or shape.Curve.TypeId == "Part::GeomEllipse" \
                                        or bool(useCameraOrientation and shape.Curve.TypeId == "Part::GeomCircle"):
                numpoints = int(shape.Curve.length()*4)
                numpoints = 10 if numpoints < 10 else numpoints
                pts = shape.discretize(numpoints)
                pts2 = []
                for p in pts:
                    pts2.append(gu.projectVectorToPlane(p, base, normal))
                startPoint = shape.valueAt(shape.ParameterRange[0])
                endPoint = shape.valueAt(shape.ParameterRange[1])
                nextStart = next.valueAt(next.ParameterRange[0])
                nextEnd = next.valueAt(next.ParameterRange[1])
                startPoint = gu.projectVectorToPlane(startPoint, base, normal)
                endPoint = gu.projectVectorToPlane(endPoint, base, normal)
                nextStart = gu.projectVectorToPlane(nextStart, base, normal)
                nextEnd = gu.projectVectorToPlane(nextEnd, base, normal)
                if wasReversed:
                    startPoint,endPoint = endPoint,startPoint
                #print(f"startPoint:{rv(startPoint)}, endPoint:{rv(endPoint)}, nextStart:{rv(nextStart)}, nextEnd:{rv(nextEnd)}")
                needsReversing = self.needsReversing(startPoint, endPoint, nextStart, nextEnd)
                #print(f"needsReversing: {needsReversing}, idx: {idx}")
                if needsReversing:
                    startPoint,endPoint = endPoint,startPoint
                    wasReversed = True
                else:
                    wasReversed = False
                if firstEdge:
                    firstStartPoint = startPoint
                elif lastEdge:
                    endPoint = firstStartPoint
                previousEndPoint = endPoint
                #print(f"making bspline from {rv(startPoint)} to {rv(endPoint)}")
                new_shape = Part.BSplineCurve()
                if wasReversed:
                   pts2.reverse()
                new_shape.interpolate(pts2)
                new_shapes.append(new_shape.toShape())

            elif shape.Curve.TypeId == "Part::GeomCircle":
                #handle arcs and circles
                startPoint = shape.valueAt(shape.ParameterRange[0])
                endPoint = shape.valueAt(shape.ParameterRange[1])
                nextStart = next.valueAt(next.ParameterRange[0])
                nextEnd = next.valueAt(next.ParameterRange[1])
                startPoint = gu.projectVectorToPlane(startPoint, base, normal)
                endPoint = gu.projectVectorToPlane(endPoint, base, normal)
                nextStart = gu.projectVectorToPlane(nextStart, base, normal)
                nextEnd = gu.projectVectorToPlane(nextEnd, base, normal)
                if wasReversed:
                    startPoint,endPoint = endPoint,startPoint
                #print(f"startPoint:{rv(startPoint)}, endPoint:{rv(endPoint)}, nextStart:{rv(nextStart)}, nextEnd:{rv(nextEnd)}")
                needsReversing = self.needsReversing(startPoint, endPoint, nextStart, nextEnd)
                #print(f"needsReversing: {needsReversing}, idx: {idx}")
                if needsReversing:
                    startPoint,endPoint = endPoint,startPoint
                    wasReversed = True
                else:
                    wasReversed = False
                if firstEdge:
                    firstStartPoint = startPoint
                elif lastEdge:
                    endPoint = firstStartPoint if endPoint.distanceToPoint(firstStartPoint) < 1e-3 else endPoint
                previousEndPoint = endPoint
                #print(f"making arc/circle from {rv(startPoint)} to {rv(endPoint)}")
                midPoint = shape.discretize(3)[1]
                midPoint = gu.projectVectorToPlane(midPoint, base, normal)
                if startPoint.distanceToPoint(endPoint) < epsilon:
                    center = gu.projectVectorToPlane(shape.Curve.Center, base, normal)
                    new_shape = Part.Circle(center, normal, shape.Curve.Radius)
                else:
                    new_shape = Part.ArcOfCircle(startPoint, midPoint, endPoint)
                new_shapes.append(new_shape.toShape())


        #print(f"new_shapes: {new_shapes}, len(new_shapes): {len(new_shapes)}")
        try:
            new_shapes = Part.sortEdges(new_shapes)
        except Exception as e:
            FreeCAD.Console.PrintWarning("MeshRemodel: Part.sortEdges function failed. {e}\n")
        #print(f"new_shapes: {new_shapes}, len(new_shapes): {len(new_shapes)}")
        wires = []
        for ns in new_shapes:
            wirename = "MR_Wire"
            wire = Part.Wire(ns)
            if wire.isClosed():
                wire = Part.makeFace([wire],"Part::FaceMakerBullseye")
                wirename = "MR_Face"
            else: #wire not closed, so close it
                p1 = wire.Vertexes[0].Point
                p2 = wire.Vertexes[-1].Point
                if p2 == p1:
                    wire = Part.Face([wire], "Part::FaceMakerBullseye")
                    wirename = "MR_Face"
                elif p2.distanceToPoint(p1) < 1e-3:
                    edges = wire.Edges
                    edges.append(Part.LineSegment(p1,p2).toShape())
                    wire = Part.Wire(Part.__sortEdges__(edges))
                    wire = Part.makeFace([wire], "Part::FaceMakerBullseye")
                    wirename = "MR_Face"
            wires.append(Part.show(wire, wirename))
        if len(wires) > 1:
            new_wires = [w.Shape.Wire1 for w in wires]
            face = Part.makeFace(new_wires, "Part::FaceMakerBullseye")
            Part.show(face,"MR_Face")
        elif len(wires) == 1:
            sel = Gui.Selection.getSelection()
            for s in sel:
                s.ViewObject.Visibility = False
            
        
    def needsReversing(self, startPoint, endPoint, nextStart, nextEnd):
        """determine if this edge needs to be reversed to be coincident with the previous shape"""
        dist1 = startPoint.distanceToPoint(nextStart)
        dist2 = endPoint.distanceToPoint(nextStart)
        dist3 = startPoint.distanceToPoint(nextEnd)
        dist4 = endPoint.distanceToPoint(nextEnd)
        lis = [dist1,dist2,dist3,dist4]

        shortest = lis.index(min(lis))
        bools = [True,False,True,False]
        #print(f"distances:{[dist1,dist2,dist3,dist4]}, shortest:{shortest}, bools[shortest]:{bools[shortest]}")
        return bools[shortest]
        
    def getKeyPoints(self):
        pts = []
        for shape in self.shapes:
            if hasattr(shape,"ParameterRange"):
                if len(shape.ParameterRange) == 2:
                    pts.append(shape.valueAt(shape.ParameterRange[0]))
                    pts.append(shape.valueAt(shape.ParameterRange[1]))
                elif len(shape.ParameterRange) >= 4:
                    pts.append(shape.valueAt(shape.ParameterRange[0],shape.ParameterRange[1]))
                    pts.append(shape.valueAt(shape.ParameterRange[2],shape.ParameterRange[3]))
            else:
                raise Exception (f"getKeyPoints(): Unsupported shape type: {shape.Curve.TypeId}")
        while pts[-1].distanceToPoint(pts[0]) < epsilon:
            pts.pop()
        return pts
        
    def IsActive(self):
        self.shapes = []
        sel = Gui.Selection.getCompleteSelection()
        if not sel:
            return False
        for s in sel:
            if not self.subElements(s):
                return False
            element = self.subElements(s)
            if hasattr(element, "X") and hasattr(element, "Y") and hasattr(element, "Z"): #vertex
                return False
            if hasattr(element, "Surface"):
                return False
            if isinstance(element,list) or isinstance(element,tuple):
                for el in element:
                    self.shapes.append(el)
            else:
                self.shapes.append(element)
        if self.shapes:
            try:
                self.shapes = Part.sortEdges(self.shapes)
            except:
                pass
            self.shapes = self.flattenList(self.shapes)
            if len(self.shapes) < 2:
                return False
        else:
            return False
        return True
    
    def subElements(self, selobj):
        if selobj.SubElementNames == ('',):
            if hasattr(selobj.Object,"Shape"):
                return selobj.Object.Shape.Edges
            return []
        return selobj.SubObjects
        
    def flattenList(self, nested_list):
        flattened_list = []
        for item in nested_list:
            if isinstance(item, list):
                flattened_list.extend(self.flattenList(item))
            else:
                flattened_list.append(item)
        return flattened_list

###################################################################
# Create a plane

class MeshRemodelCreatePlaneCommandClass(object):
    """Create a plane"""
    def __init__(self):
        self.subs = []
        
    def GetResources(self):
        return {
            'Pixmap'  : os.path.join( iconPath , 'CreatePlane.svg') ,
            'MenuText': "Create plane" ,
            'ToolTip' : fixTip(\
"""Create a 50 x 50 Part::Plane object.  Position it on the plane defined
by selected subobjects:

3 points = plane defined by the 3 points
1 edge (line) = plane's normal is the edge's direction
1 edge (circle or arc) = plane's normal is circle's normal
2 edges = plane defined by 3 points from the edges.
1 face = face normal
no selection = xy plane at origin

""")}

    def Activated(self):
        doc = FreeCAD.ActiveDocument
        self.subs = []
        sel = FreeCADGui.Selection.getCompleteSelection()
        count_vertices = 0
        count_faces = 0
        count_edges = 0
        for s in sel:
            if s.HasSubObjects and s.SubObjects:
                self.subs.append(s.SubObjects[0])
                for name in s.SubElementNames:
                    if "Vertex" in name:
                        count_vertices += 1
                    elif "Face" in name:
                        count_faces += 1
                    elif "Edge" in name:
                        count_edges += 1
        # print(f"self.subs: {self.subs}, vertices, edges, faces = {count_vertices, count_edges, count_faces}")
        base,normal = FreeCAD.Vector(0,0,0), FreeCAD.Vector(0,0,1)
        if count_vertices == 3 and count_faces + count_edges == 0:
            pts = [obj.Point for obj in self.subs if hasattr(obj,"Point")]
            base,normal = gu.planeFromPoints(pts)
        if count_faces == 1 and count_edges + count_vertices == 0:
            base,normal = self.subs[0].CenterOfGravity, self.subs[0].normalAt(0,0)
        if count_edges == 1 and count_faces + count_vertices == 0:
            edge = self.subs[0]
            if edge.Curve.TypeId == "Part::Circle" or edge.Curve.TypeId == "Part::ArcOfCircle":
                base,normal = edge.Curve.Center, edge.Curve.Axis
            elif edge.Curve.TypeId == "Part::GeomLine":
                base,normal = edge.CenterOfGravity, edge.Curve.Direction
        elif count_edges == 2 and count_faces + count_vertices == 0:
            edge1 = self.subs[0]
            edge2 = self.subs[1]
            pts = [edge1.discretize(10)[1], edge1.CenterOfMass, edge2.CenterOfMass]
            base,normal = gu.planeFromPoints(pts)
        elif global_picked and len(global_picked) >= 3:
            pts = global_picked
            base,normal = gu.planeFromPoints(pts)
        doc.openTransaction("Create plane")
        plane = doc.addObject("Part::Plane", "Plane")
        plane.Length = 50
        plane.Width = 50
        doc.recompute()
        plane.Placement.Rotation = FreeCAD.Rotation(FreeCAD.Vector(0,0,1), normal)
        #plane.Placement.translate(base)
        plane.Placement.translate(base - FreeCAD.Vector(plane.Shape.CenterOfGravity))
        doc.commitTransaction() 
        doc.recompute()
            


    def IsActive(self):
        if FreeCAD.ActiveDocument:
            return True
        return False


###################################################################

# Create a Circle from first 3 selected points

class MeshRemodelCreateCircleCommandClass(object):
    """Create Circle from first 3 selected points"""

    def __init__(self):
        self.pts = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateCircle.svg') ,
            'MenuText': "Create &circle" ,
            'ToolTip' : fixTip("Create a circle from first 3 selected points\n\
(Ctrl+Click to include Center point)\n\
(Ctrl+Shift+Click for only center)")}
 
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
        return {'Pixmap'  : os.path.join( iconPath , 'CreateArc.svg') ,
            'MenuText': "Create &arc" ,
            'ToolTip' : fixTip("Create an arc from first 3 selected points\n\
(Ctrl+Click to include Center point)\n\
(Ctrl+Shift+Click for only center)\n\
(Alt+Click for all permutations possible with the 3 selected points -- 6 arcs -- useful where you\n\
are not getting the arc orientation you were expecting -- will need to delete unwanted extra arcs)")}

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
        sel = Gui.Selection.getCompleteSelection()
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
        self.subs = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateSketch.svg') ,
            'MenuText': "Create s&ketch" ,
            'ToolTip' : fixTip("\
Create a new sketch.\n\
If no object is selected the attachment editor opens\n\
Ctrl+Click to make sketch out of selected objects\n\
Alt+Click make merged sketch out of selected objects\n\
Shift+Click 1st 3 points selected define sketch plane, points added as links to external geometry \n\
")}
 
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        if not self.objs:
            doc.openTransaction("New sketch")
            sk = doc.addObject("Sketcher::SketchObject", "Sketch")
            doc.commitTransaction()
            FreeCADGui.Selection.addSelection(sk)
            from AttachmentEditor import Commands
            Commands.editAttachment(sk)
            doc.recompute()
            return
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
            sketch.Label = 'Picked_Pts_Sketch'
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
        self.objs = Gui.Selection.getSelection()
        return True

# end create sketch class
####################################################################################

# Make a wire from selected objects
class MeshRemodelDraftUpgradeCommandClass(object):
    """Create wire from selected objects"""

    def __init__(self):
        self.objs = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'DraftUpgrade.svg') ,
            'MenuText': "Draft upgrade or downgrade" ,
            'ToolTip' : fixTip("Run Draft upgrade on selected objects\n\
Ctrl+Click to run Draft downgrade\n\
Tip: this can be used to upgrade connected edges into a wire, a wire \n\
into a face, or downgrade objects in the other direction.")}
 
    def Activated(self):
        doc = FreeCAD.ActiveDocument
        modifiers = QtGui.QApplication.keyboardModifiers()
        if (modifiers == QtCore.Qt.ControlModifier):
            doc.openTransaction("Draft downgrade")
            Draft.downgrade(self.objs)
            doc.recompute()
        else:
            selbackup = FreeCAD.Gui.Selection.getSelection()
            doc.openTransaction("Draft upgrade")
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

# end draft upgrade class

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
        return {'Pixmap'  : os.path.join( iconPath , 'AddSelectionObserver.svg') ,
            'MenuText': "Add Selection &Observer" ,
            'ToolTip' : "Allows to select vertices by Ctrl + preselection \
This is intended for use with BSpline selection of points, but \
can be used with all MeshRemodel tools that use selected points.  DOES NOT work if you \
mix points selected in this mode with normal selection mode in the same operation"}

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
        return {'Pixmap'  : os.path.join( iconPath , 'MergeSketches.svg') ,
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
        return {'Pixmap'  : os.path.join( iconPath , 'ValidateSketch.svg') ,
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

# check geometry
class MeshRemodelPartCheckGeometryCommandClass(object):
    """Check Geometry"""

    def __init__(self):
        self.objs = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CheckGeometry.svg') ,
            'MenuText': "Part CheckGeometry" ,
            'ToolTip' : "Execute Part workbench CheckGeometry tool"}
 
    def Activated(self):
        if not "Part_CheckGeometry" in Gui.listCommands():
            Gui.activateWorkbench("PartWorkbench")
            Gui.activateWorkbench("MeshRemodelWorkbench")
        Gui.runCommand("Part_CheckGeometry",0)
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
                if s.Object.isDerivedFrom("Part::Feature"):
                    self.objs.append(s.Object)
                    count += 1
        if count >= 1:
            return True
        return False

# end check geometry

##################################################################################################

# check geometry
class MeshRemodelSubShapeBinderCommandClass(object):
    """PartDesign SubShapeBinder"""

    def __init__(self):
        self.objs = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'SubShapeBinder.svg') ,
            'MenuText': "PartDesign SubShapeBinder" ,
            'ToolTip' : "Create a PartDesign SubShapeBinder object"}
 
    def Activated(self):
        if not "PartDesign_SubShapeBinder" in Gui.listCommands():
            Gui.activateWorkbench("PartDesignWorkbench")
            Gui.activateWorkbench("MeshRemodelWorkbench")
        doc = FreeCAD.ActiveDocument
        binder = doc.addObject("PartDesign::SubShapeBinder","Binder")
        sel = Gui.Selection.getCompleteSelection()
        support = []
        for s in sel:
            support.append((s.Object, s.SubElementNames))
        binder.Support = support
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
                if s.Object.isDerivedFrom("Part::Feature"):
                    self.objs.append(s.Object)
                    count += 1
        if count >= 1:
            return True
        return False

# end check geometry

##################################################################################################



# initialize

def initialize():
    if FreeCAD.GuiUp:
        Gui.addCommand("MeshRemodelCreatePointsObject", MeshRemodelCreatePointsObjectCommandClass())
        Gui.addCommand("MeshRemodelCreateWireFrameObject",MeshRemodelCreateWireFrameObjectCommandClass())
        Gui.addCommand("MeshRemodelMeshBoundaryWires",MeshRemodelMeshBoundaryWiresCommandClass())
        Gui.addCommand("MeshRemodelMeshSimpleCopy", MeshRemodelMeshSimpleCopyCommandClass())
        Gui.addCommand("MeshRemodelOffsetMesh", MeshRemodelOffsetMeshCommandClass())
        Gui.addCommand("MeshRemodelDuplicateSelectedFacets", MeshRemodelDuplicateSelectedFacetsCommandClass())
        Gui.addCommand("MeshRemodelRemovePoint", MeshRemodelRemovePointCommandClass())
        Gui.addCommand("MeshRemodelAddOrRemoveFacet", MeshRemodelAddOrRemoveFacetCommandClass())
        Gui.addCommand("MeshRemodelMovePoint", MeshRemodelMovePointCommandClass())
        Gui.addCommand("MeshRemodelExpandedMesh",MeshRemodelExpandedMeshCommandClass())
        Gui.addCommand("MeshRemodelCreateCrossSectionsObject",MeshRemodelCreateCrossSectionsCommandClass())
        Gui.addCommand("MeshRemodelAddSelectionObserver",MeshRemodelAddSelectionObserverCommandClass())
        Gui.addCommand("MeshRemodelPartSolid",MeshRemodelPartSolidCommandClass())
        Gui.addCommand("MeshRemodelSubObjectLoft", MeshRemodelSubObjectLoftCommandClass())
        Gui.addCommand("MeshRemodelCreateGridSurface",MeshRemodelCreateGridSurfaceCommandClass())
        Gui.addCommand("MeshRemodelCreatePointObject", MeshRemodelCreatePointObjectCommandClass())
        Gui.addCommand("MeshRemodelCreateCoplanarPointsObject", MeshRemodelCreateCoplanarPointsObjectCommandClass())
        Gui.addCommand("MeshRemodelCreateLine", MeshRemodelCreateLineCommandClass())
        Gui.addCommand("MeshRemodelCreatePolygon", MeshRemodelCreatePolygonCommandClass())
        Gui.addCommand("MeshRemodelCreateBSpline", MeshRemodelCreateBSplineCommandClass())
        Gui.addCommand("MeshRemodelFlattenDraftBSpline", MeshRemodelFlattenDraftBSplineCommandClass())
        Gui.addCommand("MeshRemodelCreatePlane", MeshRemodelCreatePlaneCommandClass())
        Gui.addCommand("MeshRemodelCreateCircle", MeshRemodelCreateCircleCommandClass())
        Gui.addCommand("MeshRemodelCreateArc", MeshRemodelCreateArcCommandClass())
        Gui.addCommand("MeshRemodelCreateWire", MeshRemodelCreateWireCommandClass())
        Gui.addCommand("MeshRemodelMoveAxial", MeshRemodelMoveAxialCommandClass())
        Gui.addCommand("MeshRemodelRotateObject", MeshRemodelRotateObjectCommandClass())
        Gui.addCommand("MeshRemodelGoBackSelection", MeshRemodelGoBackSelectionCommandClass())
        Gui.addCommand("MeshRemodelDraftUpgrade", MeshRemodelDraftUpgradeCommandClass())
        Gui.addCommand("MeshRemodelCreateSketch", MeshRemodelCreateSketchCommandClass())
        Gui.addCommand("MeshRemodelMergeSketches", MeshRemodelMergeSketchesCommandClass())
        Gui.addCommand("MeshRemodelValidateSketch", MeshRemodelValidateSketchCommandClass())
        Gui.addCommand("MeshRemodelPartCheckGeometry", MeshRemodelPartCheckGeometryCommandClass())
        Gui.addCommand("MeshRemodelSubShapeBinder", MeshRemodelSubShapeBinderCommandClass())
        Gui.addCommand("MeshRemodelSettings", MeshRemodelSettingsCommandClass())


initialize()
