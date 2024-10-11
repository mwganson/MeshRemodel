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
__date__    = "2024.10.10"
__version__ = "1.10.29"

import FreeCAD, FreeCADGui, Part, os, math
from PySide import QtCore, QtGui
try:
    from PySide import QtWidgets
except:
    QtWidgets = QtGui
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
        def __init__(self, total=0, txt = ""):
            self.pb = None
            self.btn = None
            self.bar = None
            self.bCanceled = False
            self.value = 0
            self.total = 0
            self.mw = FreeCADGui.getMainWindow()
            self.lastUpdate = time.time()
            if total:
                return self.makeProgressBar(total, buttonText = txt if txt else "Cancel")

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
            if self.mw.isHidden() or self.value >= self.total:
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

    def wireIsPlanar(self, wire):
        """check if the wire is planar and return bool"""
        pts = wire.discretize(5)
        pts = [Part.Vertex(p.x,p.y,p.z) for p in pts[:4]]
        return self.isCoplanar(pts[:3], pts[3].Point, tol=Part.Precision.confusion())


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

    def nearestPoint(self, pt, pts, exclude, skipOne = False):
        """ nearestPoint(pt, pts, exclude)
            pt is a vector, pts a list of vectors
            exclude is a list of vectors to exclude from process
            return nearest point to pt in pts and not in exclude
            if skipOne is true, return 2nd closest point"""
        if len(pts) == 0: #should never happen
            raise Exception("MeshRemodel GeomUtils Error: nearestPoint() pts length = 0\n")
        nearest = pts[0]
        nextNearest = nearest
        d = 10**100
        for p in pts:
            if p in exclude:
                continue
            dis = self.dist(pt, p)
            if dis < d:
                d = dis
                nextNearest = nearest
                nearest = p
        return nearest if not skipOne else nextNearest

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
        pg.RemInt("SketchRadiusPrecision") #no longer used
        keep = pg.GetBool('KeepToolbar',False)
        point_size = pg.GetFloat("PointSize", 4.0)
        line_width = pg.GetFloat("LineWidth", 5.0)
        checkUpdates = pg.GetBool("CheckForUpdates", True)
        pg.SetBool("CheckForUpdates", checkUpdates)
        coplanar_tol = pg.GetFloat("CoplanarTolerance",.01)
        wireframe_tol = pg.GetFloat("WireFrameTolerance",.01)
        items=[("","*")[keep]+"Keep the toolbar active",
            ("","*")[not keep]+"Do not keep the toolbar active",
            ("","*")[checkUpdates]+"Check for updates",
            ("","*")[not checkUpdates]+"Do not check for updates",
            "Change point size ("+str(point_size)+")",
            "Change line width ("+str(line_width)+")",

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
        elif ok and item==items[4]:
            new_point_size,ok = QtGui.QInputDialog.getDouble(window,"Point size", "Enter point size", point_size,1,50,2)
            if ok:
                pg.SetFloat("PointSize", new_point_size)
        elif ok and item==items[5]:
            new_line_width,ok = QtGui.QInputDialog.getDouble(window,"Line width", "Enter line width", line_width,1,50,2)
            if ok:
                pg.SetFloat("LineWidth", new_line_width)
        elif ok and item==items[2]:
            checkUpdates = True
            pg.SetBool("CheckForUpdates", True)
        elif ok and item==items[3]:
            checkUpdates = True
            pg.SetBool("CheckForUpdates", False)
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

Shift+Click to attempt to make faces out of nonplanar wires.  (This can sometimes
take a long time, so save your work first in case you have to force restart.)
"""} 

    def Activated(self):
        copy = self.mesh.Mesh.copy()
        wires = MeshPart.wireFromMesh(copy)
        FreeCAD.Console.PrintMessage(f"MeshRemodel: {len(wires)} wires created\n")
        doc = self.mesh.Document
        modifiers = QtGui.QApplication.keyboardModifiers()
        if modifiers & QtCore.Qt.ShiftModifier:
            makeFilledFaces = True
        else:
            makeFilledFaces = False
        if not wires:
            return
        doc.openTransaction("Mesh boundary wires")
        for idx,wire in enumerate(wires):
            FreeCAD.Console.PrintMessage(f"Making face from wire {idx+1} of {len(wires)}\n")
            FreeCADGui.updateGui()
            obj = doc.addObject("Part::Feature", f"{self.mesh.Label}_BoundaryWire")
            obj.Shape = wire
            plane = wire.findPlane()
            face = None
            if plane:
                face = Part.makeFace(wire,"Part::FaceMakerBullseye")
            elif makeFilledFaces:
                try:
                    face = Part.makeFilledFace(wire.Edges)
                except:
                    face = None
            else:
                FreeCAD.Console.PrintMessage(f"skipping wire {idx+1} of {len(wires)} because it is nonplanar, Shift+CLick to force trying to make filled face (might take a long time, so save your work first.)\n")
            if face:
                mface = MeshPart.meshFromShape(face, 0.25, 50)
                mobj = doc.addObject("Mesh::Feature",f"{self.mesh.Label}_MeshFace")
                mobj.Mesh = mface
        doc.commitTransaction()
        doc.recompute()
        FreeCAD.Console.PrintMessage("Done with Boundary Wires command\n")

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
        self.picked = []

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
        if global_picked:
            self.picked = global_picked
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
            self.mesh.Document.openTransaction("Remove points")
            for pt in self.picked:
                new_mesh = self.removeFacets(self.mesh.Mesh.copy(), pt)
                self.mesh.Mesh = new_mesh
            self.mesh.Document.commitTransaction()

    def removeFacets(self, mesh, pt):
        """remove the point, and all the facets containing the point
        and return the new mesh"""
        if not gu.findPointInMesh(mesh, pt): #removed during previous call
            return mesh
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
        # if len(sel) > 2:
        #     return False
        count = 0
        edgecount = 0
        meshcount = 0
        self.edge = None
        self.pt = None
        self.mid = None
        self.picked = []
        for s in sel:
            if hasattr(s.Object,"Mesh") or s.Object.isDerivedFrom("Mesh::Feature") :
                meshcount = 1
                self.mesh = s.Object
            if s.Object.isDerivedFrom("Part::Feature") and \
                    bool("PointsObject" in s.Object.Name or "WireFrameObject" in s.Object.Name):
                self.mesh = s.Object.Source if s.Object.Source else self.mesh
                if s.Object.Source:
                    meshcount = 1
            if s.SubElementNames and "Vertex" in s.SubElementNames[0]:
                sub = s.Object.getSubObject(s.SubElementNames[0])
                self.pt = sub.Point
                self.picked.append(self.pt)
                count += 1
            if s.SubElementNames and "Edge" in s.SubElementNames[0]:
                sub = s.Object.getSubObject(s.SubElementNames[0])
                self.edge = sub
                edgecount += 1
                self.mid = s.PickedPoints[0]
                
        if count >= 1 and meshcount == 1 and gu.findPointInMesh(self.mesh.Mesh, self.pt):
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
            pv = doc.addObject("Part::Vertex",lineName+"_Mid")
            pv.Placement.Base.x, pv.Placement.Base.y, pv.Placement.Base.z = gu.midpoint(line.firstVertex().Point,line.lastVertex().Point)
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
#TugBoat object
################################################################################
import FreeCAD, FreeCADGui, Part
from PySide import QtCore, QtWidgets, QtGui

class TugBoat:
    def __init__(self, obj):
        obj.Proxy = self
        obj.addExtension("Part::AttachExtensionPython")
        grp = "TugBoat"
        obj.addProperty("App::PropertyLinkHidden", "Vessel", grp, "Tow object to be transported")
        obj.addProperty("App::PropertyLinkHidden", "Destination", grp,"Empty = Origin, global or local depending on GlobalPlacement setting.  Otherwise can be any object with a placement property, preferably another TugBoat serving as destination point.  After the transport this TugBoat will align itself perfectly with destination TugBoat and apply the same transformation to the vessel.")
        obj.addProperty("App::PropertyEnumeration","Commands", grp, "Sit tight = do nothing, Haul vessel = do transport.  Note that after doing the transport this gets set back to Sit tight automatically.")
        obj.Commands = ["Sit tight", "Haul vessel","Rotate only","Translate only"]
        grp = "Tugboat visual"
        obj.addProperty("App::PropertyFloatConstraint","CircleScale",grp,"Scales the circles").CircleScale = (1,0.1,float("inf"),.01)
        obj.addProperty("App::PropertyFloatConstraint","LineScale",grp,"Scales the line segment").LineScale = (1,0.1,float("inf"),.01)
        obj.addProperty("App::PropertyBool","OffsetLines",grp, "Whether to have a small gap between center and axis lines").OffsetLines = True
        grp = "TugBoat settings"
        obj.addProperty("App::PropertyBool","GlobalPlacement",grp,"Whether to use global or local coordinate system, tow vessel to local origin or global origin, useful where all objects are not in the same container.").GlobalPlacement = True
        obj.addProperty("App::PropertyBool","EditAttachmentOffsetWhenMoving",grp,"Whether to adjust the attachment offset of an object in order to complete a haul/move.").EditAttachmentOffsetWhenMoving = True
        obj.addProperty("App::PropertyBool","AdjustProfileWhenMoving",grp,"Whether to move an unattached profile, such as a subshapebinder or an unattached sketch when hauling a part design feature").AdjustProfileWhenMoving = True
        obj.addProperty("App::PropertyBool","AdjustSketchAttachmentOffsetWhenMoving",grp,"When moving a PartDesign::Feature object if it is profile based, and if it has a sketch as the profile, and if the sketch is attached, then if this is true, we adjust the sketch's attachment offset property to position the feature.").AdjustSketchAttachmentOffsetWhenMoving = True
        obj.addProperty("App::PropertyString","Version",grp,"Version used to create this object").Version = __version__


    def onChanged(self, fp, prop):
        if prop == "Commands":
            if fp.Commands == "Haul vessel" or fp.Commands == "Rotate only" or fp.Commands == "Translate only":
                self.haulIt(fp)
        elif bool(prop == "Vessel" or prop == "Destination" or prop == "Placement") and fp.ViewObject:
            fp.ViewObject.signalChangeIcon()

    def isAttached(self, objToCheck):
        """checks to see if objToCheck is attached"""

        if not objToCheck.hasExtension('Part::AttachExtension'):
            return False
        supportPropertyName = "AttachmentSupport" if hasattr(objToCheck,"AttachmentSupport") else "Support"
        support = getattr(objToCheck, supportPropertyName)
        mapMode = objToCheck.MapMode
        if mapMode == "Deactivated":
            return False
        if support:
            return True

    def adjustProfile(self, fp):
        """called when a part design feature is being moved and its profile object is not attached
        this function will then move the profile object"""
        if not fp.AdjustProfileWhenMoving:
            return
        if not fp.Vessel.isDerivedFrom("PartDesign::Feature"):
            return
        if not hasattr(fp.Vessel, "Profile"):
            return
        if not fp.Vessel.Profile:
            return
        sk,subnames = fp.Vessel.Profile #sk might be sketch or binder or None
        if not sk:
            return
        if sk.isDerivedFrom("PartDesign::Feature"):
            FreeCAD.Console.PrintWarning(f"{fp.Label} captain: We will adjust the placement of {sk.Label} as part of hauling {fp.Vessel.Label}, but it's possible this might be unsuccessful.  We can move unnattached sketches or SubShapeBinders like this when they're used as profiles, but part design features should be hauled separately.")
        sk.Document.openTransaction(f"{fp.Label} haul {sk.Label}")
        sk.Placement = fp.Vessel.Placement #vessel has already been moved
        sk.Document.commitTransaction()

    def adjustSketchAttachmentOffset(self, fp):
        """ if fp.Vessel is a PartDesign feature and if it has a Profile property
        and if the profile is a sketch and if the sketch is attached and if
        fp.AdjustSketchAttachmentOffsetWhenMoving is True, then we move the object
        by adjusting its sketch's attachment offset property
        """
        if not fp.AdjustSketchAttachmentOffsetWhenMoving:
            return
        if not fp.Vessel.isDerivedFrom("PartDesign::Feature"):
            return
        if not hasattr(fp.Vessel, "Profile"):
            return
        if not fp.Vessel.Profile:
            return
        sk,subnames = fp.Vessel.Profile
        if not sk:
            return
        if not sk.isDerivedFrom("Sketcher::SketchObject"):
            return
        if not self.isAttached(sk):
            return
        oldPlm = fp.Vessel.Placement #vessel has already been moved

        att_plm = sk.Attacher.calculateAttachedPlacement(sk.Placement).multiply(sk.AttachmentOffset.inverse())

        att_offset = att_plm.inverse().multiply(oldPlm)
        sk.Document.openTransaction(f"Haul {fp.Vessel.Label} : {sk.Label} attachment offset")
        sk.AttachmentOffset = att_offset
        sk.Document.commitTransaction()

    def adjustAttachmentOffset(self, fp):
        """ adjusts the attachment offset of fp.Vessel in order to make it stay put,
        """
        if not self.isAttached(fp.Vessel) or not fp.EditAttachmentOffsetWhenMoving:
            return

        # credit to DeepSOIC for the basis of the following code, from a macro
        # he wrote to modify attachment offsets to keep a sketch  in its current position
        # I have modified his code to meet my needs here, any mistakes are mine, not his
        # since this is now being used in a different context

        oldPlm = fp.Vessel.Placement

        att_plm = fp.Vessel.Attacher.calculateAttachedPlacement(fp.Vessel.Placement).multiply(fp.Vessel.AttachmentOffset.inverse())

        att_offset = att_plm.inverse().multiply(oldPlm)
        fp.Document.openTransaction(f"Edit {fp.Vessel.Label} attachment offset")
        fp.Vessel.AttachmentOffset = att_offset
        fp.Document.commitTransaction()


    def execute(self, fp):
        radii = .5
        radii *= fp.CircleScale if fp.CircleScale else 1
        lineLength = 2.5
        lineLength *= fp.LineScale if fp.LineScale else 1
        pnt = FreeCAD.Vector()
        axisZ = FreeCAD.Vector(0,0,1)
        axisX = FreeCAD.Vector(1,0,0)
        axisY = FreeCAD.Vector(0,1,0)
        vert = Part.Vertex(pnt)
        circleZ = Part.makeCircle(radii, pnt, axisZ)
        circleY = Part.makeCircle(radii, pnt, axisY)
        circleX = Part.makeCircle(radii, pnt, axisX)
        offset = lineLength / 20 if fp.OffsetLines else 0
        upLine = Part.makePolygon([FreeCAD.Vector(0, 0, offset), FreeCAD.Vector(0, 0, lineLength)])
        rightLine = Part.makePolygon([FreeCAD.Vector(offset, 0, 0), FreeCAD.Vector(lineLength, 0, 0)])
        backLine = Part.makePolygon([FreeCAD.Vector(0, offset, 0), FreeCAD.Vector(0, lineLength, 0)])
        comp = Part.makeCompound([vert, circleX, circleY, circleZ, rightLine, backLine, upLine])
        fp.Shape = comp
        fp.ViewObject.LineColorArray = [(255,0,0,255),(0,255,0,255),(0,0,255,255)] * 2 #red, green, blue
        fp.ViewObject.LineWidth = 3
        fp.ViewObject.PointSize = 5
        fp.ViewObject.PointColor = (85, 255, 255)
        fp.positionBySupport()

    def getRecursiveGlobalPlm(self, obj):
        """ app::part containers might be nested, so recursively calculate the global placement
        obj is assumed to be a valid object that is the parent container of another object, for
        example it could be a PartDesign::Body or an App::Part, returns the product of all the
        parent object's nested placements
        """
        assert(obj is not None) #programming logical error if obj is None
        full = obj.Placement.copy()
        parent = obj.getParentGeoFeatureGroup()
        while parent:
            full = parent.Placement.multiply(full)
            parent = parent.getParentGeoFeatureGroup()
        return full

    def haulIt(self, fp):
        """do the actual move"""

        if fp.Vessel:
            if not hasattr(fp.Vessel, "Placement"):
                FreeCAD.Console.PrintError(f"{fp.Label} captain: Nothing to latch onto!  Can't move {fp.Vessel.Label} because it has no Placement property.\n")
                return

            if fp.Vessel.isDerivedFrom("App::Plane")or fp.Vessel.isDerivedFrom("App::Line") or fp.Vessel.isDerivedFrom("App::Origin"):
                FreeCAD.Console.PrintError(f"{fp.Label} captain: She won't budge!  Not moving origin object, select the container object instead, e.g. Part Design Body or App::Part object.\n")
                return

            plmDest = FreeCAD.Placement()
            if fp.Destination and not fp.GlobalPlacement:
                vParent = fp.Vessel.getParentGeoFeatureGroup()
                dParent = fp.Destination.getParentGeoFeatureGroup()
                if vParent.Name != dParent.Name:
                    FreeCAD.Console.PrintWarning(f"{fp.Label} captain: Since the vessel is in a different container than the destination, you might want to prefer setting GlobalPlacement to True.\n")

            if fp.Destination and hasattr(fp.Destination, "Placement"):
                if fp.GlobalPlacement:
                    plmDest = self.getRecursiveGlobalPlm(fp.Destination).copy()
                else:
                    plmDest = fp.Destination.Placement.copy()

            if fp.GlobalPlacement:
                fpPlm = self.getRecursiveGlobalPlm(fp)
                vesselPlm = self.getRecursiveGlobalPlm(fp.Vessel)
            else:
                fpPlm = fp.Placement.copy()
                vesselPlm = fp.Vessel.Placement.copy()

            relativePlm = fpPlm.inverse().multiply(vesselPlm) #full transformation if unchanged below

            if fp.Commands == "Translate only":
                newPlm = fpPlm.copy()
                newPlm.move(plmDest.Base - fpPlm.Base)
            elif fp.Commands == "Rotate only":
                newPlm = fpPlm.copy()
                newPlm.Rotation = vesselPlm.Rotation.slerp(plmDest.Rotation, 1) #1 means complete rotation
            else:
                newPlm = plmDest.multiply(relativePlm) #default full transformation

            #bug in isSame()? -- use the string comparison as a double check
            if not newPlm.isSame(vesselPlm, Part.Precision.confusion()) or f"{newPlm}" != f"{vesselPlm}":
                fp.Document.openTransaction(f"TugBoat haul {fp.Vessel.Label}")
                if fp.GlobalPlacement:
                    vParent = fp.Vessel.getParentGeoFeatureGroup()
                    if vParent:
                        fp.Vessel.Placement = self.getRecursiveGlobalPlm(vParent).inverse().multiply(newPlm)
                    else:
                        fp.Vessel.Placement = newPlm
                else:
                    fp.Vessel.Placement = newPlm

                if not fp.Vessel.isDerivedFrom("PartDesign::Feature"):
                    self.adjustAttachmentOffset(fp) # checks user EditAttachmentOffsetWhenMoving boolean
                else:
                    self.adjustSketchAttachmentOffset(fp) #simply returns if profile is attached
                    self.adjustProfile(fp) #handles unattached profiles

                if fp.GlobalPlacement:
                    # we must undo the parent container's global transformation here if necessary
                    # check if we have a parent, if so, get the parent's global placement
                    fpParent = fp.getParentGeoFeatureGroup()
                    if fpParent:
                        fp.Placement = self.getRecursiveGlobalPlm(fpParent).inverse().multiply(plmDest)
                    else:
                        fp.Placement = plmDest
                else:
                    fp.Placement = plmDest
                fp.Document.commitTransaction()
        else:
            FreeCAD.Console.PrintError(f"{fp.Label} captain: No vessel.  Set your vessel and try again.\n")

        fp.Commands = "Sit tight"
        fp.Document.recompute()


class TugBoatVP:
    def __init__(self, vobj):
        vobj.Proxy = self
        self.fpName = vobj.Object.Name

    def doubleClicked(self, vobj):
        fp = vobj.Object
        fp.Commands = "Haul vessel"
        fp.Document.recompute()

    def setupContextMenu(self, vobj, menu):
        def makePutTrigger(arg): return lambda : self.putInContainer(arg)
        def makePickTrigger(arg): return lambda : self.pickDestination(arg)
        def makeVesselTrigger(arg): return lambda : self.pickVessel(arg)

        fp = vobj.Object

        if fp.Vessel:
            haulMenu = menu.addMenu("Haul commands")
            haulAction = haulMenu.addAction(f"Haul {fp.Vessel.Label}")
            haulAction.triggered.connect(lambda : setattr(fp, "Commands", "Haul vessel"))
            translateAction = haulMenu.addAction(f"Translate only {fp.Vessel.Label}")
            translateAction.triggered.connect(lambda : setattr(fp, "Commands", "Translate only"))
            rotateAction = haulMenu.addAction(f"Rotate only {fp.Vessel.Label}")
            rotateAction.triggered.connect(lambda : setattr(fp, "Commands", "Rotate only"))
            if fp.Vessel.TypeId == "Part::Feature" or fp.Vessel.TypeId == "Mesh::Feature" or fp.Vessel.TypeId == "Points::Feature":
                identityAction = menu.addAction(f"Set {fp.Vessel.Label} to identity placement")
                identityAction.triggered.connect(lambda : self.setToIdentity(FreeCAD.Placement()))
                fpIdentityAction = menu.addAction(f"Set {fp.Vessel.Label} to {fp.Label} placement")
                fpIdentityAction.triggered.connect(lambda : self.setToIdentity(fp.Placement.inverse()))
        else:
            haulAction = menu.addAction("Haul vessel (no vessel configured) ")
            haulAction.setEnabled(False)

        fpParent = fp.getParentGeoFeatureGroup()
        containers = [obj for obj in fp.Document.Objects if obj.isDerivedFrom("App::Part") or obj.isDerivedFrom("PartDesign::Body")]
        if containers:
            putMenu = menu.addMenu("Put into container")
            for container in containers:
                star = "*" if container == fpParent else ""
                putAction = putMenu.addAction(f"{star}{container.Label}")
                putAction.setIcon(container.ViewObject.Icon)
                putAction.triggered.connect(makePutTrigger(container))

        tugs = [obj for obj in fp.Document.Objects if "TugBoat" in obj.Name and obj != fp]
        pickMenu = menu.addMenu("Pick destination")
        if not fp.Destination and not fp.GlobalPlacement:
            pickAction = pickMenu.addAction("*Local origin")
            pickAction.setEnabled(False)
        else:
            pickAction = pickMenu.addAction("Local origin")
        pickAction.triggered.connect(makePickTrigger("Local origin"))
        if not fp.Destination and fp.GlobalPlacement:
            pickAction = pickMenu.addAction("*Global origin")
            pickAction.setEnabled(False)
        else:
            pickAction = pickMenu.addAction("Global origin")
        pickAction.triggered.connect(makePickTrigger("Global origin"))

        pickAction = pickMenu.addAction("New TugBoat")
        pickAction.triggered.connect(makePickTrigger("New TugBoat"))

        for tug in tugs:
            if fp.Destination == tug:
                pickAction = pickMenu.addAction(f"*{tug.Label}")
                pickAction.setEnabled(False)
            else:
                pickAction = pickMenu.addAction(tug.Label)
            pickAction.setIcon(tug.ViewObject.Icon)
            pickAction.triggered.connect(makePickTrigger(tug))

        vessels = [obj for obj in fp.Document.Objects if
                    hasattr(obj,"Placement")
                    and obj!= fp and obj != fp.Destination
                    and not obj.isDerivedFrom("App::Line")
                    and not obj.isDerivedFrom("App::Plane")
                    and not obj.isDerivedFrom("App::Origin")]
        if vessels:
            vesselMenu = menu.addMenu("Pick vessel")
            for vessel in vessels:
                if fp.Vessel == vessel:
                    pickVesselAction = vesselMenu.addAction(f"*{vessel.Label}")
                    pickVesselAction.setEnabled(False)
                else:
                    pickVesselAction = vesselMenu.addAction(vessel.Label)
                pickVesselAction.setIcon(vessel.ViewObject.Icon)
                pickVesselAction.triggered.connect(makeVesselTrigger(vessel))

    def setToIdentity(self, newPlm):
        """
        Here we attempt to reset the internal placement of the object by transforming
        its geometry.  If newPlm is the identity placement, then the object is at
        the origin and the user wishes to set this as the new placement, typically.
        Otherwise, the vessel's internal placement is set to the placement of the TugBoat.
        """
        fp = FreeCAD.ActiveDocument.getObject(self.fpName)
        vessel = fp.Vessel
        if not fp.Vessel:
            FreeCAD.Console.PrintError(f"{fp.Label} captain: Can't reset internal placement of nothing.\n")
            return
        fp.Document.openTransaction(f"set new placement {fp.Vessel.Label}")
        if vessel.TypeId == "Points::Feature":
            pts = vessel.Points
            if newPlm == FreeCAD.Placement():
                pts = vessel.Points.copy()
                points.transformGeometry(vessel.Placement.toMatrix())
            else:
                pts = vessel.Points.copy()
                pts.transformGeometry(newPlm.toMatrix())
            vessel.Points = pts
            vessel.Placement = newPlm.inverse()
            fp.Document.recompute()
        elif vessel.TypeId == "Mesh::Feature":
            if newPlm == FreeCAD.Placement():
                mesh = vessel.Mesh.copy()
                mesh.transform(vessel.Placement.toMatrix())
            else:
                mesh = vessel.Mesh.copy()
                mesh.transform(newPlm.toMatrix())
            vessel.Mesh = mesh
            vessel.Placement = newPlm.inverse()
        elif vessel.TypeId == "Part::Feature":
            if newPlm == FreeCAD.Placement():
                shape = vessel.Shape.copy().transformShape(vessel.Placement.toMatrix(), True)
            else:
                shape = vessel.Shape.copy().transformShape(newPlm.toMatrix(), True)
            vessel.Shape = shape
            vessel.Placement = newPlm.inverse()
        else:
            FreeCAD.Console.PrintError(f"{fp.Label} captain: Unsupported vessel type for resetting internal placement, can't work with objects of type: {vessel.TypeId}\n")
            fp.Document.abortTransaction()
        fp.Document.commitTransaction()
        fp.Document.recompute()


    def pickVessel(self, vessel):
        fp = FreeCAD.ActiveDocument.getObject(self.fpName)
        fp.Vessel = vessel
        fp.Document.recompute()

    def pickDestination(self, destination):
        fp = FreeCAD.ActiveDocument.getObject(self.fpName)

        if destination == "Local origin":
            fp.Destination = None
            fp.GlobalPlacement = False
        elif destination == "Global origin":
            fp.Destination = None
            fp.GlobalPlacement = True
        elif destination == "New TugBoat":
            tug = fp.Document.addObject("Part::FeaturePython","TugBoat")
            TugBoat(tug)
            TugBoatVP(tug.ViewObject)
            fp.Destination = tug
            fpParent = fp.getParentGeoFeatureGroup()
            if fpParent:
                fpParent.Group = fpParent.Group + [tug]
        else:
            fp.Destination = destination
        fp.Document.recompute()

    def putInContainer(self, container):
        fp = FreeCAD.ActiveDocument.getObject(self.fpName)
        if fp in container.Group:
            return
        container.Group = container.Group + [fp]
        fp.Document.recompute()

    def getIcon(self):
        fp = FreeCAD.ActiveDocument.getObject(self.fpName)
        icon = TugBoatVP.__icon__
        if fp.Vessel:
            icon = icon.replace("none","darkgray")
            icon = icon.replace("red","blue")
        if fp.Destination:
            icon = icon.replace("burgundy", "green")
            if "TugBoat" in fp.Destination.Name:
                if fp.Destination.ViewObject:
                    fp.Destination.ViewObject.signalChangeIcon()
        tugs = [obj for obj in fp.Document.Objects if "TugBoat" in obj.Name and fp == obj.Destination]
        if tugs:
            icon = icon.replace("none","orange")
            for tug in tugs:
                if tug.ViewObject:
                    tug.ViewObject.signalChangeIcon()
        return icon

    __icon__ = """/* XPM */
static char *dummy[]={
"64 64 4 1",
". c none",
"# c burgundy",
"b c red",
"a c None",
"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"..............................aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"..............................aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"...............................aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"...............................aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"...............................aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaa.........................aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaa..................aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaa.................aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaa................aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaa..............aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaa..............aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaa.............aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaa...........aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaa............aaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaa...........aaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaa...........aaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaa..........aaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaa..........aaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaa..........aaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaaa.........aaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaaa.........aaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaaaa........aaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaaaaa.......aaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaaaaa.......aaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaaaaa.......aaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaaaaa.......aaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaa#################aaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaa###################aaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaa#####################aaaaaaaaaaaa#aaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaa#####################aaaaaaaaaaa###aaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaa###bbbbbbbbbbbbbbb###aaaaaaaaaa#####aaaaaa",
"aaaaaaaaaaaaaaaaaaaaaa###bbbbbbbbbbbbbbb###aaaaaaaaaa#####aaaaaa",
"aaaaaaaaaaaaaaaaaaaaaa###bbbbbbbbbbbbbbb###aaaaaaaaa######aaaaaa",
"aaaaaaaaaaaaaaaaaaaaaa###bbbbbbbbbbbbbbb###aaaaaaaaa######aaaaaa",
"aaaaaaaaaaaaaaaaaaaaaa###bbbbbbbbbbbbbbb###aaaaaaaa#######aaaaaa",
"aaaaaaaaaaaaaaaaaaaaaa###bbbbbbbbbbbbbbb###aaaaaaa########aaaaaa",
"aaaaaaaaaaaaaaaaaaaaaa###bbbbbbbbbbbbbbb###aaaaaa#########aaaaaa",
"aaaaaaaaaaaaaaaaaaaaaa###bbbbbbbbbbbbbbb###aaaaa##########aaaaaa",
"aaaaaaaaaaaaaaaaaaaaaa###bbbbbbbbbbbbbbb###aaaa###bbbbb###aaaaaa",
"aaaaaaaaaaaaa############bbbbbbbbbbbbbbb###aa#####bbbbb###aaaaaa",
"aaaaaaaaaaaa#############bbbbbbbbbbbbbbb###a#####bbbbbb###aaaaaa",
"aaaaaaaaaaa##############bbbbbbbbbbbbbbb########bbbbbbb###aaaaaa",
"aaaaaaaaaaa#############bbbbbbbbbbbbbbbb######bbbbbbbbb###aaaaaa",
"aaaaaaaaaaa############bbbbbbbbbbbbbbbbb#####bbbbbbbbbb###aaaaaa",
"aaaaaaaaaaa###bbbbbbbbbbbbbbbbbbbbbbbbbb###bbbbbbbbbbbb###aaaaaa",
"aaaaaaaaaaa###bbbbbbbbbbbbbbbbbbbbbbbbbbb#bbbbbbbbbbbbb###aaaaaa",
"aaaaaaaaaaa###bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb####aaaaaaa",
"aaaaaaaaaaa###bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb#####aaaaaaaa",
"aaaaaaaaaaa###bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb######aaaaaaaaa",
"aaaaaaaaaaa###bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb#######aaaaaaaaaa",
"aaaaaaaaaaa###bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb########aaaaaaaaaaa",
"aaaaaaaaaaa###bbbbbbbbbbbbbbbbbbbbbbbbbbbb#########aaaaaaaaaaaaa",
"aaaaaaaaaaa###bbbbbbbbbbbbbbbbbbbbbbbbb##########aaaaaaaaaaaaaaa",
"aaaaaaaaaaa###bbbbbbbbbbbbbbbbbbbb#############aaaaaaaaaaaaaaaaa",
"aaaaaaaaaaa###################################aaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaa################################aaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaa############################aaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaa########################aaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaa##################aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"};
"""

class MeshRemodelCreateTugBoatCommandClass:
    """TugBoat object for hauling objects around"""
    
    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'TugBoat.svg') ,
            'MenuText': "TugBoat" ,
            'ToolTip' : \
"""Use the tugboat to haul other objects around in FreeCAD.  Works on any type
of object as long as it has a Placement property, including meshes, points,
and shapes.

Usage: Create the TugBoat object.  Use the context menu to configure the vessel
and the destination.  If destination is not set, then the global origin is the
default.  Move the tugboat to the position on the vessel where you want it to 
latch onto.  The tugboat can be attached to the vessel or not, should work either
way.  Use the haul vessel option from the context menu or doubleclick the tugboat
to do the haul.
"""}

    def Activated(self):
        
        doc = FreeCAD.ActiveDocument if FreeCAD.ActiveDocument else FreeCAD.newDocument()

        doc.openTransaction("Create TugBoat")
        fp = doc.addObject("Part::FeaturePython","TugBoat")
        TugBoat(fp)
        TugBoatVP(fp.ViewObject)
        doc.commitTransaction()
        sel = Gui.Selection.getSelection()
        if len(sel) == 1:
            if "TugBoat" in sel[0].Name:
                #duplicate all properties
                tug = sel[0]
                blacklist = ["PropertiesList","ExpressionEngine","Label","Label2"]
                props = [p for p in tug.PropertiesList if not p in blacklist]
                for prop in props:
                    setattr(fp, prop, getattr(tug, prop))
            elif hasattr(sel[0], "Placement"):
                fp.Vessel = sel[0]
            vParent = sel[0].getParentGeoFeatureGroup()
            if vParent:
                vParent.Group = vParent.Group + [fp]
        doc.recompute()

    def IsActive(self):
        if FreeCAD.ActiveDocument:
            return True
        return False

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
        def report_trigger(txt):
            return lambda: self.on_report(txt)
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
            button_row_layout = QtGui.QHBoxLayout()
            
            button = QtGui.QPushButton(f'Block {i + 1}')
            tooltip_text = "\n".join(block)
            button.setToolTip(tooltip_text)
            button.clicked.connect(trigger(i+1))
            button_row_layout.addWidget(button)
            
            report_button = QtGui.QPushButton(f"> text document")
            report_button.setToolTip("send contents to a new text document object")
            report_button.clicked.connect(report_trigger(tooltip_text))
            button_row_layout.addWidget(report_button)
            
            
            scroll_layout.addLayout(button_row_layout)

        scroll_content.setLayout(scroll_layout)
        scroll_area.setWidget(scroll_content)

        layout.addWidget(scroll_area)
        self.setLayout(layout)

    def on_block_selected(self, index):
        self.selected_index = index
        self.accept()  # Close the dialog and set the result to Accepted
        
    def on_report(self, txt):
        txt = "# Copy / paste this into the python console to restore this selection block\n" + txt
        obj = FreeCAD.ActiveDocument.addObject("App::TextDocument", "Selection block")
        obj.Text = txt


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

Click the > text document button to send the block to a new text document object.  
This is useful where you wish to resume some work in another session and want this
selection block available for you to use.  Use it by copy/paste into the python 
console.  This provides persistent storage for selection blocks.

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
        blocks = self.parseTextDocuments() + blocks
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

    def parseTextDocuments(self):
        """Check for any text documents created previously and return those as blocks"""
        doc = FreeCAD.ActiveDocument
        if not doc:
            return []
        textDocs = [obj for obj in doc.Objects if obj.isDerivedFrom("App::TextDocument") and "Selection_block" in obj.Name]
        blocks = []
        for td in textDocs:
            lines = td.Text.split("\n") if td.Text else []
            lines.reverse() #because make_blocks will reverse again
            blocks.extend(self.make_blocks(lines))
        return blocks
            
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
            if hasattr(fp.ViewObject.Proxy,"onDelete"):
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
        def rowsToCols(vrows):
            if not vrows:
                return []
            
            allSameLen = all(len(row) == len(vrows[0]) for row in vrows)
            if allSameLen:
                vcols = [[vrows[row][col] for row in range(len(vrows))] for col in range(len(vrows[0]))]
            else:
                #first and last column only if not all are the same length
                vcols = [
                    [row[0] for row in vrows],
                    [row[-1] for row in vrows if row]
                ]
            return vcols

        splines = []
        vcols = rowsToCols(vrows)

        for row in vrows:
            if len(row) == 1:
                spline = Part.Vertex(row[0])
                splines.append(spline)
                continue
            spline = Part.BSplineCurve()
            if "BSpline" in fp.RowWireType:
                spline.interpolate(row, PeriodicFlag = "Closed BSpline" == fp.RowWireType)
            elif "Open Polyline" == fp.RowWireType:
                spline = Part.makePolygon(row)
            elif "Closed Polyline" in fp.RowWireType:
                spline = Part.makePolygon(row + [row[0]])
            splines.append(spline.toShape() if hasattr(spline, "toShape") else spline)

        if fp.Output == "Gordon Template":
            colSplines = []
            for col in vcols:
                colSpline = Part.BSplineCurve()
                if "BSpline" in fp.RowWireType:
                    colSpline.interpolate(col) if col else []
                    if colSpline:
                        colSpline = colSpline.toShape()
                elif "Polyline" in fp.RowWireType:
                    colSpline = Part.makePolygon(col) if col else []
                colSplines.append(colSpline)
                
            splines = colSplines + splines

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
                if len(row) == 1:
                    spline = Part.Vertex(row[0])
                    splines.append(spline)
                    if row != vrows[0] and row != vrows[-1]:
                        FreeCAD.Console.PrintWarning(f"{fp.Label}: Lofts will fail if a vertex is used in the middle instead of starting or ending wire\n")
                    continue
                if "BSpline" in fp.RowWireType:
                    spline = Part.BSplineCurve()
                    periodicFlag = "Closed" in fp.RowWireType
                    try:
                        spline.interpolate(row, PeriodicFlag = periodicFlag)
                    except:
                        FreeCAD.Console.PrintError(f"GridSurface: Error interpolating spline from row = {row}\n")
                        spline = Part.Shape()
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
                if splines:
                    shape = splines[0].toShape() if hasattr(splines[0],"toShape") else splines[0]
                    if fp.Output != "Surface":
                        fp.Shape = shape
                        return
                    else:
                        #try to make a face since there is only 1 row and output type is surface
                        if not shape.isNull() and shape.isClosed():
                            if gu.wireIsPlanar(shape):
                                face = Part.makeFace(shape, "Part::FaceMakerBullseye")
                                fp.Shape = face
                                return
                            else:
                                # FreeCAD.Console.PrintMessage(f"{fp.Label}: making filled face from nonplanar wire, check for defects.\n")
                                face = Part.makeFilledFace(shape.Edges)
                                fp.Shape = face
                                return
                        else: #not closed
                            fp.Shape = shape
                            return
                else:
                    fp.Shape = Part.Shape()
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
            if not hasattr(vert, "Row"):
                vert = self.insert_vertex(fp, vert)
                vert.recompute()
                if not vert:
                    continue
            row,col = (vert.Row, vert.Column)
            if row == -1 and col == -1:
                continue
            if hasattr(vert,"GS_Ignore") and vert.GS_Ignore:
                continue
            elif not hasattr(vert,"GS_Ignore"):
                vert.addProperty("App::PropertyBool","GS_Ignore","GridSurface").GS_Ignore = False
                vert.recompute()
            vert_dict[vert.Name] = (vert,row,col)
        return vert_dict

    def insert_vertex(self, fp, vert):
        """called by buildDictionary() when a new vertex has been inserted into the Vertices property"""
        if not hasattr(vert, "Row"):
            vert.addProperty("App::PropertyInteger","Row","GridSurface")
            vert.Row = -1
        if not hasattr(vert,"Column"):
            vert.addProperty("App::PropertyInteger","Column","GridSurface")
            vert.Column = -1
        if not hasattr(vert,"GS_Ignore"):
            vert.addProperty("App::PropertyBool","GS_Ignore","GridSurface").GS_Ignore = False
        val=0
        FreeCADGui.Selection.clearSelection()
        FreeCADGui.Selection.addSelection(fp.Document.Name, vert.Name)
        row,ok = QtGui.QInputDialog.getInt(FreeCADGui.getMainWindow(), "Row input", f"Enter row for this vertex {vert.Label} at {(vert.X, vert.Y, vert.Z)}",\
                                    val, minValue = 0, step=1)
        if not ok:
            return vert
        else:
            col,ok = QtGui.QInputDialog.getInt(FreeCADGui.getMainWindow(), "Column input", f"Enter column for this vertex {vert.Label} at {(vert.X, vert.Y, vert.Z)}",\
                                    val, minValue = 0, step=1)
        if ok:
            if not fp.ColumnMode:
                row_list = [v for v in fp.Vertices if hasattr(v, "Row") and hasattr(v, "Column") and v.Column >= col and v != vert]
                highest_col = 0
                items = ["yes, shift them all","yes, but only this row", "no, do not shift"]
                item, ok = QtGui.QInputDialog.getItem(FreeCADGui.getMainWindow(),"Shift vertices?", \
                            "Shift remaining vertices over one column?",items,0)
                if ok and item != items[2]:
                    for r in row_list:
                        r.Column += 1 if item == items[0] or r.Row == row else 0
                        if r.Column > highest_col:
                            highest_col = r.Column

                vert.Row = row
                vert.Column = col
                self.blockSignals = True
                fp.CountColumns += 1 if highest_col >= fp.CountColumns else 0
                self.blockSignals = True
                fp.CountRows += 1 if row >= fp.CountRows else 0

            else:
                column_list = [v for v in fp.Vertices if hasattr(v, "Row") and hasattr(v, "Column") and v != vert and v.Row >= row]
                highest_row = 0
                items = ["yes, shift them all","yes, but only this column", "no, do not shift"]
                item, ok = QtGui.QInputDialog.getItem(FreeCADGui.getMainWindow(),"Shift vertices?", \
                            "Shift remaining vertices over one row?",items,0)
                if ok and item != items[2]:
                    for c in column_list:
                        c.Row += 1 if item == items[0] or c.Column == col else 0
                        if c.Row > highest_row:
                            highest_row = c.Row
                vert.Row = row
                vert.Column = col
                self.blockSignals = True
                fp.CountRows += 1 if highest_row >= fp.CountRows else 0
                self.blockSignals = True
                fp.CountColumns += 1 if col >= fp.CountColumns else 0
                
        return vert
        


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
        vert.addProperty("App::PropertyBool","GS_Ignore","GridSurface").GS_Ignore = False
        # vert.setEditorMode("Row",1)
        # vert.setEditorMode("Column",1)
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

        self.ghostName = None # the name of the ghost object


    @property
    def FP(self):
        return FreeCAD.ActiveDocument.getObject(self.name) if hasattr(self,"name") else FreeCAD.ActiveDocument.getObject("GridSurface")

    def setupContextMenu(self, vobj, menu):
        def selectRowTrigger(row): return lambda: self.selectRow(row)
        def selectColTrigger(col): return lambda: self.selectCol(col)
        def outputMode(mode): return lambda: self.switchOutput(mode)
        def methodMode(mode): return lambda: self.switchMethod(mode)
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

        methodMenu = menu.addMenu("Switch method")
            
        if fp.Method != "Interpolate Points":
            methodAction = methodMenu.addAction("Interpolate Points")
            methodAction.triggered.connect(methodMode("Interpolate Points"))
            
        if fp.Method != "Lofted Edges":
            methodAction = methodMenu.addAction("Lofted Edges")
            methodAction.triggered.connect(methodMode("Lofted Edges"))
            
        if fp.Method != "Lofted Edges Ruled":
            methodAction = methodMenu.addAction("Lofted Edges Ruled")
            methodAction.triggered.connect(methodMode("Lofted Edges Ruled"))

        visibilityAction = menu.addAction("Toggle visibility of vertices")
        visibilityAction.triggered.connect(self.toggleVisibility)

        mode = "Row Mode" if fp.ColumnMode else "Column Mode"
        columnModeAction = menu.addAction(f"Switch to {mode}")
        columnModeAction.triggered.connect(self.toggleColumnMode)
        
        body = FreeCADGui.ActiveDocument.ActiveView.getActiveObject("pdbody")
        if body and not fp in body.Group:
            bodyAction = menu.addAction(f"Add to {body.Label}")
            bodyAction.triggered.connect(self.addToBody)
        elif body and fp in body.Group:
            bodyAction = menu.addAction(f"Remove from {body.Label}")
            bodyAction.triggered.connect(self.removeFromBody)
            
        ghost = fp.Document.getObject(self.ghostName) if self.ghostName else None
        if not ghost:
            transformAction = menu.addAction("Begin transforming grid")
            transformAction.triggered.connect(self.beginTransformingGrid)
        else:
            transformAction = menu.addAction("Finish transforming grid")
            transformAction.triggered.connect(self.finishTransformingGrid)
            
            cancelAction = menu.addAction("Cancel transforming grid")
            cancelAction.triggered.connect(self.cancelTransformingGrid)
            

        menu.addSeparator()

    def beginTransformingGrid(self, checked):
        """begin transform grid"""
        fp = self.FP
        vert_dict = fp.Proxy.buildDictionary(fp)
        verts = []
        for k,v in vert_dict.items():
            verts.append(v[0].Shape.copy())
        comp = Part.makeCompound(verts)
        ghost = fp.Document.getObject(f"{fp.Name}_Ghost")
        if ghost:
            fp.Document.removeObject(ghost.Name)
        ghost = Part.show(comp,f"{fp.Name}_Ghost")
        self.ghostName = ghost.Name
        ghost.addProperty("App::PropertyString","Purpose","Base").Purpose = """
This ghost object is a temporary object to be used to aid in transforming
the grid.  Transform this object using either the dragger or editing its
Placement property, attach via Part workbench menu Attachment, etc., and 
then in the GridSurface object context menu select Finish transforming grid
to delete the object and apply the transformation to the vertices in the
grid.  You can also click Cancel or simply delete this object to cancel"""        
        
        
        ghost.Shape = comp
        ghost.ViewObject.PointSize = fp.PointSize
        ghost.ViewObject.PointColor = (255,0,0)
        ghost.Document.recompute()
        FreeCADGui.Selection.clearSelection()
        FreeCADGui.Selection.addSelection(fp.Document.Name, ghost.Name)
        FreeCADGui.runCommand("Std_TransformManip",0)
        fp.ViewObject.Visibility = False

     
    def finishTransformingGrid(self,checked):
         """finish the transformation and delete ghost"""
         fp = self.FP
         vert_dict = fp.Proxy.buildDictionary(fp)
         ghost = fp.Document.getObject(self.ghostName)
         if ghost:
            fp.Document.openTransaction("Transform Grid") 
            plm = ghost.Placement
            ii = 0
            for k,v in vert_dict.items():
                vert,row,col = v
                vert.Placement.Base.x, vert.Placement.Base.y, vert.Placement.Base.z = ghost.Shape.Vertexes[ii].CenterOfGravity
                vert.X, vert.Y, vert.Z = (0,0,0)
                ii += 1
            fp.Document.removeObject(self.ghostName)
         fp.ViewObject.Visibility = True
         fp.Document.recompute()
         self.ghostName = None
         fp.Document.commitTransaction()
    
    def cancelTransformingGrid(self, checked):
        """cancel transform and delete ghost"""
        fp = self.FP
        if self.ghostName:
            ghost = fp.Document.getObject(self.ghostName)
            if ghost:
                fp.Document.removeObject(self.ghostName)
        self.ghostName = None
        fp.ViewObject.Visibility = True
        fp.Document.recompute()
    
                
    def addToBody(self, checked):
        fp = self.FP
        body = FreeCADGui.ActiveDocument.ActiveView.getActiveObject("pdbody")
        if not fp in body.Group:
            body.Group += [fp]
            for vert in fp.Vertices:
                if not vert in body.Group:
                    body.Group += [vert]
        
        
    def removeFromBody(self, checked):
        fp = self.FP
        body = FreeCADGui.ActiveDocument.ActiveView.getActiveObject("pdbody")
        if fp in body.Group:
            to_remove = [obj for obj in body.Group if obj == fp or obj in fp.Vertices]
            to_remain = [obj for obj in body.Group if obj not in to_remove]
            body.Group = to_remain

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

    def switchMethod(self, mode):
        fp = self.FP
        fp.Method = mode
        fp.recompute()

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
        
    def dropObject(self, vp, dropped):
        fp = self.FP
        if not dropped in fp.Vertices:
            fp.Vertices = fp.Vertices + [dropped]

    def canDropObjects(self):
        return True

    def canDropObject(self, dropped):
        return dropped.isDerivedFrom("Part::Vertex") or \
                bool(dropped.isDerivedFrom("Part::Feature") and \
                len(dropped.Shape.Vertexes) == 1 and \
                hasattr(dropped,"X") and hasattr(dropped,"Y") and hasattr(dropped,"Z"))


    def __getstate__(self):
        '''When saving the document this object gets stored using Python's json module.\
            If we have any unserializable stuff return them here or None'''
        return None

    def __setstate__(self,state):
        '''When restoring the serialized object from document we have the chance to set some internals here.\
                Since no data were serialized nothing needs to be done here.'''
        return None

class MeshRemodelCreateGridSurfaceCommandClass:
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
        sel = FreeCADGui.Selection.getCompleteSelection()
        picked = global_picked if global_picked else []
        if not picked:
            for s in sel:
                picked.extend(s.PickedPoints)

        fp = doc.addObject("Part::FeaturePython", "GridSurface")
        GridSurface(fp)
        GridSurfaceVP(fp.ViewObject)
        doc.recompute()
        
        if len(picked) > 1:
            fp.CountRows = 1
            fp.CountColumns = len(picked)
            fp.recompute()
            for idx, vert in enumerate(fp.Vertices):
                vert.Placement.Base.x, vert.Placement.Base.y, vert.Placement.Base.z = picked[idx].x, picked[idx].y, picked[idx].z
            doc.recompute()
            
        elif sel:
            if len(sel) == 1 and sel[0].Object.isDerivedFrom("Part::Feature") and sel[0].Object.Shape.ShapeType == "Wire":
                fp.CountRows = 1
                fp.CountColumns = len(sel[0].Object.Shape.Vertexes)
                fp.recompute()
                for idx, vert in enumerate(fp.Vertices):
                    vert.Placement.Base.x, vert.Placement.Base.y, vert.Placement.Base.z = sel[0].Object.Shape.Vertexes[idx].Point
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
            #show the center point on ctrl click or shift+ctrl click
            vert = doc.addObject("Part::Vertex", circName+"_Ctr")
            vert.Placement.Base.x, vert.Placement.Base.y, vert.Placement.Base.z = center
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
            #Part.show(Part.Point(center).toShape(),arcName+"_Ctr") #show the center point
            pv = doc.addObject("Part::Vertex",arcName+"_Ctr")
            pv.Placement.Base.x, pv.Placement.Base.y, pv.Placement.Base.z = center
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

class SketchPlus:

# class level attributes
# wire filter dialogs use this dictionary of shapes to create the ghost objects and access the wires
# each function sets its key or keys to the incoming shape for filtering in the dialogs
    wireShapeDict = {}
    ScaleUniformCenterDefault = ["Origin"]
    MirrorAxisDefault = ["None","Sketch X Axis","Sketch Y Axis","Both X and Y"]
    PathArraySettingsDefault = ["Normal To Path", "Align To Origin", "Fixed", "Tangent To Path", "Normal To Path 3D"]
    PointArraySettingsDefault = ["No Point Array","Align To Origin","Fixed"]
    PolarCenterDefault = ["Origin"]
    RotateSettingsDefault = ["No Rotation", "Center Of Gravity", "Origin"]
    ExecutionOrderDefault = ["Wire Order", "External", "Construction", "Offsetting", "Mirroring","Polar Array", \
        "Point Array", "Rectangular Array", "Rotate", "Path Array", "Scaling","Uniform Scaling", "Face Making", \
        "Wire Statuses"  ]
    RectangularRotationCenterDefault = ["Origin","Center Of Gravity"]

    def __init__(self, obj):
        obj.Proxy = self

        generalWireInstructions = \
"""
Wire filtering uses 1-based indices for the wires, example 1 = Wire1,  2 = Wire2,
etc.  An empty list [] signifies to use all the wires.  Generally, any index left
out gets used and any negative values get filtered out.

Suppose you have 3 wires: Wire1, Wire2, and Wire3.

[] = Use all 3 wires
[1] = Use Wire1, Wire2, and Wire3
[-1, 2, 3] = Use Wire2 and Wire3, but not Wire1
[-1. -2, -3] = Use none of the wires
[-1] = Use Wire2 and Wire3, but not Wire1
[3, 1, 2] = reorder the wires so Wire3 now becomes Wire1, Wire1 becomes Wire2, and
            Wire2 becomes Wire3.

Note that wires available will be different depending on the order of function
execution, which can be set with the Execution Order property.

"""
        grp = "Additional Plus Properties"

        ##### Mirroring

        obj.addProperty("App::PropertyEnumeration","MirrorAxis",grp,\
"""Axis for mirroring sketch elements.  Options are X, Y, X and Y, or Construction
Line Segment.  If you use X and Y and your elements are constrained to the X and
Y axes, it can produce a closed wire, but sometimes this fails.  I have found
that by removing and adding a constraint back can sometimes correct the issue.
""").MirrorAxis = SketchPlus.MirrorAxisDefault
        obj.addProperty("App::PropertyIntegerList","MirrorWires",grp,
f"""Wires to mirror.
{generalWireInstructions}

Unused wires remain stationary.  To remove a wire entirely, use the Wire Order property.

""")
        obj.addProperty("App::PropertyString","Version",grp,"Version of MeshRemodel used to create this SketchPlus object").Version = __version__

        ##### Scaling
        obj.addProperty("App::PropertyFloat","ScaleX",grp, "Scale in the X direction, converts to bspline curves if scaling is non-uniform").ScaleX = 1
        obj.addProperty("App::PropertyFloat","ScaleY",grp,"Scale in the Y direction, converts to bspline curves if scaling is non-uniform").ScaleY = 1
        obj.addProperty("App::PropertyIntegerList","ScaleWires",grp,"Wires to scale")


        ##### Scale Uniform
        obj.addProperty("App::PropertyFloat","ScaleUniform",grp,"Scale X and Y uniformly, does not convert to bsplines").ScaleUniform = 1

        obj.addProperty("App::PropertyEnumeration","ScaleUniformCenter",grp,\
"""Default is origin, but can be construction circle, construction point or
construction arc of circle.  Uniform scaling has the advantage that it does not
produce bsplines the way the other scaling often does (if non-uniform, x != y).
It also can make use of an arbitrary center where the other scaling feature
does not provide that option.""").ScaleUniformCenter = SketchPlus.ScaleUniformCenterDefault
        obj.addProperty("App::PropertyIntegerList","ScaleUniformWires",grp,"Wires to scale uniformly")


        ##### Offsetting
        obj.addProperty("App::PropertyDistance","Offset",grp,"Offset to be applied to the sketch").Offset = 0
        obj.addProperty("App::PropertyEnumeration","OffsetJoin",grp,"Offset join types").OffsetJoin = ["Arcs","Tangent", "Intersection"]
        obj.addProperty("App::PropertyBool","OffsetFill",grp,"Whether to fill in gap between original and offset").OffsetFill = False
        obj.addProperty("App::PropertyBool","OffsetIntersection",grp,"Affects how child edges and wires are handled, whether as a group or independently").OffsetIntersection = False
        obj.addProperty("App::PropertyBool","OffsetOpenResult",grp,"Whether to allow offset to produce open wire").OffsetOpenResult = False
        obj.addProperty("App::PropertyIntegerList","OffsetWires",grp,"Wires to offset")
        obj.addProperty("App::PropertyBool","OffsetInterimFace",grp,"Make a temporary face before offsetting, then use the generated wires afterwards.  Seems to make the offsetting more robust.").OffsetInterimFace = True


        obj.addProperty("App::PropertyEnumeration","FaceMaker",grp,\
"""Whether to make a face, and which face maker to use.  Bullseye excels at faces
within faces within faces and supports much nesting, Cheese excels at multiple
faces but only 1 layer deep, Simple only uses the outer wire.  In the actual
implementation, Simple uses Bullseye to make the face first, then another face
is made from the OuterWire of the first face, which is the return value.  If
FaceMakerSimple is used with inner wires it gives self-intersections, so that is
why Bullseye is used internally in this script.""").FaceMaker = \
            ["None","Bullseye","Cheese","Simple"]
        obj.addProperty("App::PropertyBool","ShowConstruction", grp, \
"""Whether to show construction geometry as normal geometry.  This is intended
as a way to quickly view construction geometry without needing to toggle it
to normal and back again.  For permanent use it is generally better to toggle the
construction geometry to normal mode in the sketch editor.  But this can also
be a way to dynamically change sketch elements.""").ShowConstruction = False
        obj.addProperty("App::PropertyBool","ShowExternal",grp,\
"""Whether to show external geometry as normal geometry.  This can save some
time by not requiring you to trace over the imported external geometry in cases
where you want it to be normal geometry.  Bsplines can be more of a challenge to
recreate from tracing, and this should recreate them flawlessly.  It is important
to understand that the geometry will change quite a lot as the external links
change or as the sketch is moved or rotated, so the results can sometimes be
unexpected.""").ShowExternal = False
        obj.addProperty("App::PropertyIntegerList","WireOrder",grp,\
f"""The order of the Wires in the sketch.
{generalWireInstructions}
""")

        ##### Path Array
        obj.addProperty("App::PropertyIntegerList","PathArrayWires",grp,"""
Wires to be used in the PathArray as copies""")
        obj.addProperty("App::PropertyIntegerList","PathArrayPath", grp, """
The Wire to be used as the spline for the Path Array, should only be one wire usually.""")
        obj.addProperty("App::PropertyIntegerConstraint","PathArrayCount",grp,"""
Number of copies of the path array wires to made in the path array, number of elements.""").PathArrayCount = (0,0,100000,1)
        obj.addProperty("App::PropertyBool","PathArrayKeepPath",grp,"Whether to keep the path wire after making the path array").PathArrayKeepPath = False

        obj.addProperty("App::PropertyEnumeration","PathArraySettings",grp,"Some settings for the path array").PathArraySettings = SketchPlus.PathArraySettingsDefault

        ##### Point Array
        obj.addProperty("App::PropertyIntegerList","PointArrayWires",grp,
f"""List the Wires to use for the base of the point array here.
{generalWireInstructions}
""")

        obj.addProperty("App::PropertyEnumeration","PointArraySettings",grp,
"""If "Fixed", then no rotation is done for each copy.  If "Align To Origin", then
the rotation is based on the point's angle relative to the x axis.  If "No Point
Array" then there will be no point arraying done.
""").PointArraySettings = SketchPlus.PointArraySettingsDefault
        ##### Polar Array
        obj.addProperty("App::PropertyIntegerConstraint", "PolarCount",grp, "Numer of elements in the polar array").PolarCount = (0,0,10000,1)
        obj.addProperty("App::PropertyAngle","PolarStart",grp,"Start of Polar array in degrees, default = 0 degrees").PolarStart = 0
        obj.addProperty("App::PropertyAngle","PolarEnd",grp,"End of polar array in degrees, default = 360 degrees").PolarEnd = 360
        obj.addProperty("App::PropertyFloatList", "PolarAngles",grp,
"""Optional list of custom angles in degrees.  If this is set, then PolarStart,
PolarEnd and PolarCount properties are ignored.  Include 0 if you also want the
original to be included in the array.

This can be used to rotate a circle or other single wire in place.  Just add a
construction point coincident to the center of the wire and use
that point as the PolarCenter and set a single angle in this PolarAngles property.
""").PolarAngles = []

        obj.addProperty("App::PropertyEnumeration","PolarCenter",grp,"Center of the polar array").PolarCenter = SketchPlus.PolarCenterDefault
        obj.addProperty("App::PropertyIntegerList","PolarWires",grp,
f"""Wires to use as polar array pattern sources.
{generalWireInstructions}
Note that all rogue points (non-construction mode points) are also made part of the
polar array.  There is currently no way to exclude them other than not having
rogue points in your sketch.  Rogue points can come to be in a number of ways,
including directly created non-construction mode points, construction points
converted to normal points when Show Construction = True, and external points
that are converted to normal points when ShowExternal = True.
""").PolarWires = []


        obj.addProperty("App::PropertyStringList","ExecutionOrder",grp,\
"""The order of execution of the various functions, but also a way to suppress functions.

There is a dialog in the context menu to edit this string list or you may do it manually.
Put a hyphen (-) immediately in front of the name, e.g. "-Wire Order" to disable that function.
Do not use any spaces in between the hyphen and the start of the name, e.g. "- Wire Order" is
invalid and will result in an exception being thrown and execution halted.

If a function is included twice (or more) in the list it will be executed twice, but this was
not the intention of this feature and you might not get expected results if you do this.  If there
is a property associated with the function, then that same property gets used both times.

The order of execution can be important in some cases.  Note that the wire filter
editors in the context menu will only show the available number of wires for a function
at the time the function is being executed.  So, as an example, if you have a single
wire to begin with and use a polar array to make 8 wires out of it, then the functions
that follow the polar array will have 8 wires to use.  If a mirror operation follows
the polar array, then the MirrorWires property has 8 wires to filter where if the
Mirror function comes before the polar array then MirrorWires cam only filter the
first oriiginal wire.

""").ExecutionOrder = SketchPlus.ExecutionOrderDefault


        ##### Rectangular Array
        obj.addProperty("App::PropertyIntegerConstraint","RectangularCountX",grp,"Count of elements in x direction").RectangularCountX = (1,1,100000,1)
        obj.addProperty("App::PropertyIntegerConstraint","RectangularCountY",grp,"Count of elements in y direction").RectangularCountY = (1,1,100000,1)
        obj.addProperty("App::PropertyDistance","RectangularIntervalX",grp,"Distance between elements in x direction").RectangularIntervalX = 10
        obj.addProperty("App::PropertyDistance","RectangularIntervalY",grp,"Distance between elements in y direction").RectangularIntervalY = 10
        obj.addProperty("App::PropertyIntegerList","RectangularWires",grp,"Wires to use in rectangular array")
        obj.addProperty("App::PropertyAngle","RectangularRotation",grp,"Angle by which to rotate the rectangular array as a group").RectangularRotation = 0

        obj.addProperty("App::PropertyEnumeration","RectangularRotationCenter",grp,"Center of rotation for rectangular array, also supported: construction circles, arcs, and points").RectangularRotationCenter = SketchPlus.RectangularRotationCenterDefault
        
        ##### Rotate wires
        obj.addProperty("App::PropertyIntegerList","RotateWires", grp, "Wires to rotate")
        obj.addProperty("App::PropertyFloatList","RotateAngles",grp,"Angles to rotate wires")
        obj.addProperty("App::PropertyEnumeration","RotateSettings",grp,"Rotate settings").RotateSettings = SketchPlus.RotateSettingsDefault

        grp = "Additional Status Information"
        obj.addProperty("App::PropertyString", "WireStatuses", grp, "Information about each wire")

    def onDocumentRestored(self, fp):
        
        #update legacy objects with new properties, if necessary:
        if not hasattr(fp,"RotateWires"):
            print("updating legacy object with new rotate properties")
            grp = "Additional Plus Properties"
            fp.addProperty("App::PropertyIntegerList","RotateWires", grp, "Wires to rotate")
            fp.addProperty("App::PropertyFloatList","RotateAngles",grp,"Angles to rotate wires")
            fp.addProperty("App::PropertyEnumeration","RotateSettings",grp,"Rotate settings").RotateSettings = SketchPlus.RotateSettingsDefault
            SketchPlus.wireShapeDict["RotateWires"] = Part.Shape()


        # this updates all the keys in the SketchPlus.wireShapeDict dictionary
        # done in a singleshot because otherwise there is an access violation if 
        # there are links to external geometry in the object
        QtCore.QTimer().singleShot(100, lambda: fp.recompute)

    def onChanged(self, fp, prop):
        if prop == "OffsetFill":
            if fp.OffsetFill:
                fp.ViewObject.DisplayMode = "Flat Lines" # do it here so if the user toggles back to Wireframe it will stay there
            else:
                if "Wireframe" in fp.ViewObject.getEnumerationsOfProperty("DisplayMode"):
                    fp.ViewObject.DisplayMode = "Wireframe" # facemaker function toggles it back if necessary

    def fetchConstruction(self, fp, prop, defaults, filter=["LineSegment"]):
        """get construction elements as a dictionary and put them into prop
        prop is typically an enumeration of options for centering something,
        such as a polar array or for a line segment to serve as an axis for mirroring
        """
        current = getattr(fp, prop)
        filter = [f"Part::Geom{val}" for val in filter]

        setattr(fp, prop, defaults)
        constructionElements = {}
        for idx,geo in enumerate(fp.Geometry):
                if fp.getConstruction(idx):
                    if geo.TypeId in filter:
                        elementName = f"Construction {geo.TypeId[10:]} {len(constructionElements) + 1}"
                        constructionElements[elementName] = geo
                        enums = fp.getEnumerationsOfProperty(prop)
                        if not elementName in enums:
                            setattr(fp, prop, enums + [elementName])
                    # else:
                    #     print(f"skipping {geo.TypeId}, not in {filter}")
        enums = fp.getEnumerationsOfProperty(prop)
        if current in enums:
            setattr(fp, prop, current)
        else:
            FreeCAD.Console.PrintWarning(f"{fp.Label}: {current} no longer exists, setting {prop} to defaults\n")
            setattr(fp, prop, defaults)
        return constructionElements

    def handleMirroring(self, fp, shape):
        """handle mirroring here instead of execute()"""
        constructionLines = self.fetchConstruction(fp, "MirrorAxis", SketchPlus.MirrorAxisDefault, ["LineSegment"])
        SketchPlus.wireShapeDict["MirrorWires"] = shape
        if shape.isNull():
            return shape

        normal2 = None # will become (0,1,0) if we are mirroring both x and y


        if fp.MirrorAxis == "None":
            return shape
        else:
            #shape will be transformed back to the identity placement before mirroring
            base = FreeCAD.Vector(0,0,0)
            if fp.MirrorAxis == "Sketch Y Axis":
                normal = FreeCAD.Vector(1, 0, 0)
            elif fp.MirrorAxis == "Sketch X Axis":
                normal = FreeCAD.Vector(0, 1, 0)
            elif fp.MirrorAxis == "Both X and Y":
                normal = FreeCAD.Vector(1,0,0) # do the Y axis first
                normal2 = FreeCAD.Vector(0,1,0) # then the X
            elif fp.MirrorAxis in constructionLines:
                geo = constructionLines[fp.MirrorAxis] #construction linesegment
                direction = (geo.StartPoint - geo.EndPoint).normalize()
                base = geo.StartPoint #any point will do, line treated as infinite
                normal = FreeCAD.Vector(direction.y, -direction.x, 0)

            patternWires, stationaryWires = self.filterWires(shape, fp.MirrorWires)
            if not patternWires and not fp.MirrorWires:
                patternWires = shape.Wires #use all wires if user does not filter them
                stationaryWires = []
            roguePoints = self.getRoguePoints(fp, shape)

            if not patternWires and not roguePoints:
                FreeCAD.Console.PrintError("Nothing to mirror, add some Wires to the Mirror Wires property, add some non-construction mode Points or set Mirror Axis to 'None' to avoid this message.\n")
                return shape
            patternShape = Part.makeCompound(patternWires + roguePoints)
            mirror = patternShape.mirror(base, normal)
            fusion = Part.makeCompound([patternShape,mirror])

            if normal2:
                mirror2 = fusion.mirror(base, normal2)
                fusion = Part.makeCompound([fusion, mirror2])

            # this helps connect the wires that meet at the axes, but sometimes still fails
            # Note that without this we just get a collection of edges, not wired together

            wires = []
            for edgelist in Part.sortEdges(fusion.Edges):
                wires.append(Part.Wire(edgelist))

            comp = Part.makeCompound(wires + stationaryWires + self.getRoguePoints(fp, fusion))
            return comp

    def handleScaling(self, fp, shape):
        """handle scaling, if applicable"""
        SketchPlus.wireShapeDict["ScaleWires"] = shape
        if fp.ScaleX == 1 and fp.ScaleY == 1:
            return shape
        if shape.isNull():
            return shape
        matrix = FreeCAD.Matrix()
        matrix.A11 = fp.ScaleX
        matrix.A22 = fp.ScaleY
        matrix.A33 = 1 if fp.ScaleX != fp.ScaleY else fp.ScaleX #might not make bsplines if x and y are equal

        used, unused = self.filterWires(shape, fp.ScaleWires)
        roguePoints = self.getRoguePoints(fp, shape)
        if used + roguePoints:
            usedComp = Part.makeCompound(used + roguePoints)
        else:
            return shape # there was nothing to scale

        scaledShape = usedComp.transformShape(matrix, True, True) #copy, checkScale
        if not unused:
            return usedComp
        else:
            unusedComp = Part.makeCompound(unused)
            fusion = usedComp.fuse(unusedComp)
            return fusion


    def handleUniformScaling(self, fp, shape):
        """handle uniform scaling"""
        SketchPlus.wireShapeDict["ScaleUniformWires"] = shape
        # find construction circles and add them to the enumeration
        constructionCircles = self.fetchConstruction(fp, "ScaleUniformCenter", SketchPlus.ScaleUniformCenterDefault, filter=["Circle","Point","ArcOfCircle"])
        if fp.ScaleUniform == 1.0:
            return shape
        if shape.isNull():
            return shape

        center = FreeCAD.Vector()

        if fp.ScaleUniformCenter != "Origin":
            geo = constructionCircles[fp.ScaleUniformCenter]
            center = geo.Location if hasattr(geo, "Location") else FreeCAD.Vector(geo.X, geo.Y, geo.Z)

        used, unused = self.filterWires(shape, fp.ScaleUniformWires)
        roguePoints = self.getRoguePoints(fp, shape)
        if used + roguePoints:
            usedComp = Part.makeCompound(used + roguePoints)
        else:
            return shape # there was nothing to scale
        usedComp.scale(fp.ScaleUniform, center)
        if not unused:
            return usedComp
        else:
            unusedComp = Part.makeCompound(unused)
            fusion = usedComp.fuse(unusedComp)
            return fusion


    def handleOffset(self, fp, shape):
        """Offset the shape and return the offset if Offset != 0"""
        SketchPlus.wireShapeDict["OffsetWires"] = shape
        if shape.isNull():
            return shape
        offset = fp.Offset.Value
        if offset == 0:
            return shape
        join_dict = {"Arcs": 0, "Tangent":1, "Intersection":2}
        join = join_dict[fp.OffsetJoin]
        fill = fp.OffsetFill
        openResult = fp.OffsetOpenResult
        intersection = fp.OffsetIntersection

        patternWires, stationaryWires = self.filterWires(shape, fp.OffsetWires)
        roguePoints = self.getRoguePoints(fp, shape)
        #not offsetting rogue points
        if not patternWires:
            return shape #nothing to offset
        else:
            patternWiresComp = Part.makeCompound(patternWires)

            try:
                if fp.OffsetInterimFace:
                    try:
                        compOrFace = Part.makeFace(patternWiresComp,"Part::FaceMakerBullseye")
                    except:
                        FreeCAD.Console.PrintMessage(f"{fp.Label}: unable to make interim face, proceeding without it.\n")
                        compOrFace = patternWiresComp
                else:
                    compOrFace = patternWiresComp
                offsetShape = compOrFace.makeOffset2D(offset, join = join, \
                            fill = fill, openResult = openResult, \
                            intersection = intersection)
                if not fill:
                    offsetShape = Part.makeCompound(offsetShape.Wires)

            except Exception as e:
                FreeCAD.Console.PrintError(f"{fp.Label}: Error offsetting: {e}\n")
                return shape
            if offsetShape.isNull():
                FreeCAD.Console.PrintError("Error in offsetting, output shape is null")
                return shape
            if stationaryWires + roguePoints:
                final_shape = Part.makeCompound(stationaryWires + roguePoints + [offsetShape])
            else:
                final_shape = offsetShape
            return final_shape


        return shape

    def handleFaceMaker(self, fp, shape):
        """Make a face using selected facemaker"""
        if shape.isNull():
            return shape
        if fp.FaceMaker == "None":
            return shape
        fp.ViewObject.DisplayMode = "Flat Lines"
        facemaker = f"Part::FaceMaker{fp.FaceMaker}" if fp.FaceMaker != "Simple" else "Part::FaceMakerBullseye"
        try:
            face = Part.makeFace(shape, facemaker)
            if fp.FaceMaker == "Simple":
                face = Part.makeFace(face.OuterWire, facemaker)
            return face
        except Exception as e:
            FreeCAD.Console.PrintError(f"{fp.Label}: error making face: {e}\n")
            return shape

    def getRoguePoints(self, fp, shape):
        """get vertcies that are not part of any edges"""
        def vertInList(vertex, vertLists):
            for vertList in vertLists:
                for vert in vertList:
                    if vertex.Point.isEqual(vert.Point, Part.Precision.confusion()):
                        return True
            return False
        if hasattr(Part.Wire, "getChildshapes"):
            children = [wire.getChildShapes(f"Vertex{idx+1}") for idx,wire in enumerate(shape.Wires)]
        else:
            children = []
            for idx,wire in enumerate(shape.Wires):
                inner_children = []
                for v in wire.Vertexes:
                    inner_children.append(v)
                children.append(inner_children)
        rogues = [v for v in shape.Vertexes if not vertInList(v, children)]

        return rogues

    def handleShowConstruction(self, fp, shape):
        if not fp.ShowConstruction:
            return shape
        constructionGeometry = []
        constructionPoints = []
        wires = []
        roguePoints = self.getRoguePoints(fp, shape) #non construction mode points pre-existing

        for idx,geo in enumerate(fp.Geometry):
            if fp.getConstruction(idx):
                gshape = geo.toShape()
                if "Point" in geo.TypeId:
                    constructionPoints.append(gshape)
                else:
                    constructionGeometry.append(gshape)

        if constructionGeometry:
            comp = Part.makeCompound(constructionGeometry)
        else:
            comp = Part.Shape()
        edges = shape.Edges + comp.Edges
        sortedEdges = Part.sortEdges(edges)
        wires = [Part.Wire(e) for e in sortedEdges]


        if constructionPoints or constructionGeometry or roguePoints:
            combined = wires + constructionPoints + roguePoints
            wireComp = Part.makeCompound(combined)
            return wireComp
        else:
            return shape

    def handleWireOrder(self, fp, shape):
        """handle the ordering of the wires"""
        SketchPlus.wireShapeDict["WireOrder"] = shape
        if shape.isNull():
            return shape
        used, disregarded = self.filterWires(shape, fp.WireOrder)
        roguePoints = self.getRoguePoints(fp, shape)
        if used + roguePoints:
            return Part.makeCompound(used + roguePoints)
        else:
            return Part.Shape()

    def filterWires(self, shape, order, addUnmentioned = True):
        """returns (filteredIn, filteredOut) both lists of wires.  +1 = put Wire1
        into the filteredIn list, -1 = put it into the filteredOut list.  If
        addUnmentioned = True, then unmentioned wires get added to the filteredIn
        list.
        """
        #user might have deleted some wires from sketch, so we make sure
        #to filter ot any unavailable wire indices
        available = [idx + 1 for idx in range(len(shape.Wires))]
        available += [-d for d in available]

        outs = [-wo for wo in order if wo < 0 and wo in available]
        ins = [wo for wo in order if wo > 0 and wo in available]
        if 0 in order:
            return ([], shape.Wires) #all wires unused
        filteredOut = [shape.Wires[wo-1] for wo in outs]
        filteredIn = [shape.Wires[wo-1] for wo in ins]
        if addUnmentioned:
            filteredIn += [shape.Wires[wo] for wo in range(len(shape.Wires)) if wo+1 not in outs + ins]
            ins += [wo+1 for wo in range(len(shape.Wires)) if wo+1 not in outs + ins]

        return (filteredIn, filteredOut)


    def handleExternal(self, fp, shape):
        """ create a compound from the external geometry by projecting it to the sketch plane
            the compound is combined with input shape and a compound of both is returned
        """
        if not fp.ShowExternal:
            return shape
            
        if not hasattr(fp, "ExternalGeo"):
            FreeCAD.Console.PrintError(f"{fp.Label}: Show External geometry is unavailable in this version of FreeCAD.  You must update to a newer version to enjoy this feature.  Resetting ShowExternal property to False to prevent further display of this error message.\n")
            fp.ShowExternal = False
            return shape
            
        externalGeometry = []
        externalPoints = []
        wires = []

        for idx,geo in enumerate(fp.ExternalGeo[2:]):
            if not "Point" in geo.TypeId:
                externalGeometry.append(geo.toShape())
            else:
                externalPoints.append(geo.toShape())

        if externalGeometry:
            comp = Part.makeCompound(externalGeometry)
        else:
            comp = Part.Shape()

        roguePoints = self.getRoguePoints(fp, shape)
        edges = shape.Edges + comp.Edges
        sortedEdges = Part.sortEdges(edges)
        wires = [Part.Wire(e) for e in sortedEdges]

        if externalGeometry or externalPoints or roguePoints:
            wireComp = Part.makeCompound(wires + externalPoints + roguePoints)
            return wireComp

        return shape
        
    #rotate wires
    def handleRotate(self, fp, shape):
        """Rotate one or more wires """
        SketchPlus.wireShapeDict["RotateWires"] = shape
        constructionCircles = self.fetchConstruction(fp, "RotateSettings", SketchPlus.RotateSettingsDefault, filter=["Circle","Point","ArcOfCircle"])            
        if shape.isNull():
            return shape

        if fp.RotateSettings == "No Rotation":
            return shape
        
        if not fp.RotateAngles:
            return shape
            
        roguePoints = self.getRoguePoints(fp, shape)

        patternWires, stationaryWires = self.filterWires(shape, fp.RotateWires)
        
        if not fp.RotateWires:
            patternWires = shape.Wires #use all wires as default if RotateWires is empty
            stationaryWires = []

        if not patternWires:
            return shape

        angles = fp.RotateAngles
        #use the last angle given for the remaining wires, if necessary
        while len(angles) < len(patternWires):
            angles.append(angles[-1])

        center = FreeCAD.Vector() #for default Origin
        if fp.RotateSettings not in SketchPlus.RotateSettingsDefault:
            geo = constructionCircles[fp.RotateSettings]
            center = geo.Location if hasattr(geo, "Location") else FreeCAD.Vector(geo.X, geo.Y, geo.Z)
        elif fp.RotateSettings == "Origin":
            center = FreeCAD.Vector()

        copies = []
        for idx, (angle, wire) in enumerate(zip(angles, patternWires)):
            copy = wire.copy()
            if fp.RotateSettings == "Center Of Gravity":
                center = wire.CenterOfGravity
            copy.rotate(center, FreeCAD.Vector(0,0,1), angle)
            copies.append(copy)

        copiesAndStationaryComp = Part.makeCompound(copies + stationaryWires)
        wires = []
        edges = copiesAndStationaryComp.Edges
        if edges:
            se = Part.sortEdges(edges)
            wires = [Part.Wire(s) for s in se]
        finalShape = Part.makeCompound(wires + roguePoints)
        return finalShape
        #end handleRotate()

    def handlePolar(self, fp, shape):
        """Polar array """
        SketchPlus.wireShapeDict["PolarWires"] = shape
        constructionCircles = self.fetchConstruction(fp, "PolarCenter", SketchPlus.PolarCenterDefault, filter=["Circle","Point","ArcOfCircle"])
        if shape.isNull():
            return shape
        count = len(fp.PolarAngles) if fp.PolarAngles else fp.PolarCount
        roguePoints = self.getRoguePoints(fp, shape)
        if not count:
            return shape
        patternWires, stationaryWires = self.filterWires(shape, fp.PolarWires)
        if not fp.PolarWires:
            patternWires = shape.Wires #use all wires as default if PolarWires is empty
            stationaryWires = []
        if patternWires + roguePoints:
            patternComp = Part.makeCompound(patternWires + roguePoints)
        else:
            FreeCAD.Console.PrintError("Nothing to pattern, skipping Polar Pattern function\n")
            return shape

        center = FreeCAD.Vector() #for default Origin
        if fp.PolarCenter != "Origin":
            geo = constructionCircles[fp.PolarCenter]
            center = geo.Location if hasattr(geo, "Location") else FreeCAD.Vector(geo.X, geo.Y, geo.Z)

        start = fp.PolarStart.Value
        end = fp.PolarEnd.Value
        rng = end - start
        step = rng / count
        if fp.PolarAngles:
            angles = fp.PolarAngles
        else:
            angles = [step * occurrence for occurrence in range(count)]
       # print(f"angles = {angles} from start,end,rng,step,count: {start},{end},{rng},{step},{fp.PolarCount}")
        copies = []
        for angle in angles:
            copy = patternComp.copy()
            copy.rotate(center, FreeCAD.Vector(0,0,1), angle)
            copies.append(copy)

        copiesAndStationaryComp = Part.makeCompound(copies + stationaryWires)
        roguePoints2 = self.getRoguePoints(fp, copiesAndStationaryComp)

        wires = []
        edges = copiesAndStationaryComp.Edges
        if edges:
            se = Part.sortEdges(edges)
            wires = [Part.Wire(s) for s in se]
        finalShape = Part.makeCompound(wires + roguePoints2)
        return finalShape

    def wireInWireList(self, wire, wireList):
        """check if wire is in wireList"""
        for w in wireList:
            if w.isEqual(wire):
                return True
        return False


    def handlePathArray(self, fp, shape):
        """makes a path array"""
        SketchPlus.wireShapeDict["PathArrayWires"] = shape
        SketchPlus.wireShapeDict["PathArrayPath"] = shape
        if shape.isNull():
            return shape

        if fp.PathArrayCount == 0:
            return shape

        if len(shape.Wires) == 0:
            FreeCAD.Console.PrintError(f"{fp.Label}: no wires to use in path array, skipping function.  Note that points are not patterned in this function.\n")
        if len(fp.PathArrayWires) > len(shape.Wires):
            FreeCAD.Console.PrintError("PathArrayWires property must have less than\
or equal to the available number of wires.\n")
            return shape

        roguePoints = self.getRoguePoints(fp,shape)
        patternWires, stationaryWires = self.filterWires(shape, fp.PathArrayWires)
        paths, unused = self.filterWires(shape, fp.PathArrayPath)
        if not fp.PathArrayKeepPath:
            stationaryWires = [w for w in stationaryWires if not self.wireInWireList(w, paths)]

        if not patternWires + roguePoints:
            FreeCAD.Console.PrintError("No pattern wires to use in path array, skipping\n")
            return shape

        if not paths:
            FreeCAD.Console.PrintError("No path selected to use for path array, skipping\n")
            return shape

        discretePoints = []
        for path in paths:
            extra = 1 if path.isClosed() else 0
            discretePoints.append(path.discretize(fp.PathArrayCount + extra))

        #not arraying roguePoints for this function
        copies = []
        for idx, pathPts in enumerate(discretePoints):
            path = paths[idx]  # path is a Part.Wire object
            edges = path.Edges  # Get the edges of the wire

            for i, vec in enumerate(pathPts[extra:]):
                for wire in patternWires + roguePoints:
                    copy = wire.copy()
                    vector = vec - copy.CenterOfGravity
                    # since we discretized the wire and not the edges, we have to figure
                    # out which edge contains this point
                    containing_edge = None
                    for edge in edges:
                        (dist, vectors, infos) = edge.distToShape(Part.Vertex(vec))
                        if math.fabs(dist) < Part.Precision.confusion():
                            containing_edge = edge
                            if infos[0][0] == "Edge":
                                u = infos[0][2]
                            else: #infos[0][0] == "Vertex":
                                u = infos[0][1] # 1.0 or 0.0 if this is a vertex on the edge
                            break

                    if containing_edge is None:
                        print(f"skipping vector {vec}, no containing edge found")
                        continue  # should not happen

                    if fp.PathArraySettings == "Align To Origin":
                        angle = math.degrees(math.atan2(vec.y, vec.x))
                        copy.rotate(copy.CenterOfGravity, FreeCAD.Vector(0, 0, 1), angle)

                    elif fp.PathArraySettings == "Normal To Path":
                        tangent = containing_edge.tangentAt(u)
                        normal = FreeCAD.Vector(-tangent.y, tangent.x, 0)  #rotate 90 degrees
                        angle = math.degrees(math.atan2(normal.y, normal.x))
                        copy.rotate(copy.CenterOfGravity, FreeCAD.Vector(0, 0, 1), angle)
                        
                    elif fp.PathArraySettings == "Normal To Path 3D":
                        tangent = containing_edge.tangentAt(u)
                        normal = FreeCAD.Vector(-tangent.y, tangent.x, 0)  #rotate 90 degrees
                        angle = math.degrees(math.atan2(normal.y, normal.x))
                        copy.rotate(copy.CenterOfGravity, FreeCAD.Vector(0,1,0),90)
                        copy.rotate(copy.CenterOfGravity, FreeCAD.Vector(0, 0, 1), angle)
                        copy.rotate(copy.CenterOfGravity, FreeCAD.Vector(0,0,1),90)



                    elif fp.PathArraySettings == "Tangent To Path":
                        tangent = containing_edge.tangentAt(u)
                        angle = math.degrees(math.atan2(tangent.y, tangent.x))
                        copy.rotate(copy.CenterOfGravity, FreeCAD.Vector(0, 0, 1), angle)
                    elif fp.PathArraySettings == "Fixed":
                        pass # no rotation

                    copy.Placement.translate(vector)
                    copies.append(copy)

        copiesAndStationaryComp = Part.makeCompound(copies + stationaryWires)
        roguePoints2 = self.getRoguePoints(fp, copiesAndStationaryComp)

        wires = []
        edges = copiesAndStationaryComp.Edges
        if edges:
            se = Part.sortEdges(edges)
            wires = [Part.Wire(s) for s in se]
        finalShape = Part.makeCompound(wires + roguePoints2)
        return finalShape

    def handlePointArray(self, fp, shape):
        """makes a point array of desired wires at construction point locations"""
        SketchPlus.wireShapeDict["PointArrayWires"] = shape
        if shape.isNull():
            return shape
        if fp.PointArraySettings == "No Point Array":
            return shape
        if len(shape.Wires) == 0:
            FreeCAD.Console.PrintError(f"{fp.Label}: no wires to use in point array, skipping function.  Note that points are not patterned in this function.\n")
        if len(fp.PointArrayWires) > len(shape.Wires):
            FreeCAD.Console.PrintError("PointArrayWires property must less than\
or equal to the available number of wires.\n")
            return shape

        patternWires, stationaryWires = self.filterWires(shape, fp.PointArrayWires)
        if not patternWires and not fp.PointArrayWires:
            patternWires = shape.Wires
            stationaryWires = []

        constructionPoints = []
        for geoid, geo in enumerate(fp.Geometry):
            if fp.getConstruction(geoid):
                if "Point" in geo.TypeId:
                    constructionPoints.append(geo.toShape())
        if not constructionPoints:
            FreeCAD.Console.PrintError(f"No construction mode points found in {fp.Label}.  The\
PointArray function requires there be at least one construction mode point in the sketch.  \
Skipping function\n")
            return shape


        #not arraying roguePoints for this function
        copies = []
        patternWiresComp = Part.makeCompound(patternWires)
        for vert in constructionPoints:
            for wire in patternWiresComp.Wires:
                copy = wire.copy()
                vector = vert.Point - copy.CenterOfGravity #Assume "Fixed" is alignment type"
                if fp.PointArraySettings == "Align To Origin":
                    angle = math.degrees(math.atan2(vert.Point.y, vert.Point.x))
                    copy.rotate(copy.CenterOfGravity, FreeCAD.Vector(0,0,1), angle)
                copy.Placement.translate(vector)
                copies.append(copy)
        roguePoints = self.getRoguePoints(fp,shape)
        copiesAndStationaryComp = Part.makeCompound(copies + stationaryWires + roguePoints)
        roguePoints2 = self.getRoguePoints(fp, copiesAndStationaryComp)

        wires = []
        edges = copiesAndStationaryComp.Edges
        if edges:
            se = Part.sortEdges(edges)
            wires = [Part.Wire(s) for s in se]
        finalShape = Part.makeCompound(wires + roguePoints)
        return finalShape

    def handleRectangularArray(self, fp, shape):
        """makes a rectangular array of desired wires"""
        SketchPlus.wireShapeDict["RectangularWires"] = shape
        constructionCircles = self.fetchConstruction(fp, "RectangularRotationCenter", SketchPlus.RectangularRotationCenterDefault, filter=["Circle","Point","ArcOfCircle"])
                
        if shape.isNull():
            return shape
            
        roguePoints = self.getRoguePoints(fp,shape)


        if fp.RectangularCountX == 1 and fp.RectangularCountY == 1 and fp.RectangularRotation == 0:
            return shape

        patternWires, stationaryWires = self.filterWires(shape, fp.RectangularWires)
        if not patternWires and not fp.RectangularWires:
            patternWires = shape.Wires
            stationaryWires = []
            
        if not patternWires + roguePoints:
            return shape #nothing to pattern

        center = FreeCAD.Vector() #assume "Origin" is center for now
        angle = fp.RectangularRotation.Value
        if angle != 0:
            if fp.RectangularRotationCenter == "Center Of Gravity":
                center = None #have to wait until the array is built to calculate this
            elif fp.RectangularRotationCenter not in SketchPlus.RectangularRotationCenterDefault: #construction geo
                geoCenter = constructionCircles[fp.RectangularRotationCenter].toShape()
                center = FreeCAD.Vector(geoCenter.X, geoCenter.Y, 0) if hasattr(geoCenter,"X") else geoCenter.CenterOfGravity

        #not arraying roguePoints for this function
        copies = []
        xInt = fp.RectangularIntervalX.Value
        yInt = fp.RectangularIntervalY.Value
        for xx in range(fp.RectangularCountX):
            for yy in range(fp.RectangularCountY):
                for wire in patternWires + roguePoints:
                    copy = wire.copy()
                    target = FreeCAD.Vector(xx*xInt, yy*yInt, 0)
                    copy.Placement.translate(target)
                    copies.append(copy)

        if copies and angle:
            copyComp = Part.makeCompound(copies)
            center = copyComp.CenterOfGravity if not center else center
            for copy in copies:
                copy.rotate(center, FreeCAD.Vector(0,0,1), angle)


        copiesAndStationaryComp = Part.makeCompound(copies + stationaryWires)
        roguePoints2 = self.getRoguePoints(fp, copiesAndStationaryComp)

        wires = []
        edges = copiesAndStationaryComp.Edges
        if edges:
            se = Part.sortEdges(edges)
            wires = [Part.Wire(s) for s in se]
        finalShape = Part.makeCompound(wires + roguePoints2)
        return finalShape

    def handleWireStatuses(self, fp, shape):
        """fill in the WireStatuses string property"""
        opens = [idx+1 for idx,w in enumerate(shape.Wires) if not w.isClosed()]
        closed = [idx+1 for idx,w in enumerate(shape.Wires) if w.isClosed()]
        if not opens + closed:
            status = "No wires"
        else:
            status = f"Open Wires: {opens}\n"
            status += f"Closed Wires: {closed}\n"

        fp.WireStatuses = status
        return shape

    def execute(self, fp):

        fp.recompute() # necessary, calls the c++ base class recompute handler
        if not hasattr(self,"wireShapeDict"):
            setattr(self,"wireShapeDict",{})
        fp_plm = fp.Placement
        fp.Placement = FreeCAD.Placement()
        shape = fp.Shape.copy()

        function_dict = {"Wire Order":self.handleWireOrder,
                        "External":self.handleExternal,
                        "Construction":self.handleShowConstruction,
                        "Mirroring":self.handleMirroring,
                        "Path Array":self.handlePathArray,
                        "Point Array":self.handlePointArray,
                        "Polar Array": self.handlePolar,
                        "Rectangular Array": self.handleRectangularArray,
                        "Rotate": self.handleRotate,
                        "Scaling": self.handleScaling,
                        "Uniform Scaling": self.handleUniformScaling,
                        "Offsetting": self.handleOffset,
                        "Face Making": self.handleFaceMaker,
                        "Wire Statuses":self.handleWireStatuses, }

        #for legacy objects
        if not "Rotate" in fp.ExecutionOrder and not "-Rotate" in fp.ExecutionOrder:
            fp.ExecutionOrder = fp.ExecutionOrder + ["Rotate"]

        # do not validate that the func is in the dictionary, but rather
        # just let it throw an exception if it is not
        # only thing checked for is the hyphen, which is a normal condition
        # where the user wishes to suppress the function

        for func in fp.ExecutionOrder:
            if not func.startswith("-"):
                shape = function_dict[func](fp, shape)

        fp.Shape = shape
        fp.Placement = fp_plm

    def __getstate__(self):
        '''avoids error messages about objects not being JSON serializable'''
        return None

    def __setstate__(self,state):
        '''avoid JSON serializable error messages'''
        return None

class WireFilterEditorTask:
    def __init__(self, fp, shape, prop, selectionMode = False):
        def triggerChecked(i): return lambda: self.checked(i)
        def triggerUp(i): return lambda: self.up(i)
        def triggerDown(i): return lambda: self.down(i)
        self.fp = fp
        self.shape = shape
        self.prop = prop
        self.selectionMode = selectionMode
        self.visibility = self.fp.ViewObject.Visibility
        self.camera = FreeCADGui.activeView().getCamera()
        FreeCADGui.activeDocument().activeView().viewTop()
        self.ghost = self.fp.Document.addObject("Part::Feature",f"{self.prop}_Ghost")
        self.fp.ViewObject.Proxy.ghosts.append(self.ghost.Name)
        self.ghost.ViewObject.PointSize = 10
        self.ghost.addProperty("App::PropertyString","Information","Base","information").Information =\
f"""
This is a temporary object used to aid in identifying wires during the editing
of one of the various wire filter integer list properties of the {self.fp.Label}
object.  It can be deleted if the dialog is closed, but it should be automatically
deleted upon closing the dialog.
"""
        self.fp.ViewObject.Visibility = False
        self.form = QtGui.QWidget()
        if not self.selectionMode:
            self.form.setWindowTitle(f"Edit {prop}")
            self.order = self.fixOrder(getattr(fp,prop))
        else:
            self.form.setWindowTitle(f"Select wires")
            self.order = [idx+1 for idx in range(len(self.shape.Wires))]

        defaultsBtn = QtWidgets.QPushButton("Defaults")
        defaultsBtn.clicked.connect(self.setDefaults)
        layout = QtGui.QGridLayout()
        self.form.setLayout(layout)
        layout.addWidget(defaultsBtn,0,0,1,1)
        noneBtn = QtWidgets.QPushButton("None")
        noneBtn.clicked.connect(self.noneBtnClicked)
        layout.addWidget(noneBtn,0,1,1,1)
        allBtn = QtWidgets.QPushButton("All")
        allBtn.clicked.connect(self.allBtnClicked)
        layout.addWidget(allBtn,0,2,1,1)
        layout.addWidget(QtWidgets.QSplitter(),1,0,1,3)
        self.cbs = []
        for idx, wire in enumerate(shape.Wires):
            wireNumber = self.order[idx] if self.order[idx] > 0 else -1 * self.order[idx]
            cb = QtGui.QCheckBox(f"Wire{wireNumber}")
            cb.setObjectName(f"{wireNumber}")
            cb.setChecked(wireNumber in self.order)
            self.cbs.append(cb)
            cb.clicked.connect(triggerChecked(idx))
            layout.addWidget(cb, idx+2, 1,alignment=QtCore.Qt.AlignCenter)
            up = QtGui.QPushButton()
            up.setIcon(QtGui.QIcon(FreeCADGui.getIcon("button_up")))
            up.clicked.connect(triggerUp(idx))
            if not self.selectionMode:
                layout.addWidget(up, idx+2, 0)
            down = QtGui.QPushButton()
            down.setIcon(QtGui.QIcon(FreeCADGui.getIcon("button_down")))
            down.clicked.connect(triggerDown(idx))
            if not self.selectionMode:
                layout.addWidget(down, idx+2, 2)
        self.relabel()

    def fixOrder(self, order):
        """The wire order might be incomplete, so we fill in any missing wires
        while retaining the current order"""

        available = [idx + 1 for idx in range(len(self.shape.Wires))]
        available += [-av for av in available]
        default = [i+1 for i in range(len(self.shape.Wires))]
        if 0 in order: #0 signals to filter out all wires
            return [-d for d in default]
        order = [o for o in order if not o == 0 and o in available]
        if len(list(set(order))) == len(self.shape.Wires):
            return order

        order = list(set(order)) #remove any duplicates
        default = [d for d in default if d not in order and -d not in order]
        return order + default

    def setDefaults(self):
        self.order = [d+1 for d in range(len(self.shape.Wires))]
        self.relabel()
        
    def allBtnClicked(self):
        self.order = [d+1 for d in range(len(self.shape.Wires))]
        self.relabel()
        
    def noneBtnClicked(self):
        self.order = [-(d+1) for d in range(len(self.shape.Wires))]
        self.relabel()

    def up(self, idx):
        cb = self.cbs[idx]
        wireNumber = int(cb.objectName())
        wireNumber *= -1 if not wireNumber in self.order else 1
        where = self.order.index(wireNumber)
        if where == 0:
            self.order = self.order[1:] + [self.order[0]]
        else:
             self.order[where], self.order[where - 1] = self.order[where - 1], self.order[where]
        self.relabel()

    def down(self, idx):
        cb = self.cbs[idx]
        wireNumber = int(cb.objectName())
        wireNumber *= -1 if not wireNumber in self.order else 1
        where = self.order.index(wireNumber)
        if self.order[-1] == wireNumber:
            self.order = [self.order[-1]] + self.order[:-1]
        else:
             self.order[where], self.order[where + 1] = self.order[where + 1], self.order[where]
        self.relabel()

    def relabel(self):
        for idx,o in enumerate(self.order):
            self.cbs[idx].setText(f"Wire{o if o > 0 else -o}")
            self.cbs[idx].setObjectName(f"{o if o > 0 else -o}")
            self.cbs[idx].setChecked(o > 0)
        used,unused = self.fp.Proxy.filterWires(self.shape, self.order)
        roguePoints = self.fp.Proxy.getRoguePoints(self.fp, self.shape)
        if self.ghost:
            self.ghost.Shape = Part.makeCompound(used + roguePoints) if bool(used + roguePoints) else Part.Shape()

    def checked(self, idx):
        val = self.order[idx]
        cb = self.cbs[idx]
        wireNumber = val * ( -1 if val < 0 else 1)
        val = wireNumber * (1 if cb.isChecked() else -1)
        self.order[idx] = val
        self.relabel()

    def reject(self):
        self.fp.ViewObject.Visibility = self.visibility
        FreeCADGui.activeView().setCamera(self.camera)
        self.fp.ViewObject.Proxy.closePanel()

    def accept(self):

        self.fp.ViewObject.Visibility = self.visibility
        FreeCADGui.activeView().setCamera(self.camera)
        FreeCADGui.ActiveDocument.resetEdit()
        if not self.selectionMode:
            self.fp.Document.openTransaction(f"edit {self.prop}")
            setattr(self.fp, self.prop, self.order)
            self.fp.Document.recompute()
            self.fp.Document.commitTransaction()
        else:
            used = [o for o in self.order if o > 0]
            self.fp.ViewObject.Proxy.selectWires(used)
        self.fp.ViewObject.Proxy.closePanel()

class ExecutionOrderTask:
    def __init__(self, fp):
        def triggerChecked(i): return lambda: self.checked(i)
        def triggerUp(i): return lambda: self.up(i)
        def triggerDown(i): return lambda: self.down(i)
        self.fp = fp
        self.prop = "ExecutionOrder"
        self.default = SketchPlus.ExecutionOrderDefault
        self.form = QtGui.QWidget()
        self.form.setWindowTitle(f"Edit execution order")
        #self.order will be a list of strings
        self.order = getattr(self.fp, self.prop)

        defaultsBtn = QtWidgets.QPushButton("Defaults")
        defaultsBtn.clicked.connect(self.setDefaults)
        layout = QtGui.QGridLayout()
        self.form.setLayout(layout)
        layout.addWidget(defaultsBtn,0,0,1,1)
        self.cbs = []
        for idx, funcName in enumerate(self.order):
            cb = QtGui.QCheckBox(funcName)
            cb.setChecked(not funcName.startswith("-"))
            self.cbs.append(cb)
            cb.clicked.connect(triggerChecked(idx))
            layout.addWidget(cb, idx+1, 1,alignment=QtCore.Qt.AlignCenter)
            up = QtGui.QPushButton()
            up.setIcon(QtGui.QIcon(FreeCADGui.getIcon("button_up")))
            up.clicked.connect(triggerUp(idx))
            layout.addWidget(up, idx+1, 0)
            down = QtGui.QPushButton()
            down.setIcon(QtGui.QIcon(FreeCADGui.getIcon("button_down")))
            down.clicked.connect(triggerDown(idx))
            layout.addWidget(down, idx+1, 2)
        self.relabel()

    def getStandardButtons(self):
        return (QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel | QtWidgets.QDialogButtonBox.Apply)

    def setDefaults(self):
        self.order = SketchPlus.ExecutionOrderDefault
        self.relabel()

    def up(self, idx):
        cb = self.cbs[idx]
        funcName = cb.text()
        where = self.order.index(funcName)
        if where == 0:
            self.order = self.order[1:] + [self.order[0]]
        else:
             self.order[where], self.order[where - 1] = self.order[where - 1], self.order[where]
        self.relabel()

    def down(self, idx):
        cb = self.cbs[idx]
        funcName = cb.text()
        where = self.order.index(funcName)
        if self.order[-1] == funcName:
            self.order = [self.order[-1]] + self.order[:-1]
        else:
             self.order[where], self.order[where + 1] = self.order[where + 1], self.order[where]
        self.relabel()

    def relabel(self):
        for idx,o in enumerate(self.order):
            self.cbs[idx].setText(o)
            self.cbs[idx].setChecked(not o.startswith('-'))

    def checked(self, idx):
        funcName = self.order[idx]
        cb = self.cbs[idx]
        funcName = "-" + funcName if not cb.isChecked() else funcName[1:]
        self.order[idx] = funcName
        self.relabel()

    def clicked(self, button):
        if button == QtWidgets.QDialogButtonBox.Apply:
           self.fp.Document.openTransaction("Apply execution order edit")
           self.fp.ExecutionOrder = self.order
           self.fp.Document.recompute()
           self.fp.Document.commitTransaction()

    def reject(self):
        self.fp.ViewObject.Proxy.closePanel()
        
    def accept(self):
        FreeCADGui.ActiveDocument.resetEdit()
        self.fp.Document.openTransaction("Edit execution order")
        setattr(self.fp, self.prop, self.order)
        self.fp.Document.recompute()
        self.fp.Document.commitTransaction()
        self.fp.ViewObject.Proxy.closePanel()
        
class AddExternalTask:
    def __init__(self, fp):
        self.fp = fp
        self.form = QtWidgets.QWidget()
        self.form.setWindowTitle("Add external geometry")
        layout = QtWidgets.QVBoxLayout()
        msg = \
"""
Selection options:

Whole object = add all vertices
Face = add all edges of that face
Edges or Vertices = add those subobjects.

OK = add selected
"""
        label = QtWidgets.QLabel(msg)
        layout.addWidget(label)
        self.form.setLayout(layout)
        FreeCADGui.Selection.clearSelection()
        
    def addExternalLink(self, objName, subName):
        """add the external link"""
        try:
            self.fp.addExternal(objName, subName)
            FreeCAD.Console.PrintMessage(f"{self.fp.Label}: added {objName}:{subName} as external geometry.\n")
        except:
            FreeCAD.Console.PrintMessage(f"{self.fp.Label}: error adding {objName}, {subName}, skipping\n")
        
    def addExternalLinks(self):
        """add the external links selected by the user"""
        def isSameEdge(edge1, edge2):
            # if edge1.Length != edge2.Length:
            #     return False
            if len(edge1.Vertexes) != len(edge2.Vertexes):
                return False
            for idx,vert in enumerate(edge1.Vertexes):
                if not edge1.Vertexes[idx].Point.isEqual(edge2.Vertexes[idx].Point, Part.Precision.confusion()):
                    return False
            return True
                
        sel = FreeCADGui.Selection.getCompleteSelection()
        if not sel:
            return
        pbOuter = gu.MRProgress(len(sel)) if len(sel) > 100 else None
        for s in sel:
            if pbOuter and pbOuter.isCanceled():
                break
            obj = s.Object
            while obj.isDerivedFrom("App::Link"):
                obj = obj.LinkedObject
            shape = Part.getShape(obj)
            if obj.Document != s.Object.Document:
                FreeCAD.Console.PrintError("Linked object must reside in same document.\n")                
            if not shape or shape.isNull():
                continue
            subNames = s.SubElementNames
            if not subNames[0]: #get all the vertices
                pb = gu.MRProgress(len(shape.Vertexes),"Stop") if len(shape.Vertexes) > 100 else None
                for idx,vert in enumerate(shape.Vertexes):
                    if pb and pb.isCanceled():
                        break
                    self.addExternalLink(obj.Name, f"Vertex{idx+1}")

            elif "Face" in subNames[0]:
                face = s.SubObjects[0]
                pb = gu.MRProgress(len(face.Edges),"Stop") if len(face.Edges) > 100 else None
                if len(face.Edges) > 100:
                    pb = gu.MRProgress(len(face.Edges),"Stop")
                for idx, edge in enumerate(face.Edges):
                    if pb and pb.isCanceled():
                        break
                    for idx2, edge2 in enumerate(shape.Edges):
                        if isSameEdge(edge, edge2):
                            self.addExternalLink(obj.Name, f"Edge{idx2+1}")
                            break
            else:
                self.addExternalLink(obj.Name, subNames[0])
        
        
    def accept(self):
        self.fp.Document.openTransaction("Add external links")
        self.addExternalLinks()
        self.fp.Document.commitTransaction()
        FreeCADGui.Control.closeDialog()

    def reject(self):
        FreeCADGui.Control.closeDialog()
        
#for handling app::propertyfloatlist properties in the dialogs
class CustomListWidget(QtWidgets.QListWidget):
    def __init__(self, parent=None, float_list= []):
        super(CustomListWidget, self).__init__(parent)
        self.add_default_items(float_list)
        self.itemClicked.connect(self.handle_item_click)
        self.itemDoubleClicked.connect(self.edit_item_value)
        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.showContextMenu)
        
    def edit_item_value(self, item):
        if item.text() != "+":
            current_value = float(item.text())
            new_value, ok = QtGui.QInputDialog.getDouble(self, "Edit Value", 
                                                   "Enter new value:", 
                                                   current_value, 
                                                   -10000, 10000, 2)
            if ok:
                item.setText(str(new_value))

    def add_default_items(self, floats):
        for f in floats:
            item = QtWidgets.QListWidgetItem(f"{f}")
            self.addItem(item)
       
        plus_item = QtWidgets.QListWidgetItem("+")
        plus_item.setFlags(plus_item.flags() & ~QtCore.Qt.ItemIsSelectable)
        self.addItem(plus_item)

    def showContextMenu(self, pos: QtCore.QPoint):
        # Find the item at the clicked position
        item = self.itemAt(pos)

        if item:
            idx = self.row(item)
            context_menu = QtWidgets.QMenu(self)
            delete_action = QtGui.QAction("Delete", self)
            delete_action.triggered.connect(lambda: self.removeItem(item))
            if item.text() == "+":
                delete_action.setEnabled(False)
            context_menu.addAction(delete_action)
            add_action = context_menu.addAction("Add item")
            add_action.triggered.connect(self.addNewItem)
            edit_action = context_menu.addAction("Edit item")
            edit_action.triggered.connect(lambda: self.edit_item_value(item))
            if item.text() == "+":
                edit_action.setEnabled(False)
            move_up_action = context_menu.addAction("Move up")    
            move_up_action.triggered.connect(lambda: self.move_up(item))
            if item.text() == "+":
                move_up_action.setEnabled(False)
            move_down_action = context_menu.addAction("Move down")
            move_down_action.triggered.connect(lambda: self.move_down(item))
            if item.text() == "+":
                move_down_action.setEnabled(False)
                
            context_menu.exec_(self.mapToGlobal(pos))
            
    def addNewItem(self):
        new_item = QtWidgets.QListWidgetItem("0.0")
        self.insertItem(self.count() - 1, new_item)
        self.setCurrentItem(new_item)
        
    def move_up(self, item):
        row = self.row(item)
        if row > 0:
            previous = row
        else:
            row == 0
            if self.count() > 2:
                previous = self.count() - 1
        removed = self.takeItem(row)
        self.insertItem(previous - 1, removed)
        
    def move_down(self, item):
        row = self.row(item)
        if row > self.count() - 1:
            next = 0
        else:
            next = row + 1
            if self.item(next).text() == "+":
                next = 0

        removed = self.takeItem(row)
        self.insertItem(next, removed)
        
        
    def handle_item_click(self, item):
        if item.text() == "+":
            new_item = QtWidgets.QListWidgetItem("0.0")
            self.insertItem(self.count() - 1, new_item)
            self.setCurrentItem(new_item)

    def removeItem(self, item: QtWidgets.QListWidgetItem):
        row = self.row(item)
        self.takeItem(row) 
    
    def value(self):
        float_list = []
        for index in range(self.count()):
            item = self.item(index)
            if item.text() == "+":
                continue
            try:
                float_value = float(item.text())
                float_list.append(float_value)
            except ValueError:
                print(f"Error: Item '{item.text()}' is not a valid float.")
        return float_list

class FunctionTask:
    def __init__(self, fp, funcType):
        self.fp = fp
        self.form = QtWidgets.QWidget()
        self.funcType = funcType
        self.form.setWindowTitle(f"Setup {funcType} function")
        self.props = self.getProps(self.funcType)
        layout = QtWidgets.QVBoxLayout()
        self.form.setLayout(layout)
        if self.funcType == "Point":
            constructionPoints = []
            for idx,geo in enumerate(self.fp.Geometry):
                if self.fp.getConstruction(idx):
                    if "Point" in geo.TypeId:
                        constructionPoints.append(geo)
            
            layout.addWidget(QtWidgets.QLabel(f"Construction points count: {len(constructionPoints)}"))

        for prop in self.props:
            hbox = QtWidgets.QHBoxLayout()
            hbox.addWidget(QtWidgets.QLabel(prop[0]))
            hbox.addWidget(prop[2])
            layout.addLayout(hbox)
            
    def accept(self):
        for prop in self.props:
            val = None
            task = prop[4]
            propName = prop[0]
            if hasattr(task, "order"):
                val = task.order
                QtCore.QTimer().singleShot(10, task.reject)
            else:
                val = prop[3]()
            setattr(self.fp, propName, val)
        self.fp.Document.recompute()
        self.fp.ViewObject.Proxy.closePanel()
        
    def reject(self):
        for prop in self.props:
            task = prop[4]
            if hasattr(task, "reject"):
                task.reject()
        self.fp.ViewObject.Proxy.closePanel()
        
    def getPathArrayWires(self):
        """we don't handle patharraywires here because the path array requires
        2 different wire arrays, one for the path, and one for the wires to pattern"""
        return getattr(self.fp, "PathArrayWires")
        
    def editPathArrayWires(self):
        msg_box = QtWidgets.QMessageBox()
        msg_box.setWindowTitle("Exiting dialog")
        msg_box.setIcon(QtWidgets.QMessageBox.Question)
        msg_box.setText("Apply any changes to the SketchPlus object properties before closing?")
        msg_box.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        msg_box.setDefaultButton(QtWidgets.QMessageBox.No)
        yes_button = msg_box.button(QtWidgets.QMessageBox.Yes)
        no_button = msg_box.button(QtWidgets.QMessageBox.No)
        yes_button.setText("Apply")
        no_button.setText("Discard")
        result = msg_box.exec_()
        if result == QtWidgets.QMessageBox.Yes:
            self.accept()
        else:
            self.reject()
        QtCore.QTimer().singleShot(100, self.doEdit)
        
    def doEdit(self):
        self.fp.ViewObject.Proxy.editWireFilter("PathArrayWires")
        
    def getProps(self, funcType):
        """returns a list of tuples:
        (propName, propType, widget, getter, task)"""
        def makeWidget(propName, propType):
            widget = None
            task = None
            if propType == "App::PropertyDistance":
                widget = QtWidgets.QDoubleSpinBox()
                widget.setRange(-float("inf"),float("inf"))
                widget.setValue(getattr(self.fp, propName))
                getter = widget.value
            elif propType == "App::PropertyBool":
                widget = QtWidgets.QCheckBox()
                widget.setChecked(getattr(self.fp, propName))
                getter = widget.isChecked
            elif propType == "App::PropertyIntegerConstraint":
                widget = QtWidgets.QSpinBox()
                if "Rectangular" in propName:
                    widget.setRange(1, 100000)
                elif "Polar" in propName or "Path" in propName:
                    widget.setRange(0, 100000)
                widget.setValue(getattr(self.fp, propName))
                getter = widget.value
            elif propType == "App::PropertyFloatList":
                widget = CustomListWidget(parent=self.form, float_list = getattr(self.fp, propName))
                getter = widget.value
            elif propType == "App::PropertyFloat":
                widget = QtWidgets.QDoubleSpinBox()
                widget.setRange(-float("inf"),float("inf"))
                widget.setValue(getattr(self.fp, propName))
                getter = widget.value
            elif propType == "App::PropertyAngle":
                widget = QtWidgets.QDoubleSpinBox()
                widget.setRange(-float("inf"),float("inf"))
                widget.setValue(getattr(self.fp, propName))
                getter = widget.value
            elif propType == "App::PropertyEnumeration":
                widget = QtWidgets.QComboBox()
                widget.addItems(self.fp.getEnumerationsOfProperty(propName))
                widget.setEditable(False)
                widget.setCurrentText(getattr(self.fp, propName))
                getter = widget.currentText
            elif propType == "App::PropertyIntegerList":
                if propName != "PathArrayWires":
                    if not propName in SketchPlus.wireShapeDict: #doc reloaded
                        self.fp.recompute()
                        self.fp.Document.recompute()
                    task = WireFilterEditorTask(self.fp, SketchPlus.wireShapeDict[propName], propName)
                    widget = task.form
                    getter = task.order
                else:
                    widget = QtWidgets.QPushButton("Edit PathArrayWires")
                    widget.setToolTip("Closes this dialog and opens the PathArrayWires filter editor")
                    widget.clicked.connect(self.editPathArrayWires)
                    getter = self.getPathArrayWires
            else:
                FreeCAD.Console.PrintError(f"unsupported property type: {propType}, propName = {propName}")
            widget.setObjectName(f"{propName}")
            return (widget, getter, task)
            
        propNames = [prop for prop in self.fp.PropertiesList if prop.startswith(funcType)]
        if "Scale" in funcType and not "Uniform" in funcType:
            propNames = [p for p in propNames if not "Uniform" in p]
        elif "ScaleUniform" in funcType:
            propNames = [p for p in propNames if "Uniform" in p]
        propTypes  = [self.fp.getTypeIdOfProperty(p) for p in propNames]
        propTuples = []
        for idx,propName in enumerate(propNames):
            widget,getter,task = makeWidget(propNames[idx], propTypes[idx])
            propTuples.append((propNames[idx], propTypes[idx], widget, getter, task))
        return propTuples
        
class SketchPlusVP:
    def __init__(self, vobj):
        vobj.Proxy = self
        self.ghosts = []

    def attach(self, vobj):
        self.Object = vobj.Object

    def setupContextMenu(self, vobj, menu):
        def makeSelTrigger(i): return lambda : self.selectWire(self.Object, i)
        def makeWFTrigger(key): return lambda : self.editWireFilter(key)
        def makeFuncTrigger(func): return lambda: self.setupFunction(func)
        fp = self.Object
        if not hasattr(self,"ghosts"):
            self.ghosts = []
        
        funcMenu = menu.addMenu("Setup array or function")
        funcs = sorted(["Show","Path", "Point", "Polar", "Rectangular", "Rotate", "Offset", "Scale",
                "ScaleUniform","Mirror","FaceMaker"])
        for func in funcs:
            funcAction = funcMenu.addAction(func)
            funcAction.triggered.connect(makeFuncTrigger(func))
        
        if fp.Shape.Wires:
            selWireAction = menu.addAction("Select Wires")
            selWireAction.triggered.connect(self.setupWireSelection)

        wireFilterMenu = None
        if hasattr(fp.Proxy,"wireShapeDict"):
            for k,v in sorted(SketchPlus.wireShapeDict.items()):
                if v and not v.isNull():
                    if not wireFilterMenu:
                        wireFilterMenu = menu.addMenu("Edit wire filter properties")
                    wfAction = wireFilterMenu.addAction(k)
                    wfAction.triggered.connect(makeWFTrigger(k))

        execAction = menu.addAction("Edit execution order")
        execAction.triggered.connect(self.editExecutionOrder)
        
        addExternalAction = menu.addAction("Add external geometry")
        addExternalAction.triggered.connect(self.addExternal)
        
        hideDependentAction = menu.addAction("Hide dependencies")
        hideDependentAction.triggered.connect(self.hideDependencies)
        
        oefaction = menu.addAction("Check with OpenEdgeFinder")
        oefaction.triggered.connect(self.openEdgeFinder)
        
    def openEdgeFinder(self):
        cmd = MeshRemodelOpenEdgeFinderCommandClass()
        cmd.Activated()
        self.Object.Document.recompute()
        
    def closePanel(self):
        fp = self.Object
        if hasattr(self, "ghosts"):
            for ghost in self.ghosts:
                obj = fp.Document.getObject(ghost)
                if obj:
                    fp.Document.removeObject(obj.Name)
            self.ghosts = []
        else:
            print("SketchPlus object has no ghosts")
        FreeCADGui.Control.closeDialog()
        fp.ViewObject.Visibility = True

    def openPanel(self, panel):
        import time
        start = time.time()
        timedOut = False
        while FreeCADGui.Control.activeDialog():
            FreeCAD.Console.PrintMessage(".")
            FreeCADGui.updateGui()
            time.sleep(0.1)
            now = time.time()
            if now - start > 15:
                timedOut = True
                break
        if not timedOut and not FreeCADGui.Control.activeDialog():
            FreeCADGui.Control.showDialog(panel)
        else:
            FreeCAD.Console.PrintError("Another task panel is already open\n")

        
    def setupFunction(self, funcType):
        fp = self.Object
        panel = FunctionTask(fp, funcType)
        self.openPanel(panel)
        
    def hideDependencies(self):
        fp = self.Object
        for obj in fp.OutListRecursive:
            if hasattr(obj, "ViewObject") and obj.ViewObject.Visibility == True:
                obj.ViewObject.Visibility = False
        
    def addExternal(self):
        fp = self.Object
        panel = AddExternalTask(fp)
        self.openPanel(panel)
        
    def editExecutionOrder(self):
        """Edit the execution order string"""
        fp = self.Object
        panel = ExecutionOrderTask(fp)
        self.openPanel(panel)

    def setupWireSelection(self):
        fp = self.Object
        shape = fp.Shape
        panel = WireFilterEditorTask(fp, shape, "", selectionMode = True)
        self.openPanel(panel)

    def editWireFilter(self, key):
        fp = self.Object
        shape = SketchPlus.wireShapeDict[key]
        panel = WireFilterEditorTask(fp, shape, key)
        self.openPanel(panel)

    def selectWires(self, order):
        FreeCADGui.Selection.clearSelection()
        for o in order:
            self.selectWire(self.Object, o-1)

    def selectWire(self, fp, idx):
        """select the edges of Wire{idx+1} or Wires[idx]"""
        wireEdges = fp.Shape.Wires[idx].Edges
        shapeEdges = fp.Shape.Edges
        subnames = []
        for we in wireEdges:
            for ii,se in enumerate(shapeEdges):
                if se.isEqual(we):
                    #print(f"Match found: Edge{ii+1} is in Wire{idx+1}")
                    subnames.append(f"Edge{ii+1}")

        for subname in subnames:
            FreeCADGui.Selection.addSelection(fp.Document.Name, fp.Name, subname)

    def __getstate__(self):
        '''When saving the document this object gets stored using Python's json module.\
                Since we have some un-serializable parts here -- the Coin stuff -- we must define this method\
                to return a tuple of all serializable objects or None.'''
        return None
    def __setstate__(self,state):
        '''When restoring the serialized object from document we have the chance to set some internals here.\
                Since no data were serialized nothing needs to be done here.'''
        return None



# Make a sketch plus object
class MeshRemodelCreateSketchCommandClass(object):
    """Create sketch plus object"""
    #class variables accessed as MeshRemodelCreateSketchCommandClass.sketchPlus, etc.
    sketchPlus = None
    activeObject = None
    activeObjectOriginVisibility = False

    def __init__(self):
        self.subs = []

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'CreateSketch.svg') ,
            'MenuText': "Create s&ketch plus" ,
            'ToolTip' : fixTip("\
Create a new SketchPlus object\n\
")}

    def Activated(self):
        
        def addToActiveObject(obj):
            """check for existing active Body or App::Part and add object to it if found"""

            body = FreeCADGui.ActiveDocument.ActiveView.getActiveObject("pdbody")
            part = FreeCADGui.ActiveDocument.ActiveView.getActiveObject("part")
            MeshRemodelCreateSketchCommandClass.activeObject = body if body else part if part else None
            if body:
                if not obj in body.Group:
                    body.Group += [obj]
                    MeshRemodelCreateSketchCommandClass.activeObjectOriginVisibility = \
                                                    body.Origin.ViewObject.Visibility
                    body.Origin.ViewObject.Visibility = True
            elif part:
                if not obj in part.Group:
                    part.Group += [obj]
                    MeshRemodelCreateSketchCommandClass.activeObjectOriginVisibility = \
                                                    part.Origin.ViewObject.Visibility
                    part.Origin.ViewObject.Visibility = True
        
        def attached():
            """a callback for the attachment function, when OK or Cancel is clicked, opens sketch in editor"""
            obj = MeshRemodelCreateSketchCommandClass.sketchPlus
            if hasattr(obj,"ViewObject"): #might be None
                QtCore.QTimer.singleShot(100, obj.ViewObject.doubleClicked) # give time for other dialog to close
            if MeshRemodelCreateSketchCommandClass.activeObject:
                MeshRemodelCreateSketchCommandClass.activeObject.Origin.Visibility = \
                        MeshRemodelCreateSketchCommandClass.activeObjectOriginVisibility
        
        def attachment(obj):
            """open attachment dialog if there are selections"""
            MeshRemodelCreateSketchCommandClass.sketchPlus = obj
            take_selection = True
            sel = FreeCADGui.Selection.getCompleteSelection()
            if not sel and not MeshRemodelCreateSketchCommandClass.activeObject:
                return
            if sel and not hasattr(sel[0].Object,"Shape"):
                FreeCAD.Console.PrintMessage(f"{sketchPlus.Label}: {sel[0].Object.Label} is not a part feature, attachment dialog is being skipped.\n")
                if MeshRemodelCreateSketchCommandClass.activeObject:
                   MeshRemodelCreateSketchCommandClass.activeObject.Origin.ViewObject.Visibility = \
                        MeshRemodelCreateSketchCommandClass.activeObjectOriginVisibility
                return
            if sel and sel[0].Object in obj.InListRecursive:
                take_selection = False
            import AttachmentEditor.TaskAttachmentEditor as TaskAttachmentEditor
            taskd = TaskAttachmentEditor.AttachmentEditorTaskPanel(obj,
                                                                    take_selection = take_selection,
                                                                    create_transaction= True,
                                                                    callback_OK = attached,
                                                                    callback_Cancel = attached,
                                                                    callback_Apply = None)
            FreeCADGui.Control.showDialog(taskd)
            #taskd.getCurrentMode()
            suggested = None
            item = None
            if taskd.last_sugr is not None:
                if taskd.last_sugr['message'] == 'OK':
                    suggested = taskd.last_sugr['bestFitMode']
        
            if suggested:
                suggested = taskd.attacher.getModeInfo(suggested)['UserFriendlyName']
                listOfModes = taskd.form.listOfModes
                item = listOfModes.findItems(suggested, QtCore.Qt.MatchExactly)
                if item:
                    listOfModes.setCurrentItem(item[0])
        
        
        doc = FreeCAD.ActiveDocument if FreeCAD.ActiveDocument else FreeCAD.newDocument()
        obj = doc.addObject("Sketcher::SketchObjectPython","SketchPlus")
        SketchPlus(obj)
        SketchPlusVP(obj.ViewObject)
        addToActiveObject(obj)
        sel = FreeCADGui.Selection.getSelection()
        skipAttach = False
        if len(sel) == 1:
            if sel[0].isDerivedFrom("Sketcher::SketchObject"):
                items = ["Duplicate selected sketch", "Do not duplicate"]
                item, ok = QtGui.QInputDialog.getItem(FreeCADGui.getMainWindow(), "Duplicate selected sketch", f"Do you wish to duplicate the selected sketch: {sel[0].Label}?", items)
                if ok and item == items[0]:
                    skipAttach = True
                    sk = sel[0]
                    enumerationList = ["MirrorAxis", "PolarCenter"," ScaleUniformCenter",
                                        "PointArraySettings","RectangularRotationCenter"]
                    blacklist = ["Constraints","Geometry","ExpressionEngine","Proxy","Version","Visibility","Label","Label2",
                                "PolarCenter","ScaleUniformCenter","MirrorAxis","PointArraySettings",
                                "RectangularRotationCenter", "ExternalGeometry" ]
                    for prop in sk.PropertiesList:
                        if prop in blacklist:
                            continue
                        if not "ReadOnly" in sk.getEditorMode(prop):
                            # print(f"copying {prop}")
                            setattr(obj,prop,getattr(sk,prop))
                            # print(f"{prop} copied")
                    doc.recompute()
                    for geoid,geo in enumerate(sk.Geometry):
                        obj.addGeometry(geo)
                        obj.setConstruction(geoid, sk.getConstruction(geoid))
                    obj.Constraints = sk.Constraints
                    for prop in enumerationList:
                        if not hasattr(sk, prop):
                            continue
                        if not getattr(sk, prop) in obj.getEnumerationsOfProperty(prop):
                            setattr(obj, prop, sk.getEnumerationsOfProperty(prop))
                        setattr(obj,prop, getattr(sk, prop))
                    doc.recompute()
                    for extGeo in sk.ExternalGeometry:
                        link, subs = extGeo
                        for sub in subs:
                            obj.addExternal(link.Name, sub)
                    doc.recompute()
                    skParent = sk.getParentGeoFeatureGroup()
                    if skParent:
                        if hasattr(skParent,"Group"):
                            if not obj in skParent.Group:
                                skParent.Group = skParent.Group + [obj]
    
                if MeshRemodelCreateSketchCommandClass.activeObject:
                    MeshRemodelCreateSketchCommandClass.activeObject.Origin.ViewObject.Visibility = \
                            MeshRemodelCreateSketchCommandClass.activeObjectOriginVisibility
    
        if not skipAttach:
            attachment(obj)
        doc.recompute()


    def IsActive(self):
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
        self.lastObj = None
        self.lastSubs = []

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
        
    def addNearestPoint(self, skipOne = False):
        """called when user has Ctrl+(Shift+)Clicked on the Cancel objserver button on the statusbar"""
        if not self.lastObj or not global_picked:
            FreeCAD.Console.PrintError("No last object, add a selected vertex from an object and try again.\n")
            return
        
        verts = self.lastObj.Shape.Vertexes
        pts = [v.Point for v in verts]
        pt = global_picked[-1]
        nearest = gu.nearestPoint(pt, pts, global_picked, skipOne)
     #   FreeCAD.Console.PrintMessage(f"nearest = {nearest}\n")
        if nearest.isEqual(pt, Part.Precision.confusion()):
            FreeCAD.Console.PrintMessage("No more vertices to add for this object\n")
            return
        idx = pts.index(nearest)
       # FreeCAD.Console.PrintMessage(f"idx = {idx}\n")
        self.setPreselection(self.lastObj.Document.Name, self.lastObj.Name, f"Vertex{idx + 1}")
        
    def setPreselection(self,doc,obj,sub):                # Preselection object
        modifiers = QtGui.QApplication.keyboardModifiers()
        if not modifiers & QtCore.Qt.ControlModifier:
            return
        self.lastObj = FreeCAD.ActiveDocument.getObject(obj)
        self.lastSubs.append(sub)
        if self.isVertexMode() and "Vertex" in str(sub):
            Gui.Selection.addSelection(doc,obj,str(sub))
            idx = int(sub[6:])
            thisobj = FreeCAD.ActiveDocument.getObject(obj)
            p = thisobj.Shape.Vertexes[idx-1].Point
            if not gu.hasPoint(p,global_picked,.0001):
                global_picked.append(p)
                self.lastSubs.append(sub)
               # FreeCAD.Console.PrintMessage(f"Adding {sub} to selection\n")
                self.writeToConsole(doc, obj, sub)
            # else:
            #     FreeCAD.Console.PrintMessage(f"skipping {sub} as it is already in global_picked\n")

        elif self.isEdgeMode() and "Edge" in str(sub):
            Gui.Selection.addSelection(doc,obj,str(sub))
            
    def writeToConsole(self, doc, obj, sub):
        """write to the python console as a comment so the go back selection tool
        can work with these selections"""
        Gui.doCommandGui(f"# Gui.Selection.addSelection('{doc}','{obj}','{sub}')\n")

    def removeLastPoint(self):
        """remove the last point added to global_picked"""
        if not global_picked:
            return
        if not self.lastObj:
            return
        if not self.lastSubs:
            return
        global_picked.pop()
        Gui.doCommandGui(f"Gui.Selection.removeSelection('{self.lastObj.Document.Name}', '{self.lastObj.Name}', '{self.lastSubs[-1]}')\n")
        self.lastSubs.pop()
        

    def addSelection(self,doc,obj,sub,pnt):               # Selection object
        if not "Vertex" in sub:
            return
        v = FreeCAD.ActiveDocument.getObject(obj).getSubObject(sub)
        point = v.Point if hasattr(v,"Point") else None
        if not point:
            # FreeCAD.Console.PrintMessage(f"not adding {v} to global_picked list, obj = {obj}, sub = {sub}\n")
            return
        if not gu.hasPoint(point, global_picked, 0.0001):
            global_picked.append(point)
            self.lastObj = FreeCAD.ActiveDocument.getObject(obj)
            if self.lastSubs and self.lastSubs[-1] != sub:
                self.lastSubs.append(sub)
           # FreeCAD.Console.PrintMessage(f"added {point} = sub {sub}\n")
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
                            if len(self.lastSubs) >= ii and sub == self.lastSubs[ii]:
                                self.lastSubs.remove(self.lastSubs[ii])
                            FreeCAD.Console.PrintMessage(f"Removing {sub} from selection\n")
                    except:
                        pass
        pass

    def setSelection(self,doc):                           # Selection in ComboView
        #App.Console.PrintMessage("setSelection"+ "\n")
        pass

    def clearSelection(self,doc):                         # If click on the screen, clear the selection
        #FreeCAD.Console.PrintMessage("clearSelection"+ "\n")  # If click on another object, clear the previous object
        global_picked.clear()
        self.lastObj = None
        self.lastSubs = []


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
            'ToolTip' : \
"""
Allows to use Ctrl + Preselect to add vertices to the selection.  The selected
vertices should work with all MeshRemodel commands, but not necessarily with
other workbenches.

If the selection observer is already installed:
Click = remove observer
Ctrl+Click = add nearest point to last point
Ctrl+Shift+Click = add 2nd nearest point to last point
Alt+Click = remove last added point
"""}

    def on_clicked(self):
        modifiers = QtGui.QApplication.keyboardModifiers()
        if modifiers & QtCore.Qt.ControlModifier:
            skipOne = False
            if modifiers & QtCore.Qt.ShiftModifier:
                skipOne = True
                FreeCAD.Console.PrintMessage("skipping one = true\n")
            self.sel.addNearestPoint(skipOne)
        elif modifiers & QtCore.Qt.AltModifier:
            self.sel.removeLastPoint()
        else:
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
            modifiers = QtGui.QApplication.keyboardModifiers()
            if modifiers & QtCore.Qt.ControlModifier: #add nearest point
                self.on_clicked()
                return
            if modifiers & QtCore.Qt.AltModifier: #remove last point
                self.on_clicked()
                return
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

# end subshapebinder
################################################################################



class OpenEdgeFinder:
    def __init__(self, obj):
        obj.Proxy = self
        grp = "OpenEdgeFinder"
        obj.addProperty("App::PropertyFloat","Tolerance",grp,"Tolerance to use when finding coincident points.").Tolerance = Part.Precision.confusion()
        obj.addProperty("App::PropertyLink","Sketch",grp,"Sketch to find open edges, also works with other object types").Sketch = None
        obj.addProperty("App::PropertyVectorList","Points",grp,"List of points in linked (Sketch) object.").Points = []
        obj.addProperty("App::PropertyVectorList","IssuePoints",grp,"List of the points with issues").IssuePoints = []
        obj.addProperty("App::PropertyColor","PointColor",grp,"Color for open points").PointColor = (255,0,0)
        obj.addProperty("App::PropertyFloat","PointSize",grp,"Sphere radius for open end markers, 0 = attempt auto scaling").PointSize = 0
        obj.addProperty("App::PropertyLinkHidden","Container",grp,"Body or Part parent object, if any")
        obj.addProperty("App::PropertyString","Version",grp,"Version used to create this object").Version = __version__
        obj.addProperty("App::PropertyBool","FindTSections",grp,"If True, also check for T sections with more than 2 points at a location").FindTSections=True
        obj.addProperty("App::PropertyBool","FindCrossSections", grp, "If True, also check for cross sections").FindCrossSections = True
        obj.addProperty("App::PropertyBool","FindOverlaps",grp,"If True, also check for and show overlapping edges").FindOverlaps = True
        obj.addProperty("App::PropertyBool","FindStrayPoints",grp,"If True, find non-construction mode points").FindStrayPoints = True
        obj.addProperty("App::PropertyFloat","SphereRadius",grp,"Readonly, sphere radius used, which is either PointSize or calculated based on sketch size if PointSize = 0 for auto scaling")
        obj.setEditorMode("SphereRadius", 1)
        obj.addProperty("App::PropertyBool","FindOpenEnds",grp,"If true look for open ended edges").FindOpenEnds = True
        obj.addProperty("App::PropertyString","Report",grp,"Report on findings of sketch validation")
        obj.addProperty("App::PropertyInteger","CountIssues",grp,"Count of how many issues the sketch has").CountIssues = 0
        obj.addProperty("App::PropertyVectorList","OpenPoints",grp,"coordinate list for all the open edge ends found")

    def equal(self, fp, p1, p2):
        """return true if p1 and p2 are within tolerance distance of each other"""
        return p1.isEqual(p2, fp.Tolerance)

    def in_list(self, fp, pt, plist):
        """return true if pt is in plist"""
        for p in plist:
            if self.equal(fp, p, pt):
                return True
        return False


    def near_points(self, fp, pt, pts):
        """return a list of points in pts that are within tolerance distnace of pt"""
        nears = [p for p in pts if self.equal(fp, p, pt)]
        return nears

    def is_on_multiple_edges(self, fp, pt):
        """return true if pt is on multiple edges, note that only open points are passed into this function"""
        on_edges = [edge for edge in fp.Sketch.Shape.Edges if edge.isInside(pt, fp.Tolerance, True)]
        return len(on_edges) > 1

    def stray_points(self, fp):
        """return all non-construction mode points in the sketch"""
        if not hasattr(fp.Sketch, 'Geometry'):
            return [] # not a sketch

        pt_list = []
        for idx,geo in enumerate(fp.Sketch.Geometry):
            if geo.TypeId == "Part::GeomPoint":
                if not fp.Sketch.getConstruction(idx):
                    pt_list.append(FreeCAD.Vector(geo.X, geo.Y, geo.Z))
        return pt_list


    def intersections(self, fp, edges):
        """return a list of self-intersections as vectors among the edges, overlapping edges, if any"""

        count = 0
        try:
            fp.Sketch.Shape.check(True)
        except ValueError as e:
           # print(f"{type(e)}, count = {count}")
            count += 1
        # print(f"count = {count}")
        if not count:
            return ([],[])
        intersections = []
        overlaps = []
        all_vertexes = [v.Point for v in fp.Sketch.Shape.Vertexes]
        for i, edge1 in enumerate(edges):
            for j, edge2 in enumerate(edges):
                if i == j:
                    continue
                intpts = []
                try:
                    intpts = edge1.Curve.intersect(edge2.Curve, fp.Tolerance)
                except Exception as e:
                    if not i in overlaps:
                        if edge1.common(edge2).Length: #do not count line segments on the same line, but not connected as overlaps
                            overlaps.append(i)

                vectors = [FreeCAD.Vector(p.X,p.Y,p.Z) for p in intpts]
                for v in vectors:
                    if edge1.isInside(v, fp.Tolerance, False) and \
                    edge2.isInside(v,fp.Tolerance,False):
                        if not self.in_list(fp, v, intersections):
                            if not self.in_list(fp, v, all_vertexes):
                                intersections.append(v)

        #intersections = [FreeCAD.Vector(p[0], p[1], p[2]) for p in intersections]
        return (intersections, overlaps)

    def scale_spheres(self, fp):
        if fp.PointSize != 0:
            return fp.PointSize #no scaling unless set to 0
        bb = fp.Sketch.Shape.BoundBox
        longest = max(bb.XLength, bb.YLength, bb.ZLength)
        return longest/50

    def execute(self, fp):
        if not fp.Sketch:
            fp.Shape = Part.Shape()
            fp.Label = f"OpenEdgeFinder -- no sketch"
            return
        fp.Label = f"OpenEdgeFionder: checking {fp.Sketch.Label}"
        edges = fp.Sketch.Shape.Edges if hasattr(fp.Sketch,"Shape") else []
        if not edges:
            edges = fp.Sketch.Edges if hasattr(fp.Sketch,"Edges") else []
        if not edges and not fp.FindStrayPoints:
            fp.Shape = Part.Shape()
            return
        # print(f"edges = {edges}")
        pts = []
        for edge in edges:
            first = edge.valueAt(edge.FirstParameter)
            last = edge.valueAt(edge.LastParameter)
            if not self.equal(fp, last, first):
                pts.append(first)
                pts.append(last)

        opens = []
        if fp.FindOpenEnds:
            for p in pts:
                if len(self.near_points(fp, p, pts)) < 2:
                    if not self.in_list(fp, p, opens):
                        opens.append(p)

        if not hasattr(fp,"OpenPoints"):
            fp.addProperty("App::PropertyVectorList","OpenPoints","OpenEdgeFinder","Open edge end points found")
        fp.OpenPoints = opens
        tsections = [p for p in pts if len(self.near_points(fp, p, pts)) > 2] if fp.FindTSections else []
        # extra_ts are only used for the Report, these are open ends that are actually t-sections
        extra_ts = [op for op in opens if self.is_on_multiple_edges(fp, op)] if fp.FindTSections else []
        intersections, overlaps = self.intersections(fp, edges) if fp.FindCrossSections else ([],[])
        strays = self.stray_points(fp) if fp.FindStrayPoints else []
        fp.Points = pts + strays
        issues = opens + tsections + intersections + strays
        count_issues = len(issues) + int(len(overlaps) / 2) #count overlapped edges as 1 issue per every 2 overlapped edges

        if not count_issues:
            fp.Label = f"{count_issues} issues in {fp.Sketch.Label}"
        else:
            s = "s"if count_issues != 1 else ""
            fp.Label = f"{count_issues} issue{s} in {fp.Sketch.Label}"
        fp.IssuePoints = issues
        radius = self.scale_spheres(fp)
        fp.SphereRadius = radius
        verts = [Part.makeSphere(radius, pt) for pt in issues]
        overlapped_edges = [fp.Sketch.Shape.Edges[ii] for ii in overlaps]
        verts += overlapped_edges if fp.FindOverlaps else []
        if len(verts) == 1:
            comp = verts[0]
        else:
            comp = Part.makeCompound(verts) if verts else Part.Shape()

        if hasattr(fp.ViewObject,"ShapeAppearance"):
            fp.ViewObject.ShapeAppearance = (FreeCAD.Material(DiffuseColor=fp.PointColor,), )
        else:
            fp.ViewObject.DiffuseColor = fp.PointColor
        fp.ViewObject.LineColor = fp.PointColor
        fp.ViewObject.LineWidth = fp.Sketch.ViewObject.LineWidth * 4
        fp.Shape = comp
        fp.CountIssues = count_issues
        if fp.CountIssues:
            fp.ViewObject.signalChangeIcon()
        report = f"""Total issues: {count_issues}
Open ended edges: {len(opens)-len(extra_ts)}
TSections: {len(tsections)+len(extra_ts)}
Cross Sections: {len(intersections)}
Overlapped edges: {overlaps}
Stray points: {len(strays)}"""
        fp.Report = report


class OpenEdgeFinderVP:
    def __init__(self, obj):
        obj.Proxy = self
        self.name = obj.Object.Name

    def attach(self, vobj):
        self.Object = vobj.Object

    def doubleClicked(self, vobj):
        fp = vobj.Object
        if fp.Sketch:
            fp.Sketch.ViewObject.doubleClicked()

    def setupContextMenu(self, vobj, menu):
        fp = vobj.Object
        if fp.Sketch:
            text = f"Edit {fp.Sketch.Label}"
            action = menu.addAction(text)
            action.triggered.connect(lambda: fp.Sketch.ViewObject.doubleClicked())
        text = "Print report to Report View"
        action = menu.addAction(text)
        action.triggered.connect(self.printReport)
    def printReport(self):
        fp = FreeCAD.ActiveDocument.getObject(self.name)
        msg = []
        msg.append(f"\n\nOpenEdgeFinder {fp.Version} report, now part of MeshRemodel workbench in Extras menu")
        msg.append(f"Sketch: {fp.Sketch.Label}\n{fp.Report}")
        if fp.OpenPoints:
            distances = []
            for p in fp.OpenPoints:
                dist_to_points = [math.fabs(p.distanceToPoint(p2)) for p2 in fp.OpenPoints if math.fabs(p.distanceToPoint(p2)) != 0.0]
                minimum = min(dist_to_points)
                if not minimum in distances:
                    distances.append(minimum)
            distances_strings = [f"{d} mm or {d*1e6} nanometers" if d < 1e4 else f"{d} mm" for d in distances]

            msg.append(f"Shortest distances between nearest open points:")
            for dis in distances_strings:
                msg.append(f"   {dis}")
        msg = "\n".join(msg)
        FreeCAD.Console.PrintMessage(msg)


    def dropObject(self, vp, dropped):
        fp = vp.Object
        fp.Sketch = dropped
        # for some reason dropped object gets kicked out of container, so put it back in
        if fp.Container:
            if not dropped in fp.Container.Group:
                fp.Container.Group += [dropped]

    def canDropObjects(self):
        return True

    def canDropObject(self, dropped):
        return hasattr(dropped, "ViewObject")

    def getIcon(self):
        fp = self.Object
        return __icon__ if not fp.CountIssues else __icon__.replace("None","red")

    def __getstate__(self):
        '''When saving the document this object gets stored using Python's json module.\
            If we have any unserializable stuff return them here or None'''
        return None

    def __setstate__(self,state):
        '''When restoring the serialized object from document we have the chance to set some internals here.\
                Since no data were serialized nothing needs to be done here.'''
        return None


def put_in_body(doc, fp, sk):
    """put object in same body or part as sketch"""
    bodies = [obj for obj in doc.Objects if obj.isDerivedFrom("PartDesign::Body") and sk in obj.Group]
    if bodies:
        bodies[0].Group = bodies[0].Group + [fp]
        fp.Container = bodies[0]
    else:
        parts = [obj for obj in doc.Objects if obj.isDerivedFrom("App::Part") and sk in obj.Group]
        if parts:
            parts[0].Group = parts[0].Group + [fp]
            fp.Container = parts[0]

__icon__ = """
/* XPM */
static char *dummy[]={
"64 64 5 1",
"# c #6b0000",
"a c #aa0000",
"b c #cacaca",
"c c limegreen",
". c None",
"................................................................",
"................................................................",
"................................................................",
"................................................................",
"................................................................",
"................................................................",
"................................................................",
"................................................................",
"......#####.....................................................",
".....#aaaaa#....................................................",
"....##aaaaaa#...................................................",
"...#aaaaaaaaa#bb................................................",
"...#aaa###aaa#bbbbbbbbbbbb......................................",
"...#aaa#b#aaa#bbbbbbbbbbbbbbbbbbbbbbb...........................",
"...#aaa###aaa#ccccccccbbbbbbbbbbbbbbbbbbbbbbbbb.................",
"...#aaaaaaaaa#bbbbbbbccccccccccccbbbbbbbbbbbbbbbbbbbb...........",
"....#aaaaaaa#...bbbbbbbbbbbbbbbbcccccccccccbbbbbbbbbbb..........",
".....#aaaaa#..............bbbbbbbbbbbbbbbbcccccccccccbb.........",
"......#####..........................bbbbbbbbbbbbbccbb..........",
"..............................................bbbccbbb..........",
"............................................bbbcccbbb...........",
"...........................................bbbccbbbb............",
"..........................................bbbccbbb..............",
".........................................bbbccbbb...............",
".......................................bbbcccbbb................",
"......................................bbbccbbbb.................",
".....................................bbbccbbb...................",
"....................................bbbccbbb....................",
"..................................bbbcccbbb.....................",
".................................bbbccbbbb......................",
"................................bbbccbbb........................",
"...............................bbbccbbb.........................",
".............................bbbcccbbb..........................",
"............................bbbccbbbb...........................",
"...........................bbbccbbb.............................",
"..........................bbbccbbb..............................",
"........................bbbcccbbb...............................",
".......................bbbccbbbb................................",
"......................bbbccbbb..................................",
".....................bbbccbbb...................................",
"...................bbbcccbbb....................................",
"..................bbbccbbbb.....................................",
".................bbbccbbb.......................................",
".................bbccbbbb.......................................",
"................bbccbbbbbbbbb...................................",
".................bbccccbbbbbbbbb................................",
".................bbbbbcccccbbbbbbbbb............................",
"....................bbbbbbcccccbbbbbbbbb........................",
".......................bbbbbbbcccccbbbbbbbbb.......#####........",
"...........................bbbbbbbccccbbbbbbbbbb..#aaaaa#.......",
"...............................bbbbbbcccccbbbbbbb##aaaaaa#......",
"...................................bbbbbbcccccbb#aaaaaaaaa#.....",
".......................................bbbbbbccc#aaa###aaa#.....",
"..........................................bbbbbb#aaa#b#aaa#.....",
"..............................................bb#aaa###aaa#.....",
"................................................#aaaaaaaaa#.....",
".................................................#aaaaaaa#......",
"..................................................#aaaaa#.......",
"...................................................#####........",
"................................................................",
"................................................................",
"................................................................",
"................................................................",
"................................................................"};
"""


class MeshRemodelOpenEdgeFinderCommandClass:
    
    def __init__(self):
        pass

    def GetResources(self):
        return {'Pixmap'  : os.path.join( iconPath , 'OpenEdgeFinder.svg') ,
            'MenuText': "Open edge finder" ,
            'ToolTip' : "Validate sketches and other objects."}
 
    def Activated(self):
    
        doc = FreeCAD.ActiveDocument
        if not doc:
            FreeCAD.Console.PrintError("No open document.  Open a document, select object to check, run macro.\n")
        else:

            sketches = Gui.Selection.getSelection()
            for sk in sketches:
                if not sk.isDerivedFrom("Part::Feature"):
                    continue
                fp = doc.addObject("Part::FeaturePython","OpenEdgeFinder")
                OpenEdgeFinder(fp)
                OpenEdgeFinderVP(fp.ViewObject)
                fp.Sketch = sk
                put_in_body(doc, fp, sk)
            doc.recompute()
   
    def IsActive(self):
        if not FreeCAD.ActiveDocument:
            return False
        sel = Gui.Selection.getSelection()
        if len(sel) == 0:
            return False
        if sel[0].isDerivedFrom("Part::Feature"):
            return True
        return False





##################################################################################################
class MeshRemodelGroupCommandPointsObjects:
    def GetCommands(self):
        return tuple([
        "MeshRemodelCreatePointsObject",
        "MeshRemodelCreateWireFrameObject",
        "MeshRemodelCreateCoplanarPointsObject",
        
        ])

    def GetDefaultCommand(self): # return the index of the tuple of the default command. This method is optional and when not implemented '0' is used  
        return 0

    def GetResources(self):
        return { 'MenuText': 'MeshRemodel mesh proxy objects', 'ToolTip': 'Proxy objects for dealing with meshes'}
        
    def IsActive(self): # optional
        return True


class MeshRemodelGroupCommandExtras:
    def GetCommands(self):
        
        return tuple([
                    "MeshRemodelDraftUpgrade",
                    "MeshRemodelCreateSketch",
                    "MeshRemodelMergeSketches",
                    "MeshRemodelFlattenDraftBSpline",
                    "MeshRemodelOpenEdgeFinder",
                    "MeshRemodelValidateSketch",
                    "MeshRemodelPartCheckGeometry",
                    "MeshRemodelSubShapeBinder",
                    "MeshRemodelSettings"]) # a tuple of command names that you want to group

    def GetDefaultCommand(self): # return the index of the tuple of the default command. This method is optional and when not implemented '0' is used  
        return 0

    def GetResources(self):
        return { 'MenuText': 'MeshRemodel Extras', 'ToolTip': 'Extra convenient commands'}
        
    def IsActive(self): # optional
        return True
        


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
        Gui.addCommand("MeshRemodelCreateTugBoat", MeshRemodelCreateTugBoatCommandClass())
        Gui.addCommand("MeshRemodelMoveAxial", MeshRemodelMoveAxialCommandClass())
        Gui.addCommand("MeshRemodelRotateObject", MeshRemodelRotateObjectCommandClass())
        Gui.addCommand("MeshRemodelGoBackSelection", MeshRemodelGoBackSelectionCommandClass())
        Gui.addCommand("MeshRemodelDraftUpgrade", MeshRemodelDraftUpgradeCommandClass())
        Gui.addCommand("MeshRemodelCreateSketch", MeshRemodelCreateSketchCommandClass())
        Gui.addCommand("MeshRemodelMergeSketches", MeshRemodelMergeSketchesCommandClass())
        Gui.addCommand("MeshRemodelValidateSketch", MeshRemodelValidateSketchCommandClass())
        Gui.addCommand("MeshRemodelOpenEdgeFinder",MeshRemodelOpenEdgeFinderCommandClass())
        Gui.addCommand("MeshRemodelPartCheckGeometry", MeshRemodelPartCheckGeometryCommandClass())
        Gui.addCommand("MeshRemodelSubShapeBinder", MeshRemodelSubShapeBinderCommandClass())
        Gui.addCommand("MeshRemodelSettings", MeshRemodelSettingsCommandClass())
        Gui.addCommand("MeshRemodelGroupCommandExtras",MeshRemodelGroupCommandExtras())
        Gui.addCommand("MeshRemodelGroupCommandPointsObjects", MeshRemodelGroupCommandPointsObjects())


initialize()
