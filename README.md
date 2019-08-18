# MeshRemodel Workbench
<img src="Resources/icons/MeshRemodelLogo.png" alt="icon">
<br/>
This workbench is currently under construction<br/>
<br/>

## Toolbar Icon
Download the <a href = "https://github.com/mwganson/MeshRemodel/blob/master/Resources/icons/MeshRemodelSVGLogo.svg">SVG Toolbar Icon</a><br/>

## Installation
Can be installed via the AddonManager in Tools menu -> AddonManager if you have a very recent build.  In the AddonManger click Configure, then add https://github.com/mwganson/MeshRemodel to the list of custom repositories.  After restarting the AddonManager you should find MeshRemodel in the list of workbenches you can install.
<br/>
## Overview
Use this workbench to aid in remodeling imported mesh objects.  The workflow is generally to first create a points object, then use those points in creating circles, arcs, lines, and polygons, which are then combined into wires, and then made into sketches.  Select 3 points to create a circle or arc, 2 points to create a line, 3 or more points to create a polygon / polyline.  These tools also work with any selectable points in the 3d view, not just points objects created with the workbench.<br/>
<br/>
## Create Points Object
<img src="Resources/icons/CreatePointsObject.png" alt="create points object"><br/>
Select the mesh object in the tree, then use this command to create a points object containing all the vertices of the selected mesh object.  The points object is a compound consisting of Part Point (vertex) objects, one per vertex in the selected mesh.  The purpose of this object is to provide selectable points in the 3d view.  We can use these selectable points with the other tools in the workbench to create the lines, circles, arcs, and polygons needed to remodel the mesh.<br/>
<br/>
## Create Line
<img src="Resources/icons/CreateLine.png" alt = "create line"><br/>
Select 2 points in the 3d view to enable this command.  It creates a Part Line object using the 2 selected points as a reference.  Note: any 2 selectable points may be used, even points that are part of an edge (the edge is shown selected in the 3d view, but the point at the mouse position is the one used), sketch, 3d object vertex, circle, etc.  If you can see it in the 3d view and select it, then it most likely can be used with this function.<br/>
<br/>
Use Ctrl+Click to include a point at the midpoint of the line.  Use Ctrl+Shift+Click if only the midpoint is desired.<br/>
In the report view you will find some basic information about the line, including its length and coordinates of its midpoint.<br/>
<br/>
## Create Polygon
<img src="Resources/icons/CreatePolygon.png" alt = "create polygon"><br/>
Select 3 or more points in the 3d view to enable this command.  It creates a polygon from the selected points.  Note: this is not a regular polygon, meaning the side lengths are not necessarily all equal to each other.  The order of selection is important.  By default the polygon will be closed, but you can prevent this with Shift+Click.  Use Ctrl+Click to include a point at the center of mass of the polygon.  Use Ctrl+Shift+Click of only the center of mass point is desired.  The points need not all lie on the same plane, but if they are not all on the same plane you will not be able to create a sketch from this later.<br/>
<br/>
In the report view you will find some basic information about the polygon, including its overall length and coordinates of its center of mass.<br/>
<br/>
## Create Circle
<img src="Resources/icons/CreateCircle.png" alt = "create circle"><br/>
Select 3 points in the 3d view to enable this command.  It creates a circle from those 3 points.  Use Ctrl+Click to add a point at the center of the new circle.  Use Ctrl+Shift+Click if only the center is desired.<br/>
<br/>
In the report view you will find some basic information about the circle, including its radius and coordinates of the center.<br/>
<br/>
## Create Arc
<img src="Resources/icons/CreateArc.png" alt = "create arc"><br/>
Select 3 points in the 3d view to enable this command.  It creates a Part Arc (internally using Part.ArcOfCircle() function) from those 3 selected points.  Use Ctrl+Click to include a point at the center of the arc.  Use Ctrl+Shift+Click if only the center is desired.<br/>
<br/>
In the report view you will find some basic information about the arc, including its radius and coordinates of its center.<br/>
<br/>
## Create Wire
<img src="Resources/icons/CreateWire.png" alt = "create wire"><br/>
Select 2 or more objects to enable this command.  It uses Draft.upgrade() to connect the objects into a single wire.  It is here as a convenience.  Note: the selected objects should all be connected together, but need not necessarily form a closed loop.  For example, you might have an arc and 2 lines connected one to each end of the arc.  You should not include circles unless you wish to connect them to other objects (not common).  The idea here to create wires from connected lines, open polygons, and arcs, then use these new wires, along with existing (coplanar) circles and closed polygons to create a sketch with the Create Sketch tool.  This is typically an intermediate step in creating a sketch.<br/>
<br/>
## Create Sketch
<img src="Resources/icons/CreateSketch.png" alt = "create sketch"><br/>
Select one or more wire objects to enable this command.  It uses Draft.makeSketch() to create a sketch from the selected objects.  It is here as a convenience.  If the selected objects are all coplanar this should usually result in a single sketch.  Some constraints are applied, but there will always be some constraining necessary.<br/>
<br/>
## Merge Sketches
<img src="Resources/icons/MergeSketches.png" alt = "merge sketches"><br/>
Select 2 or more sketches to enable this command.  This uses Sketcher workbench merge sketches command.  It is here as a convenience. 
<br/>

#### Release notes:<br/>
* 2019.08.18 (version 1.01)<br/>
** Add some simple information for created objects in the report view<br/>
* 2019.08.16 (version 1.00)<br/>
** Initial version
