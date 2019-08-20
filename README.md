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
## Create Coplanar Points Object
<img src="Resources/icons/CreateCoplanar.png" alt = "create coplanar"><br/>
Select 3 points from the points object in the 3d view to enable this command.  It creates a new points object filtered to contain only those points that are coplanar with the 3 selected points.  These are packed into a compound, which is then exploded in order to support (Shift+B) block selection.  The exploding is done as a separate document transaction, thus enabling you to undo this (Ctrl+Z) if you would prefer it to not be exploded.  You can also use Shift+Click to bypass the explode operation.<br/>
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
Select 3 or more points in the 3d view to enable this command.  It creates a polygon from the selected points.  Note: this is not a regular polygon, meaning the side lengths are not necessarily all equal to each other.  The order of selection is important.  By default the polygon will be closed, but you can prevent this with Shift+Click.<br/>
<br/>
When selecting using Shift+B, block selection, the points will generally need to be sorted or else you will get a polygon which zig zags all about.  To enable sorting using Alt+Click.  The sorting algorithm takes the first selected point, then finds the nearest point among the other selected points, and puts it 2nd.  Then it uses the 2nd point to find the next nearest point, and puts it 3rd, and so on.<br/>
<br/>
The polygon object created is a compound made up of individual Part Lines.  If you do not want the compound, but would instead prefer to have individual line objects, press Undo (Ctrl+Z) after making the polygon.  This will enable you to delete any lines you would prefer not to have, for example if you get a closed polygon, but would prefer it not to be closed or if some lines get crossed, etc.<br/>
<br/>
## Create BSpline
<img src="Resources/icons/CreateBSpline.png" alt = "create bspline"><br/>
Select 3 or more points in the 3d view to enable this command.  It creates a BSpline from the selected points.  The order of selection is important.  By default the BSpline will be closed, but you can prevent this with Shift+Click.  The points need not all lie on the same plane, but if they are not all on the same plane you will not be able to create a sketch from this later.<br/>
<br/>
This command supports block selections (Shift+B, draw rectangle).  Generally, the points will need to be sorted when using that block selection method.  Use Alt+Click to sort.  See the section on Create polygon for details on the sorting algorithm used.<br/>
<br/>
## Create Circle
<img src="Resources/icons/CreateCircle.png" alt = "create circle"><br/>
Select 3 (or more) points in the 3d view to enable this command.  It creates a circle from those first 3 selected points.  Any points after the first 3 are ignored, but you are allowed to select more.  This is to support easier block selection mode (Shift+B, draw rectangle).  Use Ctrl+Click to add a point at the center of the new circle.  Use Ctrl+Shift+Click if only the center is desired.<br/>
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
<br/>
## Validate Sketch
<img src="Resources/icons/ValidateSketch.png" alt = "validate sketch"><br/>
Opens Sketch workbench validate sketch tool.  Enabled only if you have 1 sketch selected.  It is here as a convenience.<br/>
<br/>
## Settings
<img src="Resources/icons/Settings.png" alt="settings"><br/>
### Keep toolbar active
This setting will keep the toolbar active after leaving the workbench, but you have to open the workbench at least once each session.  After that, if this is enabled, when you leave the workbench the toolbar will appear in the new workbench after a couple seconds.<br/>
<br/>
### Point size
This sets the point size on all points created with the workbench.  It does not affect objects already created after the setting is changed, only those created after the setting is changed.  Does not affect wire objects or sketch objects, only the centers of arcs, circles, and polygons, and the midpoints of lines.  Default 4.0<br/>
<br/>
### Line width
This sets the line width on all lines created with the workbench.  It does not affect objects already created after the setting is changed, only those created after the setting is changed.  Does not affect wire objects or sketch objects, only the edges of arcs, circles, lines, and polygons.  Default: 5.0<br/>
<br/>
### Sketch radius precision
This sets the precision to use when constraining radii (for circles and arcs) when creating sketches.  These are integer values from -1 to 12.  If -1, then no constraining of any radii occurs.  If 0, then radii are constrained to maximum precision.  If > 0, then radius constraints are rounded to that many digits precision, e.g. 1 results in 1.5, 2 in 1.49, 3 in 1.498, etc.  Additionally, if > 0, then equality constraints are used where possible (where radii round to the same value). Default: 1<br/>
<br/>
#### Release notes:<br/>
* 2019.08.20 (version 1.28)<br/>
** add sketch radius precision to settings
* 2019.08.20 (version 1.27)<br/>
** add shift+click option in creating coplanar points that are not exploded
* 2019.08.20 (version 1.26)<br/>
** fix coplanar point filter when original mesh was rotated
* 2019.08.20 (version 1.25)<br/>
** make sketch now constrains radii to precision = 0.1 mm
** try to place sketch even with first selected object center of mass
* 2019.08.20 (version 1.24)<br/>
** remove ? from title bar in settings dialog
* 2019.08.19 (version 1.23)<br/>
** show version info in settings dialog
* 2019.08.19 (version 1.22)<br/>
** convenience link to sketcher validate tool
* 2019.08.19 (version 1.21)<br/>
** fix for polygons needing coincidence constraints
* 2019.08.19 (version 1.20)<br/>
** add coplanar points object creation<br/>
** revert automatic selection of bspline and polygons for easier redo<br/>
** make polygons compound lines, explodable via undo (Ctrl+Z)<br/>
** support block select (Ctrl+B) for polygons, bsplines, and circles<br/>
* 2019.08.18 (version 1.10)<br/>
** add bspline creation
* 2019.08.18 (version 1.02)<br/>
** add settings for point size and line width
* 2019.08.18 (version 1.01)<br/>
** Add some simple information for created objects in the report view<br/>
* 2019.08.16 (version 1.00)<br/>
** Initial version
