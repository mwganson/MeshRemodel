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
Use this workbench to aid in remodeling imported mesh objects.  The preferred workflow is to select the mesh, then click create points object.  This creates an object with selectable points at all of the mesh vertices.  The next step is to choose a face you wish to remodel.  Select 3 points on that face, then click the create coplanar points tool.  This creates an empty sketch mapped using 3 points make a plane mode, and into which all of the coplanar points are links to external geometry.  Use the sketcher tools, such as circle from 3 points or arc from 3 points, etc. to remodel the profile as a sketch.  If desired, it is possible to delete all of the links to external geometry later using the sketcher validation tool.  You would then be able to rearrange the elements relative to their rotation and/or to the sketch's origin, if desired.<br/>
<br/>
An alternative workflow is to first create a points object, then use those points in creating circles, arcs, lines, and polygons, which are then combined into wires, and then made into sketches.  Select 3 points to create a circle or arc, 2 points to create a line, 3 or more points to create a polygon / polyline.  These tools also work with any selectable points in the 3d view, not just points objects created with the workbench.  The problem with this is sometimes the objects created are no coplanar, and so you end up with problems in future operations, such as extrude or trying to make them into a sketch.<br/>
<br/>
## Create Points Object
<img src="Resources/icons/CreatePointsObject.png" alt="create points object"><br/>
Select the mesh object in the tree, then use this command to create a points object containing all the vertices of the selected mesh object.  The points object is a compound consisting of Part Point (vertex) objects, one per vertex in the selected mesh.  The purpose of this object is to provide selectable points in the 3d view.  We can use these selectable points with the other tools in the workbench to create the lines, circles, arcs, and polygons needed to remodel the mesh.<br/>
<br/>
## Create Point Object
<img src="Resources/icons/CreatePointObject.png" alt="create point object"><br/>
Select a vertex in the 3d view, then use this command to create a point object containing that vertex.  The point object is a Part::Vertex that we can use in some operations, such as Part::Loft.<br/>
<br/>
## Create Coplanar Points Object
<img src="Resources/icons/CreateCoplanar.png" alt = "create coplanar"><br/>
Select 3 points from the points object in the 3d view to enable this command.  It creates a new points object filtered to contain only those points that are coplanar with the 3 selected points.  Then an empty sketch is created, and added to that empty sketch are links to external geometry for all of the points in the new coplanar points object.  It is advised to recreate the profile inside the sketch using those external links and the sketcher tools.<br/>
<br/>
In order to filter the original points object into a set of coplanar points aligned on the plane defined by the 3 selected points an internal isCoplanar algorithm is used.  An alternative (Alt+Click) is to use the isPlaner() algorithm from the Draft workbench.  The difference is the Draft workbench algorithm is more restrictive, and so you end up with fewer points than when using the default internal algorithm at its default tolerance level.  There is a settings option for changing the tolerance level.  The smaller the number the fewer points get produced.  The filtering is done by using the 3 selected points and each other point in turn to create a tetrahedron.  If the 4 points are coplanar, then the tetrahedron should have volume ~= zero.  Default tolerance is 0.001 mm^3.
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
It is recommended to remodel the object inside the sketcher using the sketcher line tools instead of this, but it is here for those who wish to use it instead.  The reason it is recommened to use the sketcher is some of the polygons created might not be coplanar, and might give problems with future operations.<br/>
<br/>
Select 3 or more points in the 3d view to enable this command.  It creates a polygon from the selected points.  Note: this is not a regular polygon, meaning the side lengths are not necessarily all equal to each other.  The order of selection is important.  By default the polygon will be closed, but you can prevent this with Shift+Click.<br/>
<br/>
When selecting using Shift+B, block selection, the points will generally need to be sorted or else you will get a polygon which zig zags all about.  To enable sorting using Alt+Click.  The sorting algorithm takes the first selected point, then finds the nearest point among the other selected points, and puts it 2nd.  Then it uses the 2nd point to find the next nearest point, and puts it 3rd, and so on.<br/>
<br/>
The polygon object created is made up of individual Part Lines.  This will enable you to delete any lines you would prefer not to have, for example if you get a closed polygon, but would prefer it not to be closed or if some lines get crossed, etc.  Use the Create wire tool to upgrade the individual lines to a single wire object, if you prefer.<br/>
<br/>
## Create BSpline
<img src="Resources/icons/CreateBSpline.png" alt = "create bspline"><br/>
Select 3 or more points in the 3d view to enable this command.  It creates a BSpline from the selected points.  The order of selection is important.  By default the BSpline will be closed, but you can prevent this with Shift+Click.  The points need not all lie on the same plane, but if they are not all on the same plane you will not be able to create a sketch from this later.  Sometimes points that appear to lie on the same plane are not actually on the same plane.  It is better to create the bspline in the sketcher.<br/>
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
Creates a sketch, optionally attached to 3 points on a plane if 3 points are selected.  This does not create any links to external geometry.  See Create coplanar points command if you want to automatically import all coplanar points that lie on this same plane.<br/>
<br/>
Use Ctrl+Click to make a sketch out of selected circles, polygons, etc.  If a circle or arc is the first selected object, it will map the sketch concentrically to that circle or arc.  Note: there is a known issue using this method that sometimes objects that appear to be coplanar might not actually be coplanar.  It is recommended to remodel using the sketcher with links to external geometry to the points objects instead of this method. Uses method of creating a single sketch from all selected objects.<br/>
<br/>
Use Alt+Click to create multiple sketches, one from each object selected, and then merge them all together into a single sketch, deleting the temporary sketches afterward.  This can sometimes resolve coplanar issues.<br/>
<br/>
## Merge Sketches
<img src="Resources/icons/MergeSketches.png" alt = "merge sketches"><br/>
Select 2 or more sketches to enable this command.  This uses Sketcher workbench merge sketches command.  It is here as a convenience. 
<br/>
<br/>
## Validate Sketch
<img src="Resources/icons/ValidateSketch.png" alt = "validate sketch"><br/>
Opens Sketch workbench validate sketch tool.  Enabled only if you have 1 sketch selected.  It is here as a convenience.  Occasionally, sketches will have missing coincidence constraints.  That tool is good for fixing that issue.  It can also be used to easily remove all links to external geometry.<br/>
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
This sets the precision to use when constraining radii (for circles and arcs) when creating sketches.  These are integer values from -1 to 12.  If -1, then no constraining of any radii occurs.  If 0, then radii are constrained to maximum precision.  If > 0, then radius constraints are rounded to that many digits precision, e.g. 1 results in 1.5, 2 in 1.49, 3 in 1.498, etc. Default: 1<br/>
<br/>
### Coplanar tolerance
This sets the tolerance to use when determining which points lie on the same plane as the 3 selected points that define the plane.  Higher numbers mean less restrictive results, producing more points, not all of which might be accepted as coplanar in later operations.  This is not an issue when modeling within the sketcher using links to external geometry.  It is recommened to not change the default unless you are missing some points that you think should be included or perhaps you are getting points that should not be included.  The tolerance number represents the volume of a tetrahedron created using the 3 selected points and the point currently under consideration in cubic mm.  Default: 0.001 mm^3
#### Release notes:<br/>
* 2020.08.06 (version 1.43)<br/>
** added create a point command
* 2020.08.05 (version 1.42)<br/>
** add downgrade option to upgrade (shift+click to downgrade)
* 2020.08.05 (version 1.41)<br/>
** select the lines and arcs used to create a wire after wire creation (updgrade) so it is easier to delete them
** select the lines after creating a polygon so it is easier to upgrade them to a wire
* 2020.08.05 (version 1.4)<br/>
** improved settings dialog to show current settings
** no longer using Draft coplanar check
** no longer creating sketch by default when creating coplanar points object
** default coplanar tolerance now 0.0001 mm (was .001 mm)
* 2019.08.28 (version 1.31)<br/>
** open sketcher workbench when needed to prevent errors
* 2019.08.22 (version 1.30)<br/>
** Lots of changes here.  We've moved from creating objects in the 3d view to be used later for creating a sketch to creating an empty sketch with links to external geometry automatically added.  Elements are then created directly in the sketch using sketcher tools.  This is to avoid some coplanar issues that can arise out of creating elements in the 3d view with the tools in this workbench.<br/>
* 2019.08.20 (version 1.292)<br/>
** option to attach sketch concentrically to first selected object, if that is a circle or arc with Shift+Click.
** bug fix where sketch is sometimes reported as noncoplanar, but lost ability to place equality constraints on some circles/arcs
* 2019.08.20 (version 1.291)<br/>
** fix bug in calls to gu.isColinear()
* 2019.08.20 (version 1.29)<br/>
** reorganize code internally, consolidating geometry utilities into a single class<br/>
** accessible via python:<br/>
**  from MeshRemodelCmd import MeshRemodelGeomUtils<br/>
**  gu = MeshRemodelGeomUtils()<br/>
**  gu.circumcenter(A,B,C)<br/>
**  gu.circumradius(A,B,C)<br/>
**  gu.inradius(A,B,C)<br/>
**  gu.incenter(A,B,C)<br/>
**  gu.dist(A,B)<br/>
**  gu.getDistance3d(Ax,Ay,Az,Bx,By,Bz)<br/>
**  gu.isCoplanar(A,B,C)<br/>
**  gu.isColinear(A,B,C)<br/>
**  gu.nearestPoint(pt, pts, exclude)<br/>
**  gu.sortPoints(pts)<br/>
**  gu.midpoint(A,B)<br/>
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
