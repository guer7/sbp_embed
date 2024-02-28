# Script to generate quadrilaterals on the complex part of the domain.
using GLMakie, HOHQMesh

h = 4

square_project = newProject("seabottom","seabottom_h" * string(h))

west = newEndPointsLineCurve("outer", [0, 0.0, 0.0], [0, -100.0, 0.0])
south = newEndPointsLineCurve("outer", [0, -100.0, 0.0], [100.0, -100.0, 0.0])
east = newEndPointsLineCurve("outer", [100.0, -100.0, 0.0], [100.0, 0.0, 0.0])
north = newEndPointsLineCurve("interface", [100.0, 0.0, 0.0], [0, 0.0, 0.0])

addCurveToOuterBoundary!(square_project,north)
addCurveToOuterBoundary!(square_project,west)
addCurveToOuterBoundary!(square_project,south)
addCurveToOuterBoundary!(square_project,east)

edge1 = newEndPointsLineCurve("inner", [45,-20, 0.0], [35, -40, 0.0])
edge2 = newEndPointsLineCurve("inner", [35, -40, 0.0], [55, -30, 0.0])
edge3 = newEndPointsLineCurve("inner", [55, -30, 0.0], [45, -20, 0.0])
addCurveToInnerBoundary!(square_project, edge1, "inner2")
addCurveToInnerBoundary!(square_project, edge2, "inner2")
addCurveToInnerBoundary!(square_project, edge3, "inner2")

circle1 = newCircularArcCurve("inner",       # curve name
                              [20, -25, 0.0], # circle center
                              10.0,             # circle radius
                              0.0,             # start angle
                              360.0,           # end angle
                              "degrees")       # angle units

addCurveToInnerBoundary!(square_project, circle1, "inner1")

spline_data = [ [0.0      18.0  -50.0 0.0]
                [0.2     13.0  -80.0 0.0]
                [0.4     38.0  -90.0 0.0]
                [0.6     78.0  -85.0 0.0]
                [0.8     68.0  -65.0 0.0]
                [1.0    18.0  -50.0 0.0] ]

spline = newSplineCurve("inner", 6, spline_data)
addCurveToInnerBoundary!(square_project, spline, "inner3")

edge1 = newEndPointsLineCurve("inner", [70.0, -30.0, 0.0], [65.0, -25.0, 0.0])
edge2 = newEndPointsLineCurve("inner", [65.0, -25.0, 0.0], [65.0, -15.0, 0.0])
edge3 = newEndPointsLineCurve("inner", [65.0, -15.0, 0.0], [70.0, -10.0, 0.0])
edge4 = newEndPointsLineCurve("inner", [70.0, -10.0, 0.0], [80.0, -10.0, 0.0])
edge5 = newEndPointsLineCurve("inner", [80.0, -10.0, 0.0], [85.0, -15.0, 0.0])
edge6 = newEndPointsLineCurve("inner", [85.0, -15.0, 0.0], [85.0, -25.0, 0.0])
edge7 = newEndPointsLineCurve("inner", [85.0, -25.0, 0.0], [80.0, -30.0, 0.0])
edge8 = newEndPointsLineCurve("inner", [80.0, -30.0, 0.0], [70.0, -30.0, 0.0])
addCurveToInnerBoundary!(square_project, edge1, "inner5")
addCurveToInnerBoundary!(square_project, edge2, "inner5")
addCurveToInnerBoundary!(square_project, edge3, "inner5")
addCurveToInnerBoundary!(square_project, edge4, "inner5")
addCurveToInnerBoundary!(square_project, edge5, "inner5")
addCurveToInnerBoundary!(square_project, edge6, "inner5")
addCurveToInnerBoundary!(square_project, edge7, "inner5")
addCurveToInnerBoundary!(square_project, edge8, "inner5")

r = newRefinementCenter("ref", "sharp", [45.0,-20.0,0.0], 3.0, 10.0)
addRefinementRegion!(square_project, r)

addBackgroundGrid!(square_project, [h, h, 0.0])

setPolynomialOrder!(square_project, 15)

plotProject!(square_project, MODEL+GRID)
@info "Press enter to generate the mesh and update the plot."
readline()

generate_mesh(square_project)
@info "Press enter to generate the mesh and update the plot."
readline()