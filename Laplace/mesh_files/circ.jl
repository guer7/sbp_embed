# Script to generate quadrilaterals of varying size on the unit circle.
using GLMakie, HOHQMesh

hvec = [0.5,0.4,0.3,0.2,0.1,0.08,0.06]

for h in hvec

	circ_project = newProject("circle","circ_h" * string(h))

	xEqn = "x(t) = 1*cos(2*pi*t)"
	yEqn = "y(t) = 1*sin(2*pi*t)"
	zEqn = "z(t) = 0.0"

	circ = newParametricEquationCurve("outer", xEqn, yEqn, zEqn)

	addCurveToOuterBoundary!(circ_project,circ)

	addBackgroundGrid!(circ_project, [2.0, 2.0, 0.0])
	setBackgroundGridSize!(circ_project, h, h)

	setPolynomialOrder!(circ_project, 15)
	generate_mesh(circ_project)
end