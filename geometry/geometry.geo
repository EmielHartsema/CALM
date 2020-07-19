Merge "geometry.msh";
Line(5) = {6, 8};
//+
Line(6) = {6, 4};
//+
Line(7) = {8, 3};
//+
Line(8) = {4, 3};
//+
Circle(9) = {12, 11, 10};
//+
Circle(10) = {10, 11, 12};
//+
Curve Loop(1) = {6, 8, -7, -5};
//+
Curve Loop(2) = {10, 9};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("slipwall") = {6, 7};
//+
Physical Curve("Inflow") = {5};
//+
Physical Curve("outflow") = {8};
//+
Physical Curve("wall") = {10, 9};
//+
Physical Surface("domain") = {1};
