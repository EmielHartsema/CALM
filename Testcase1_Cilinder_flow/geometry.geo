// Gmsh project created on Thu Aug 27 12:47:27 2020
//SetFactory("OpenCASCADE");

// define parameters
diam = DefineNumber[ 0.5, Name "Parameters/diam" ]; // diameter of the circle
width = DefineNumber[ 2, Name "Parameters/width" ]; // width of the channel
length = DefineNumber[ 8, Name "Parameters/length" ]; // length of the channel
cos45 = DefineNumber[ 0.70710678118, Name "Parameters/cos45" ]; // cosine of 45 degrees

// add points cilinder
Point(1) = {width/2, 0, 0, 1.0};
Point(2) = {width/2+diam/2*cos45, diam/2*cos45, 0, 1.0};
Point(3) = {width/2+diam/2*cos45,-diam/2*cos45, 0, 1.0};
Point(4) = {width/2-diam/2*cos45,-diam/2*cos45, 0, 1.0};
Point(5) = {width/2-diam/2*cos45, diam/2*cos45, 0, 1.0};

// add points perimiter
Point(6) = {0,width/2,0,1.0};
Point(7) = {0,-width/2,0,1.0};
Point(8) = {length,width/2,0,1.0};
Point(9) = {length,-width/2,0,1.0};

// construct lines
// cilinder
Circle(1) = {5, 1, 2};
Circle(2) = {2, 1, 3};
Circle(3) = {3, 1, 4};
Circle(4) = {4, 1, 5};
//perimeter
Line(5) = {6, 7};
Line(6) = {6, 8};
Line(7) = {8, 9};
Line(8) = {9, 7};

// generate surface
Curve Loop(1) = {6, 7, 8, -5};
Curve Loop(2) = {1, 2, 3, 4};
Plane Surface(1) = {1, 2};

// specify boundary tags
Physical Curve("inflow") = {5};
Physical Curve("outflow") = {7};
Physical Curve("slip") = {6, 8};
Physical Curve("no_slip") = {4, 1, 2, 3};

//specify fluid domain
Physical Surface("fluid") = {1};

// specify boundaries between boundaryconditions
Physical Point("slip") = {6, 7, 8, 9};
Physical Point("no_slip") = {5, 2, 3, 4};
