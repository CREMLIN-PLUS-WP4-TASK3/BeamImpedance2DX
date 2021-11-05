a=0.001;
mesh_beam = a/15;
mesh = a/4;

DefineConstant[geomtype={0, Choices{0="sibc", 1="metal"}, Name"Geometry type"}];

Point(1) = {0,0,0,mesh_beam};
Point(2) = {a,0,0,mesh_beam};
Point(3) = {0,a,0,mesh_beam};
Point(4) = {-a,0,0,mesh_beam};
Point(5) = {0,-a,0,mesh_beam};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Point(30) = {-0.003,0.03,0,mesh};
Point(31) = {-0.05,0.03,0,mesh};
Point(32) = {-0.05,0.05,0,mesh};
Point(33) = {0.05,0.05,0,mesh};
Point(34) = {0.05,0.03,0,mesh};
Point(35) = {0.003,0.03,0,mesh};
Point(36) = {0.003,-0.03,0,mesh};
Point(37) = {0.05,-0.03,0,mesh};
Point(38) = {0.05,-0.05,0,mesh};
Point(39) = {-0.05,-0.05,0,mesh};
Point(40) = {-0.05,-0.03,0,mesh};
Point(41) = {-0.003,-0.03,0,mesh};

Line(30) = {30,31};
Line(31) = {31,32};
Line(32) = {32,33};
Line(33) = {33,34};
Line(34) = {34,35};
Line(35) = {35,36};
Line(36) = {36,37};
Line(37) = {37,38};
Line(38) = {38,39};
Line(39) = {39,40};
Line(40) = {40,41};
Line(41) = {41,30};

Line Loop(30) = {30,31,32,33,34,35,36,37,38,39,40,41};
Plane Surface(30) = {30,1};

Physical Surface(1) = {1}; // beam
Physical Surface(2) = {30}; // vacuum

If(geomtype)
Line(50) = {31,40};
Line(51) = {37,34};

Line Loop(50) = {40,41,30,50};
Line Loop(51) = {34,35,36,51};
Plane Surface(50) = {50};
Plane Surface(51) = {51};

Physical Surface(3) = {50,51}; // metal
EndIf
