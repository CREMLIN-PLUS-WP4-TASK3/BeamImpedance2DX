Geometry.OCCTargetUnit = "M";
Mesh.Algorithm=2; // 1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms
Mesh.Algorithm3D=1; // 1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT
Mesh.Format=1; // 1: msh, 2: unv, 10: auto, 16: vtk, 19: vrml, 21: mail, 26: pos stat, 27: stl, 28: p3d, 30: mesh, 31: bdf, 32: cgns, 33: med, 34: diff, 38: ir3, 39: inp, 40: ply2, 41: celum, 42: su2, 47: tochnog, 49: neu, 50: matlab
Mesh.RecombinationAlgorithm=1; // 0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad

r1 = 0.001;
r2 = 0.010;
r3 = 0.011;
r4 = 0.05;

DefineConstant[geomtype={0, Choices{0="sibc", 1="metal"}, Name"Geometry type"}];

Macro NewCircle
    // Input: newcircle_center - center point of a circle
    // Input: newcircle_r - radius of a circle
    // Output: newcircle_ll - line loop of a circle
    newcircle__center[] = Point{newcircle_center};
    newcircle__p1 = newp; Point(newcircle__p1) = {newcircle_r+newcircle__center[0],newcircle__center[1],0};
    newcircle__p2 = newp; Point(newcircle__p2) = {newcircle__center[0],newcircle_r+newcircle__center[1],0};
    newcircle__p3 = newp; Point(newcircle__p3) = {-newcircle_r+newcircle__center[0],newcircle__center[1],0};
    newcircle__p4 = newp; Point(newcircle__p4) = {newcircle__center[0],-newcircle_r+newcircle__center[1],0};
    newcircle__l1 = newl; Circle(newcircle__l1) = {newcircle__p1,newcircle_center,newcircle__p2};
    newcircle__l2 = newl; Circle(newcircle__l2) = {newcircle__p2,newcircle_center,newcircle__p3};
    newcircle__l3 = newl; Circle(newcircle__l3) = {newcircle__p3,newcircle_center,newcircle__p4};
    newcircle__l4 = newl; Circle(newcircle__l4) = {newcircle__p4,newcircle_center,newcircle__p1};
    newcircle_ll = newll; Line loop(newcircle_ll) = {newcircle__l1,newcircle__l2,newcircle__l3,newcircle__l4};
Return

Macro NewCircleTransfinite
    // Input: newcircletransfinite - transfinite value
    Transfinite Curve{newcircle__l1,newcircle__l2,newcircle__l3,newcircle__l4} = newcircletransfinite/4;
Return

newcircle_center = newp; Point(newcircle_center) = {0,0,0};

newcircle_r = r1;
newcircletransfinite = 160;
Call NewCircle;
Call NewCircleTransfinite;
ll1 = newcircle_ll;
s1 = news; Surface(s1) = {ll1};
Physical Surface(1) = {s1}; // beam

newcircle_r = r1+0.001;
newcircletransfinite = 50;
Call NewCircle;
Call NewCircleTransfinite;
ll1_2 = newcircle_ll;
s1_2 = news; Surface(s1_2) = {ll1_2,ll1};

newcircle_r = r2-0.001;
newcircletransfinite = 150;
Call NewCircle;
Call NewCircleTransfinite;
ll2_2 = newcircle_ll;
s2_2 = news; Surface(s2_2) = {ll2_2,ll1_2};

newcircle_r = r2;
newcircletransfinite = 1000;
Call NewCircle;
Call NewCircleTransfinite;
ll2 = newcircle_ll;
s2 = news; Surface(s2) = {ll2,ll2_2};
Physical Surface(2) = {s2,s1_2,s2_2}; // inner vacuum

If(geomtype)

newcircle_r = r3;
newcircletransfinite = 1000;
Call NewCircle;
Call NewCircleTransfinite;
ll3 = newcircle_ll;
s3 = news; Surface(s3) = {ll3,ll2};
Physical Surface(3) = {s3}; // metal

newcircle_r = r3+0.002;
newcircletransfinite = 150;
Call NewCircle;
Call NewCircleTransfinite;
ll3_2 = newcircle_ll;
s3_2 = news; Surface(s3_2) = {ll3_2,ll3};

newcircle_r = r4;
newcircletransfinite = 100;
Call NewCircle;
Call NewCircleTransfinite;
ll4 = newcircle_ll;
s4 = news; Surface(s4) = {ll4,ll3_2};
Physical Surface(4) = {s4,s3_2}; // outer vacuum

EndIf
