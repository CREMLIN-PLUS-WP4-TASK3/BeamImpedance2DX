Geometry.OCCTargetUnit = "M";
Mesh.Algorithm=2; // 1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms
Mesh.Algorithm3D=1; // 1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT
Mesh.Format=1; // 1: msh, 2: unv, 10: auto, 16: vtk, 19: vrml, 21: mail, 26: pos stat, 27: stl, 28: p3d, 30: mesh, 31: bdf, 32: cgns, 33: med, 34: diff, 38: ir3, 39: inp, 40: ply2, 41: celum, 42: su2, 47: tochnog, 49: neu, 50: matlab
Mesh.RecombinationAlgorithm=1; // 0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad

d0 = 0.0025;
d1 = 0.0089;
d2 = 0.01525;
d3 = 0.0165;
dx1 = 0;
dy1 = 0;
mesh1 = 0.01;
mesh2 = 0.01;
mesh3 = 0.01;

p0 = newp; Point(p0) = {0,0,0,mesh1};
p1 = newp; Point(p1) = {d0,dy1,0,mesh1};
p2 = newp; Point(p2) = {dx1,d0,0,mesh1};
p3 = newp; Point(p3) = {-d0,dy1,0,mesh1};
p4 = newp; Point(p4) = {dx1,-d0,0,mesh1};
l1 = newl; Circle(l1) = {p1,p0,p2};
l2 = newl; Circle(l2) = {p2,p0,p3};
l3 = newl; Circle(l3) = {p3,p0,p4};
l4 = newl; Circle(l4) = {p4,p0,p1};
ll1 = newll; Line loop(ll1) = {l1,l2,l3,l4};
s1 = news; Surface(s1) = {ll1};

p00 = newp; Point(p00) = {0,0,0,mesh1};
p10 = newp; Point(p10) = {d1,dy1,0,mesh1};
p20 = newp; Point(p20) = {dx1,d1,0,mesh1};
p30 = newp; Point(p30) = {-d1,dy1,0,mesh1};
p40 = newp; Point(p40) = {dx1,-d1,0,mesh1};
l10 = newl; Circle(l10) = {p10,p00,p20};
l20 = newl; Circle(l20) = {p20,p00,p30};
l30 = newl; Circle(l30) = {p30,p00,p40};
l40 = newl; Circle(l40) = {p40,p00,p10};
ll10 = newll; Line loop(ll10) = {l10,l20,l30,l40};
s10 = news; Surface(s10) = {ll10,ll1};

p01 = newp; Point(p01) = {0,0,0,mesh1};
p11 = newp; Point(p11) = {d2,dy1,0,mesh1};
p21 = newp; Point(p21) = {dx1,d2,0,mesh1};
p31 = newp; Point(p31) = {-d2,dy1,0,mesh1};
p41 = newp; Point(p41) = {dx1,-d2,0,mesh1};
l11 = newl; Circle(l11) = {p11,p01,p21};
l21 = newl; Circle(l21) = {p21,p01,p31};
l31 = newl; Circle(l31) = {p31,p01,p41};
l41 = newl; Circle(l41) = {p41,p01,p11};
ll11 = newll; Line loop(ll11) = {l11,l21,l31,l41};
s11 = news; Surface(s11) = {ll11,ll10};

p02 = newp; Point(p02) = {0,0,0,mesh1};
p12 = newp; Point(p12) = {d3,dy1,0,mesh1};
p22 = newp; Point(p22) = {dx1,d3,0,mesh1};
p32 = newp; Point(p32) = {-d3,dy1,0,mesh1};
p42 = newp; Point(p42) = {dx1,-d3,0,mesh1};
l12 = newl; Circle(l12) = {p12,p02,p22};
l22 = newl; Circle(l22) = {p22,p02,p32};
l32 = newl; Circle(l32) = {p32,p02,p42};
l42 = newl; Circle(l42) = {p42,p02,p12};
ll12 = newll; Line loop(ll12) = {l12,l22,l32,l42};
s12 = news; Surface(s12) = {ll12,ll11};

Transfinite Curve{l1,l2,l3,l4} = 10;
Transfinite Curve{l10,l20,l30,l40} = 70;
Transfinite Curve{l11,l21,l31,l41} = 70;
Transfinite Curve{l12,l22,l32,l42} = 40;

Physical Surface(1) = {s1};
Physical Surface(2) = {s10, s12};
Physical Surface(3) = {s11};
