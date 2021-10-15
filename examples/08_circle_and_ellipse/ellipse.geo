Geometry.OCCTargetUnit = "M";
Mesh.Algorithm=2; // 1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms
Mesh.Algorithm3D=1; // 1: Delaunay, 3: Initial mesh only, 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT
Mesh.Format=1; // 1: msh, 2: unv, 10: auto, 16: vtk, 19: vrml, 21: mail, 26: pos stat, 27: stl, 28: p3d, 30: mesh, 31: bdf, 32: cgns, 33: med, 34: diff, 38: ir3, 39: inp, 40: ply2, 41: celum, 42: su2, 47: tochnog, 49: neu, 50: matlab
Mesh.RecombinationAlgorithm=1; // 0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad

d1 = 0.01;
d2x = 0.04;
d2y = 0.03;
dx1 = 0;
dx2 = 0;
dy1 = 0;
dy2 = 0;
mesh1 = 0.0007;
mesh2 = 0.002;

p0 = newp; Point(p0) = {dx1,dy1,0,mesh1};
p1 = newp; Point(p1) = {d1+dx1,dy1,0,mesh1};
p2 = newp; Point(p2) = {dx1,d1+dy1,0,mesh1};
p3 = newp; Point(p3) = {-d1+dx1,dy1,0,mesh1};
p4 = newp; Point(p4) = {dx1,-d1+dy1,0,mesh1};
l1 = newl; Circle(l1) = {p1,p0,p2};
l2 = newl; Circle(l2) = {p2,p0,p3};
l3 = newl; Circle(l3) = {p3,p0,p4};
l4 = newl; Circle(l4) = {p4,p0,p1};
ll1 = newll; Line loop(ll1) = {l1,l2,l3,l4};
s1 = news; Surface(s1) = {ll1};


p00 = newp; Point(p00) = {dx2,dy2,0,mesh2};
p10 = newp; Point(p10) = {d2x+dx2,dy2,0,mesh2};
p20 = newp; Point(p20) = {dx2,d2y+dy2,0,mesh2};
p30 = newp; Point(p30) = {-d2x+dx2,dy2,0,mesh2};
p40 = newp; Point(p40) = {dx2,-d2y+dy2,0,mesh2};
l10 = newl; Ellipse(l10) = {p10,p00,p10,p20};
l20 = newl; Ellipse(l20) = {p20,p00,p10,p30};
l30 = newl; Ellipse(l30) = {p30,p00,p10,p40};
l40 = newl; Ellipse(l40) = {p40,p00,p10,p10};
ll10 = newll; Line loop(ll10) = {l10,l20,l30,l40};
s10 = news; Surface(s10) = {ll10,ll1};

Physical Surface(1) = {s1};
Physical Surface(2) = {s10};
