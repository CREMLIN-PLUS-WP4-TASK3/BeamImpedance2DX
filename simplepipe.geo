lveryfine=2e-3;
lfine = 1e-2;
lrough = 5;

a=0.01;
b=0.04;


Point(1) = {0, 0, 0, lveryfine};

Point(2) = {a, 0, 0, lveryfine};
Point(3) = {b,0, 0, lfine};

Point(4) = {0,a, 0, lveryfine};
Point(5) = {0,b, 0, lfine};

Point(6) = {-a,0, 0, lveryfine};
Point(7) = {-b,0, 0, lfine};

Point(8) = {0,-a, 0, lveryfine};
Point(9) = {0,-b, 0, lfine};



Line(1) = {1,2};
Line(2) = {2,3};

Line(3) = {1,4};
Line(4) = {4,5};

Line(5) = {7,6};
Line(6) = {6,1};

Line(7) = {9,8};
Line(8) = {8,1};



Circle(9) = {3,1,5};
Circle(10) = {5,1,7};
Circle(11) = {7,1,9};
Circle(12) = {9,1,3};

Circle(13) = {2,1,4};
Circle(14) = {4,1,6};
Circle(15) = {6,1,8};
Circle(16) = {8,1,2};

Line Loop(30)={1,13,-3};
Line Loop(31)={3,14,6};
Line Loop(32)={-6,15,8};
Line Loop(33)={-8,16,-1};


//Line Loop(30)={13,14,15,16};


Line Loop(34)={9,10,11,12};

Plane Surface(1) = {30};
Plane Surface(2) = {31};
Plane Surface(3) = {32};
Plane Surface(4) = {33};
Plane Surface(5) = {34,30,31,32,33};

//Line {1} In Surface {30};


Physical Surface(1) = {1,2,3,4};   //beam
//Physical Surface(2) = {2};
//Physical Surface(3) = {3};
//Physical Surface(4) = {4};
Physical Surface(2) = {5};   //vacuum

//Physical Line(6)={1};


