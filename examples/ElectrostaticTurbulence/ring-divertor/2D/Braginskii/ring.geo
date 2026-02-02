xc = 2.0;
yc = 0.0;
r0 = 0.5;
r1 = 1.0;

xd = 0.25;
ncirc = 1.65*r0/0.01;
nrad = (r1-r0)/0.01;

Point(1) = {xc,    yc,    0};
Point(2) = {xc+r0, yc,    0};
Point(3) = {xc,    yc+r0, 0};
Point(4) = {xc-r0, yc,    0};
Point(5) = {xc,    yc-r0, 0};

Point(6) = {xc+r1, yc,    0};
Point(7) = {xc,    yc+r1, 0};
Point(8) = {xc-r1, yc,    0};
Point(9) = {xc,    yc-r1, 0};
Point(10)= {xc,    yc,    0};


Point(11) = {xc-0.707106781*r1, yc-0.707106781*r1, 0};
Point(12) = {xc+0.707106781*r1, yc-0.707106781*r1, 0};
Point(13) = {xc-0.707106781*(r1-xd), yc-0.707106781*(r1-xd), 0};
Point(14) = {xc+0.707106781*(r1-xd), yc-0.707106781*(r1-xd), 0};
//Point(11) = {xc-0.707106781*r1, yc-0.707106781*r1, 0};
//Point(12) = {xc+0.707106781*r1, yc-0.707106781*r1, 0};
//Point(13) = {xc-0.707106781*r1+xd, yc-0.707106781*r1, 0};
//Point(14) = {xc+0.707106781*r1-xd, yc-0.707106781*r1, 0};
//Point(15) = {xc, yc-r1+xd*0.707106781, 0};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Circle(5) = {6,10,7};
Circle(6) = {7,10,8};
Circle(7) = {8,10,11};
Circle(8) = {12,10,6};
Circle(9) = {13, 10, 14};
//Circle(10) = {14, 10, 15};

Transfinite Curve{1:10} = ncirc;

Line(11) = {2, 6};
Line(12) = {3, 7};
Line(13) = {4, 8};
//Line(14) = {5, 15};
Transfinite Curve{11:13} = nrad Using Progression 1;

Line(15) = {11, 13};
Line(16) = {12, 14};
Transfinite Curve{15:16} = nrad Using Progression 1;


Curve Loop(11) = {11, 5, -12, -1};
Plane Surface(12) = {11};
Transfinite Surface{12};
Recombine Surface{12};

Curve Loop(13) = {12, 6, -13, -2};
Plane Surface(14) = {13};
Transfinite Surface{14};
Recombine Surface{14};

Curve Loop(15) = {13, 7, 15, 9, -16, 8, -11, -4, -3};
Plane Surface(16) = {15};
//Transfinite Surface{16};
Recombine Surface{16};

//Curve Loop(17) = {14, -10, -16, 8, -11, -4};
//Plane Surface(18) = {17};
//Transfinite Surface{18};
//Recombine Surface{18};

Physical Surface(0) = {12,14,16};
Physical Line(1) = {1,2,3,4};
Physical Line(2) = {5,6,7,8};
Physical Line(3) = {9};
Physical Line(4) = {15};
Physical Line(5) = {16};//+
Physical Curve(17) = {6, 5, 8, 7};
//+
Physical Curve(18) = {2, 1, 4, 3};
//+
Physical Curve(19) = {9};
//+
Physical Curve(20) = {15};
//+
Physical Curve(21) = {16};
