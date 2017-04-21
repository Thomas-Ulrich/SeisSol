lc= 5000;
lc_fault = 200;

xlim=70000;
ylim=70000;

zlayerArray_1=0;
zlayerArray_2=-1000;
zlayerArray_3=-1800;
zlayerArray_4=-2000;
zlayerArray_5=-2100;
zlayerArray_6=-3000;
zlayerArray_7=-3400;
zlayerArray_8=-3500;
zlayerArray_9=-3900;
zlayerArray_10=-5800;
zlayerArray_11=-8300;
zlayerArray_12=-12700;
zlayerArray_13=-14100;
zlayerArray_14=-15000;
zlayerArray_15=-17100;
zlayerArray_16=-17500;
zlayerArray_17=-20300;
zlayerArray_18=-20400;
zlayerArray_19=-xlim;

nLayers = 19;
//id fault layer (14)-2
nLayersFault=12;
i0=0;

//horizontal layers
For i In {1:nLayers}
   zlayer = zlayerArray~{i};
   Point(1+i0) = {-xlim, -ylim, zlayer, lc};
   Point(2+i0) = { xlim,   -ylim, zlayer, lc};
   Point(3+i0) = { xlim,   0, zlayer, lc};
   Point(4+i0) = {-xlim,  0, zlayer, lc};
   If (zlayer>=-15e3)
      Point(5+i0) = {10e3,  0, zlayer, lc};
      Point(6+i0) = {-30e3,   0, zlayer, lc};
   EndIf
   Line(1+i0) ={1+i0,2+i0};
   Line(2+i0) ={2+i0,3+i0};

   If (zlayer<-15e3)
      Line(3+i0) ={3+i0,4+i0};
      Line(4+i0) ={4+i0,1+i0};
      Line Loop(i0) = {1+i0, 2+i0, 3+i0, 4+i0};
   EndIf
   If (zlayer>=-15e3)
      Line(3+i0) ={3+i0,5+i0};
      Line(4+i0) ={5+i0,6+i0};
      Line(5+i0) ={6+i0,4+i0};
      Line(6+i0) ={4+i0,1+i0};
      Line Loop(i0) = {1+i0, 2+i0, 3+i0, 4+i0,5+i0,6+i0};
   EndIf
   Plane Surface(i0) = {i0};
   i0=i0+6;
EndFor



//vertical sides
v0=10000;
i0=0;
For i In {2:nLayers}
   zlayer = zlayerArray~{i};
   Line(1+i0+v0) ={1+i0,1+i0+6};
   Line(2+i0+v0) ={2+i0,2+i0+6};
   Line(3+i0+v0) ={3+i0,3+i0+6};
   Line(4+i0+v0) ={4+i0,4+i0+6};
   If (zlayer>=-15e3)
      Line(5+i0+v0) ={5+i0,5+i0+6};
      Line(6+i0+v0) ={6+i0,6+i0+6};
   EndIf
   Line Loop(i0+1) = {1+i0, 2+i0+v0, -(1+i0+6), -(1+i0+v0)};
   Plane Surface(i0+1) = {i0+1};
   Line Loop(i0+2) = {2+i0, 3+i0+v0, -(2+i0+6), -(2+i0+v0)};
   Plane Surface(i0+2) = {i0+2};
   
   If (zlayer<-15e3)
      If (Fabs(zlayer+17100)<1e-3)
         Line Loop(i0+3) = {3+i0, 4+i0, 5+i0, 4+i0+v0, -(3+i0+6), -(3+i0+v0)};
         Plane Surface(i0+3) = {i0+3};
         Line Loop(i0+4) = {6+i0, 1+i0+v0, -(4+i0+6), -(4+i0+v0)};
         Plane Surface(i0+4) = {i0+4};
      EndIf
      If (Fabs(zlayer+17100)>=1e-3)
         Line Loop(i0+3) = {3+i0, 4+i0+v0, -(3+i0+6), -(3+i0+v0)};
         Plane Surface(i0+3) = {i0+3};
         Line Loop(i0+4) = {4+i0, 1+i0+v0, -(4+i0+6), -(4+i0+v0)};
         Plane Surface(i0+4) = {i0+4};
      EndIf
   EndIf
   
   If (zlayer>=-15e3)
      Line Loop(i0+3) = {3+i0, 5+i0+v0, -(3+i0+6), -(3+i0+v0)};
      Plane Surface(i0+3) = {i0+3};
      Line Loop(i0+4) = {4+i0, 6+i0+v0, -(4+i0+6), -(5+i0+v0)};
      Plane Surface(i0+4) = {i0+4};
      Line Loop(i0+5) = {5+i0, 4+i0+v0, -(5+i0+6), -(6+i0+v0)};
      Plane Surface(i0+5) = {i0+5};
      Line Loop(i0+v0) = {6+i0, 1+i0+v0, -(6+i0+6), -(4+i0+v0)};
      Plane Surface(i0+v0) = {i0+v0};
   EndIf
   i0=i0+6;
EndFor


//volumes
v0=10000;
i0=0;
For i In {2:nLayers}
   zlayer = zlayerArray~{i};

   If (zlayer<-15e3)
      Surface Loop(i0+v0+1) = {i0+1, i0+2, i0+3, i0+4, i0+6, i0};
      Volume(i0/6+1) = {i0+v0+1};
   EndIf

   If (zlayer>=-15e3)
      Surface Loop(i0+v0+1) = {i0+1, i0+2, i0+3, i0+4, i0+5, i0+v0, i0+6, i0};
      Volume(i0/6+1) = {i0+v0+1};
   EndIf
   Physical Volume(i0/6+1) = {i0/6+1};

   i0=i0+6;
EndFor

/*
Physical Surface(101) = {0};
Physical Surface(103) = {4:4+nLayersFault*6:6};
//wrong
Physical Surface(105) = {1:1+nLayersFault*6:6,2:2+nLayersFault*6:6,3:3+nLayersFault*6:6,5:5+nLayersFault*6:6,10000:10000+nLayersFault*6:6,
7+nLayersFault*6:103:6,8+nLayersFault*6:104:6,9+nLayersFault*6:105:6,10+nLayersFault*6:106:6,108};
*/

//Defining Ruled surface for defining mesh coarsening
Line(10107) = {83,5};
Line(10108) = {6,84};
Line Loop(10109)={4,10108,-82,10107};
Ruled Surface(10110) = {10109};

/*
Field[1] = Attractor;
Field[1].FacesList = {10110};
Field[1].FieldX = -1;
Field[1].FieldY = -1;
Field[1].FieldZ = -1;
Field[1].NNodesByEdge = 20;
Field[2] = MathEval;
Field[2].F = Sprintf("0.05*F1 +(F1/2.5e3)^2 + %g", lc_fault);
//Field[2].F = Sprintf("0.1*F1 +(F1/5.0e3)^2 + %g", lc_fault);
Background Field = 2;
*/
//Physical Surface(103) = {4:4+nLayersFault*6:6};
Physical Surface(10111) = {0,1,2,3,5,6,7,8,9,11,12,13,14,15,17,18,19,20,21,23,24,25,26,27,29,30,31,32,33,35,36,37,38,39,41,42,43,44,45,47,48,49,50,51,53,54,55,56,57,59,60,61,62,63,65,66,67,68,69,71,72,73,74,75,77,78,79,80,81,82,84,85,86,87,88,90,91,92,93,94,96,97,98,99,100,102,103,104,105,106,108,10000,10006,10012,10018,10024,10030,10036,10042,10048,10054,10060,10066,10072};

//{2, 8, 14, 20, 26, 32, 38, 44, 50, 56, 62, 68, 74, 80, 86, 92, 98, 104, 3, 9, 15, 21, 27, 33, 39, 45, 51, 57, 63, 69, 75, 1, 7, 13, 19, 25, 4, 31, 37, 43, 10, 16, 22, 49, 28, 34, 40, 46, 55, 52, 58, 61, 64, 67, 73, 70, 79, 76, 85, 91, 97, 0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 103, 84, 90, 96, 102, 10000, 5, 10006, 11, 10012, 17, 10018, 23, 10024, 29, 10030, 35, 10036, 41, 10042, 47, 10048, 53, 10054, 59, 10060, 65, 10066, 71, 10072, 77, 81, 82, 87, 88, 93, 94, 99, 100, 108, 105, 106};
//Volume(10112) = {10073};
