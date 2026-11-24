// A 2D flat-bottomed tank, half of it (symmetry plane at x = 0), for the
// celSetup userMesh example (WP G.3). Gmsh built-in kernel; 0.4 m half
// width, 0.6 m tall, outlet 0.06 m wide in the middle of the bottom, one
// cell thick (0.01 m) with front and back for an empty (2D) mesh.
SetFactory("Built-in");
Mesh.MshFileVersion = 2.2;
w = 0.4; h = 0.6; o = 0.03; t = 0.01;
Point(1) = {0, 0, 0}; Point(2) = {o, 0, 0}; Point(3) = {w, 0, 0};
Point(4) = {w, h, 0}; Point(5) = {o, h, 0}; Point(6) = {0, h, 0};
Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 5};
Line(5) = {5, 6}; Line(6) = {6, 1}; Line(7) = {2, 5};
// 5 mm across the outlet, growing to the wall; curves 4 and 6 run the
// other way, so their progression is inverted
Transfinite Curve {1, 5} = 7;
Transfinite Curve {2} = 31 Using Progression 1.06;
Transfinite Curve {4} = 31 Using Progression 1 / 1.06;
Transfinite Curve {3, 7} = 31 Using Progression 0.98;
Transfinite Curve {6} = 31 Using Progression 1 / 0.98;
Curve Loop(1) = {1, 7, 5, 6}; Plane Surface(1) = {1};
Curve Loop(2) = {2, 3, 4, -7}; Plane Surface(2) = {2};
Transfinite Surface {1, 2}; Recombine Surface {1, 2};
ex1[] = Extrude {0, 0, t} { Surface{1}; Layers{1}; Recombine; };
ex2[] = Extrude {0, 0, t} { Surface{2}; Layers{1}; Recombine; };
// ex[]: top, volume, sides in curve-loop order
Physical Surface("outlet") = {ex1[2]};
Physical Surface("bottom") = {ex2[2]};
Physical Surface("wall") = {ex2[3]};
Physical Surface("top") = {ex1[4], ex2[4]};
Physical Surface("axis") = {ex1[5]};
Physical Surface("frontAndBack") = {1, 2, ex1[0], ex2[0]};
Physical Volume("internal") = {ex1[1], ex2[1]};
