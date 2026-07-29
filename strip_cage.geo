
// =====================================================================
// strip_cage.geo  --  graded-boundary version
//
// FIXES IN THIS VERSION
//   A/B bugfixes:
//     * lc_arc was defined on line 7 using strip_t, which was not
//       declared until line 12  ->  Gmsh error, lc_arc evaluated to 0.
//       All size parameters are now declared before first use.
//     * "Strip_w" (capital S) on the n_flat line was an undefined
//       variable -> n_flat evaluated to 1 -> degenerate Transfinite
//       curve -> parser abort and a collapsed 2822-element mesh.
//     * The n_flat fix had only been applied to the LEFT strip; the
//       right strip was still hardcoded to 4.  That alone guaranteed a
//       left/right asymmetric mesh.  Both sides are now identical.
//
//   C - graded transition (the actual fix):
//     Each flat edge is split at its midpoint into two half-edges that
//     are graded with a geometric Progression.  Edge length starts at
//     h_junc (a gentle 2:1 step off the cap arcs) and grows to h_max at
//     the flat midpoint.  This removes the ~100-170x edge-length step
//     that previously existed where the flat met the cap arc.
// =====================================================================

// ---------- 1. Geometry ----------
pitch    = 0.254;
n_strips = 118;
gap      = 3.0;
strip_w  = 0.20;     // 2 mm strip width  (long dimension, along y)
strip_t  = 0.0025;   // 25 um strip thickness -> 12.5 um cap radius
z_thick  = 0.01;
y_offset = 14.859;

// ---------- 2. Mesh size control ----------
// Everything below is DERIVED from the geometry above, so the mesh
// stays consistent if strip_w / strip_t / pitch are ever changed.

n_arc  = 6;                                // nodes per 90 deg cap arc (matches field_cage)
r_cap  = strip_t/2.0;                      // 12.5 um cap radius
h_arc  = (Pi/2.0)*r_cap/(n_arc - 1.0);     // resulting edge length on the caps (~3.93 um)

// h_junc = size of the FIRST flat element, expressed as a multiple of h_arc.
//          This multiple IS the residual step at the arc/flat junction.
//          It was 167x before; 6x here.  Lower = better quality, more nodes.
h_junc = 6.0*h_arc;
gap_v  = pitch - strip_w;                  // 540 um vertical gas gap between adjacent strips
n_gap  = 2;                                // elements demanded across that gap
h_max  = gap_v/n_gap;                      // coarsest flat element, at the flat midpoint

// ---- Measured trade-off (full 118 strips, -3 -order 2) ----
//   h_junc  n_gap | 3D nodes | worst gamma | bad elems near strips
//    8x      1.35 |  121,433 |   1.17e-04  |  41.2%
//    6x      2    |  148,445 |   1.17e-04  |  34.0%   <-- current
//    4x      3    |  190,721 |   1.17e-04  |  36.5%
//   (ungraded A+B |   88,662 |   2.61e-09  |  40.9%)
// Drop to 8x / 1.35 if ElmerSolver runs short of memory.

L_half = (strip_w - strip_t)/2.0;          // half the flat edge length
R_grow = h_max/h_junc;                     // total size ratio to span
pr     = (L_half - h_junc)/(L_half - h_junc*R_grow);   // required Progression
n_half = Ceil(1.0 + Log(R_grow)/Log(pr)) + 1.0;        // nodes per half-flat

lc_bulk     = 0.5;
lc_strip    = h_max;    // point size at the strips now AGREES with the boundary
                        // discretisation there (field_cage's lc_wire=10um likewise
                        // matches its 7.85um arc segments -- that agreement is the
                        // thing that makes field_cage mesh cleanly)
lc_electrode = 0.01;

// ---------- 3. Gas vessel ----------
x_min = -22.5; x_max = 22.5;
y_min = -22.5; y_max = 22.5;

Point(1) = {x_min, y_min, 0, lc_bulk};
Point(2) = {x_max, y_min, 0, lc_bulk};
Point(3) = {x_max, y_max, 0, lc_bulk};
Point(4) = {x_min, y_max, 0, lc_bulk};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};

// ---------- 4. M-THGEM and Cathode ----------
Point(200) = {-gap/2.0, -16.01, 0, lc_electrode};
Point(201) = { gap/2.0, -16.01, 0, lc_electrode};
Point(202) = { gap/2.0, -15.99, 0, lc_electrode};
Point(203) = {-gap/2.0, -15.99, 0, lc_electrode};
Line(20) = {200,201}; Line(21) = {201,202}; Line(22) = {202,203}; Line(23) = {203,200};
Curve Loop(20) = {20,21,22,23};

Point(210) = {-gap/2.0, 15.99, 0, lc_electrode};
Point(211) = { gap/2.0, 15.99, 0, lc_electrode};
Point(212) = { gap/2.0, 16.01, 0, lc_electrode};
Point(213) = {-gap/2.0, 16.01, 0, lc_electrode};
Line(24) = {210,211}; Line(25) = {211,212}; Line(26) = {212,213}; Line(27) = {213,210};
Curve Loop(21) = {24,25,26,27};

// ---------- 5. Strips ----------
// Point stride is 20 (16 original + 4 new flat midpoints).
// Curve stride is 16 (12 original + 4 extra from splitting the flats).
p = 1000; l = 1000; cl = 5000;
strip_loops[] = {};

For i In {1:n_strips}
    y_pos = ((i - 1) * pitch) - y_offset;

    // ---- Left strip ----
    Point(p)   = {-gap/2.0 - strip_t/2.0, y_pos - strip_w/2.0 + strip_t/2.0, 0, lc_strip}; // Outer Bottom
    Point(p+1) = {-gap/2.0,               y_pos - strip_w/2.0 + strip_t/2.0, 0, lc_strip}; // Center Bottom
    Point(p+2) = {-gap/2.0 + strip_t/2.0, y_pos - strip_w/2.0 + strip_t/2.0, 0, lc_strip}; // Inner Bottom
    Point(p+3) = {-gap/2.0,               y_pos - strip_w/2.0,               0, lc_strip}; // Tip Bottom

    Point(p+4) = {-gap/2.0 + strip_t/2.0, y_pos + strip_w/2.0 - strip_t/2.0, 0, lc_strip}; // Inner Top
    Point(p+5) = {-gap/2.0,               y_pos + strip_w/2.0 - strip_t/2.0, 0, lc_strip}; // Center Top
    Point(p+6) = {-gap/2.0 - strip_t/2.0, y_pos + strip_w/2.0 - strip_t/2.0, 0, lc_strip}; // Outer Top
    Point(p+7) = {-gap/2.0,               y_pos + strip_w/2.0,               0, lc_strip}; // Tip Top

    // flat-edge midpoints (new)
    Point(p+16) = {-gap/2.0 + strip_t/2.0, y_pos, 0, lc_strip}; // inner flat midpoint
    Point(p+17) = {-gap/2.0 - strip_t/2.0, y_pos, 0, lc_strip}; // outer flat midpoint

    Circle(l)   = {p,   p+1, p+3};   // Outer Bottom Arc
    Circle(l+1) = {p+3, p+1, p+2};   // Inner Bottom Arc
    Line(l+2)   = {p+2,  p+16};      // Inner Flat, lower half  (fine -> coarse)
    Line(l+3)   = {p+16, p+4};       // Inner Flat, upper half  (coarse -> fine)
    Circle(l+4) = {p+4, p+5, p+7};   // Inner Top Arc
    Circle(l+5) = {p+7, p+5, p+6};   // Outer Top Arc
    Line(l+6)   = {p+6,  p+17};      // Outer Flat, upper half  (fine -> coarse)
    Line(l+7)   = {p+17, p};         // Outer Flat, lower half  (coarse -> fine)

    Curve Loop(cl) = {l, l+1, l+2, l+3, l+4, l+5, l+6, l+7};
    strip_loops[] += {cl};

    // ---- Right strip (mirror image, identical mesh controls) ----
    Point(p+8)  = {gap/2.0 - strip_t/2.0, y_pos - strip_w/2.0 + strip_t/2.0, 0, lc_strip}; // Inner Bottom
    Point(p+9)  = {gap/2.0,               y_pos - strip_w/2.0 + strip_t/2.0, 0, lc_strip}; // Center Bottom
    Point(p+10) = {gap/2.0 + strip_t/2.0, y_pos - strip_w/2.0 + strip_t/2.0, 0, lc_strip}; // Outer Bottom
    Point(p+11) = {gap/2.0,               y_pos - strip_w/2.0,               0, lc_strip}; // Tip Bottom

    Point(p+12) = {gap/2.0 + strip_t/2.0, y_pos + strip_w/2.0 - strip_t/2.0, 0, lc_strip}; // Outer Top
    Point(p+13) = {gap/2.0,               y_pos + strip_w/2.0 - strip_t/2.0, 0, lc_strip}; // Center Top
    Point(p+14) = {gap/2.0 - strip_t/2.0, y_pos + strip_w/2.0 - strip_t/2.0, 0, lc_strip}; // Inner Top
    Point(p+15) = {gap/2.0,               y_pos + strip_w/2.0,               0, lc_strip}; // Tip Top

    Point(p+18) = {gap/2.0 + strip_t/2.0, y_pos, 0, lc_strip}; // outer flat midpoint
    Point(p+19) = {gap/2.0 - strip_t/2.0, y_pos, 0, lc_strip}; // inner flat midpoint

    Circle(l+8)  = {p+8,  p+9,  p+11}; // Inner Bottom Arc
    Circle(l+9)  = {p+11, p+9,  p+10}; // Outer Bottom Arc
    Line(l+10)   = {p+10, p+18};       // Outer Flat, lower half  (fine -> coarse)
    Line(l+11)   = {p+18, p+12};       // Outer Flat, upper half  (coarse -> fine)
    Circle(l+12) = {p+12, p+13, p+15}; // Outer Top Arc
    Circle(l+13) = {p+15, p+13, p+14}; // Inner Top Arc
    Line(l+14)   = {p+14, p+19};       // Inner Flat, upper half  (fine -> coarse)
    Line(l+15)   = {p+19, p+8};        // Inner Flat, lower half  (coarse -> fine)

    Curve Loop(cl+1) = {l+8, l+9, l+10, l+11, l+12, l+13, l+14, l+15};
    strip_loops[] += {cl+1};

    // ---- Mesh controls: applied IDENTICALLY to both strips ----
    // caps
    Transfinite Curve {l, l+1, l+4, l+5, l+8, l+9, l+12, l+13} = n_arc;
    // flat halves running fine -> coarse
    Transfinite Curve {l+2, l+6, l+10, l+14} = n_half Using Progression pr;
    // flat halves running coarse -> fine (mirrored grading)
    Transfinite Curve {l+3, l+7, l+11, l+15} = n_half Using Progression 1.0/pr;

    p += 20; l += 16; cl += 2;
EndFor

// ---------- 6. Surface and extrusion ----------
Plane Surface(1) = {1, strip_loops[], 20, 21};

// Layers{1} is REQUIRED here, and it does NOT cost you tetrahedra.
// Without it, Gmsh must fill the 100 um slab by unstructured 3D boundary
// recovery; with the graded strip boundary present that fails outright
// ("Could not recover boundary mesh: error 2" / "No elements in volume 1").
// With Layers{1} and no Recombine, the extrusion is structured in z and each
// prism is split into tets -- the mesh is pure tet10, exactly what
// ComponentElmer expects.  Verified: 3D element types = ['tet10'].
ext[] = Extrude {0, 0, z_thick} { Surface{1}; Layers{1}; };

// ---------- 7. Physical groups ----------
Physical Volume(1) = {ext[1]};

Physical Surface(1) = {
    Surface In BoundingBox{-22.6, -22.6, -0.1, -22.4,  22.6, z_thick + 0.1}, // Left
    Surface In BoundingBox{ 22.4, -22.6, -0.1,  22.6,  22.6, z_thick + 0.1}, // Right
    Surface In BoundingBox{-22.6, -22.6, -0.1,  22.6, -22.4, z_thick + 0.1}, // Bottom
    Surface In BoundingBox{-22.6,  22.4, -0.1,  22.6,  22.6, z_thick + 0.1}  // Top
};

Physical Surface(2) = Surface In BoundingBox{-gap/2.0-0.1, -16.1, -0.1, gap/2.0+0.1, -15.9, z_thick+0.1}; // M-THGEM
Physical Surface(3) = Surface In BoundingBox{-gap/2.0-0.1,  15.9, -0.1, gap/2.0+0.1,  16.1, z_thick+0.1}; // Cathode

For i In {1:n_strips}
    y_pos = ((i - 1) * pitch) - y_offset;
    Physical Surface(3 + i) = Surface In BoundingBox{-gap/2.0 - strip_t/2.0 - 0.01, y_pos - strip_w/2.0 - 0.01, -0.1,
                                                      gap/2.0 + strip_t/2.0 + 0.01, y_pos + strip_w/2.0 + 0.01, z_thick + 0.1};
EndFor

// ---------- 8. Meshing flags ----------

// ---------- Boundary-layer refinement around the strip planes ----------
// The strip staircase creates a field ripple that decays with length
// pitch/2pi = 404 um.  Resolving it needs h << 404 um within ~1 mm of the
// strip plane.  Built from Fabs(Fabs(x)-gap/2) so it is EXACTLY symmetric in
// x -- no hard edges anywhere (unlike a Box field).
h_near = 0.020;    // element size inside the boundary layer  (200 um)
r_band = 0.05;     // half-width of the refined band (cm)
grow   = 2.0;      // how fast it relaxes outside the band

// ---- MEASURED (full 118 strips, -3 -order 2) ----
//  h_near  r_band | 3D nodes | worst gamma | elems<0.01 | h at x=1.49
//   400um   0.5mm |    ~332k |      -      |     -      |    ~70 um
//   300um   0.5mm |    ~379k |      -      |     -      |    ~55 um
//   200um   0.5mm |     584k |   1.10e-02  |   0.0000%  |     36 um   <-- current
//   150um   1.0mm |   1,146k |      -      |     -      |     -       (OOM risk)
//  (previous version |   148k |   1.17e-04  |  39.3%     |    500 um)
// If ElmerSolver runs out of memory, step h_near up to 0.030 then 0.040.
Field[1] = MathEval;
Field[1].F = Sprintf("%g + %g*( (Sqrt((Fabs(Fabs(x)-%g))^2 + ((Fabs(y)-15.1+Fabs(Fabs(y)-15.1))/2)^2) - %g) + Fabs(Sqrt((Fabs(Fabs(x)-%g))^2 + ((Fabs(y)-15.1+Fabs(Fabs(y)-15.1))/2)^2) - %g) )/2",
                     h_near, grow, gap/2.0, r_band, gap/2.0, r_band);
Background Field = 1;
Mesh.CharacteristicLengthMax = lc_bulk;

Mesh.Algorithm = 6;
// Keep this at 0, as field_cage does.  Setting it to 1 was tested and is a
// trap here: it propagates the 3.9 um cap size across the whole 45 cm vessel
// and the 2D mesh never finishes.  It is not needed, because lc_strip is now
// tied to h_max -- the point size and the coarse end of the boundary grading
// agree, so the gas mesh leaves the strip at the right size on its own.
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.SecondOrderLinear = 1;

