# Borderline rescued-residue gallery (multi-ref helix-anchored superposition)
# gray90 = rest | orange spheres = Cu | dashed lines = saved-residue CA -> Cu
# cyan sticks = the other 5 mapped Cu-site anchor positions (mostly His; some
#   sites are already substituted, e.g. s03 has 0 His left at the mapped anchors)
# GREEN sticks = residue 'saved' when the 3.0 A Cu-anchor cutoff is relaxed
# (label N_<cutoff>Årescue = recovered at that cutoff; notsaved5Å = still gapped >5 A)
#
# Inspect one: set grid_mode,0; disable all; enable s05*; zoom s05* and elem Cu, 12
reinitialize
bg_color white
set cartoon_transparency, 0.1
set sphere_scale, 0.4, elem Cu
set label_size, 16
set label_color, black
set dash_color, gray40

# --- s01_3p5_A0A177T6E5  ref=ref_PmTYR qtm=0.529  saved: CuA3 TYR140 3.01Å
load "1_3.5Årescue_A0A177T6E5.cif", s01_3p5_A0A177T6E5
set grid_slot, 1, s01_3p5_A0A177T6E5
set grid_slot, 1, d01
hide everything, s01_3p5_A0A177T6E5
show cartoon, s01_3p5_A0A177T6E5
color gray90, s01_3p5_A0A177T6E5
show spheres, s01_3p5_A0A177T6E5 and elem Cu
color orange, s01_3p5_A0A177T6E5 and elem Cu
select cx, s01_3p5_A0A177T6E5 and resi 146+135+295+317+291
show sticks, cx and not name C+N+O
color cyan, cx
select gx, s01_3p5_A0A177T6E5 and resi 140
show sticks, gx and not name C+N+O
color green, gx
label gx and name CA, "%s%s" % (resn, resi)
distance d01, gx and name CA, s01_3p5_A0A177T6E5 and elem Cu

# --- s02_3p5_A0A161WK35  ref=ref_2Y9W_Abisporus qtm=0.499  saved: CuA1 GLN99 3.03Å
load "2_3.5Årescue_A0A161WK35.cif", s02_3p5_A0A161WK35
set grid_slot, 2, s02_3p5_A0A161WK35
set grid_slot, 2, d02
hide everything, s02_3p5_A0A161WK35
show cartoon, s02_3p5_A0A161WK35
color gray90, s02_3p5_A0A161WK35
show spheres, s02_3p5_A0A161WK35 and elem Cu
color orange, s02_3p5_A0A161WK35 and elem Cu
select cx, s02_3p5_A0A161WK35 and resi 86+93+286+257+261
show sticks, cx and not name C+N+O
color cyan, cx
select gx, s02_3p5_A0A161WK35 and resi 99
show sticks, gx and not name C+N+O
color green, gx
label gx and name CA, "%s%s" % (resn, resi)
distance d02, gx and name CA, s02_3p5_A0A161WK35 and elem Cu

# --- s03_4p0_A0A074W6F3  ref=ref_PmTYR qtm=0.491  saved: CuB3 ILE246 3.54Å
load "3_4.0Årescue_A0A074W6F3.cif", s03_4p0_A0A074W6F3
set grid_slot, 3, s03_4p0_A0A074W6F3
set grid_slot, 3, d03
hide everything, s03_4p0_A0A074W6F3
show cartoon, s03_4p0_A0A074W6F3
color gray90, s03_4p0_A0A074W6F3
show spheres, s03_4p0_A0A074W6F3 and elem Cu
color orange, s03_4p0_A0A074W6F3 and elem Cu
select cx, s03_4p0_A0A074W6F3 and resi 162+145+153+252+268
show sticks, cx and not name C+N+O
color cyan, cx
select gx, s03_4p0_A0A074W6F3 and resi 246
show sticks, gx and not name C+N+O
color green, gx
label gx and name CA, "%s%s" % (resn, resi)
distance d03, gx and name CA, s03_4p0_A0A074W6F3 and elem Cu

# --- s04_4p0_A0A9R1RCW5  ref=ref_5CE9_Jregia qtm=0.558  saved: CuA1 PRO182 3.51Å
load "4_4.0Årescue_A0A9R1RCW5.cif", s04_4p0_A0A9R1RCW5
set grid_slot, 4, s04_4p0_A0A9R1RCW5
set grid_slot, 4, d04
hide everything, s04_4p0_A0A9R1RCW5
show cartoon, s04_4p0_A0A9R1RCW5
color gray90, s04_4p0_A0A9R1RCW5
show spheres, s04_4p0_A0A9R1RCW5 and elem Cu
color orange, s04_4p0_A0A9R1RCW5 and elem Cu
select cx, s04_4p0_A0A9R1RCW5 and resi 168+186+342+308+312
show sticks, cx and not name C+N+O
color cyan, cx
select gx, s04_4p0_A0A9R1RCW5 and resi 182
show sticks, gx and not name C+N+O
color green, gx
label gx and name CA, "%s%s" % (resn, resi)
distance d04, gx and name CA, s04_4p0_A0A9R1RCW5 and elem Cu

# --- s05_4p5_A0A0D2UIJ1  ref=ref_1BT3_Ibatatas qtm=0.505  saved: CuA3 VAL281 4.01Å
load "5_4.5Årescue_A0A0D2UIJ1.cif", s05_4p5_A0A0D2UIJ1
set grid_slot, 5, s05_4p5_A0A0D2UIJ1
set grid_slot, 5, d05
hide everything, s05_4p5_A0A0D2UIJ1
show cartoon, s05_4p5_A0A0D2UIJ1
color gray90, s05_4p5_A0A0D2UIJ1
show spheres, s05_4p5_A0A0D2UIJ1 and elem Cu
color orange, s05_4p5_A0A0D2UIJ1 and elem Cu
select cx, s05_4p5_A0A0D2UIJ1 and resi 171+182+286+290+320
show sticks, cx and not name C+N+O
color cyan, cx
select gx, s05_4p5_A0A0D2UIJ1 and resi 281
show sticks, gx and not name C+N+O
color green, gx
label gx and name CA, "%s%s" % (resn, resi)
distance d05, gx and name CA, s05_4p5_A0A0D2UIJ1 and elem Cu

# --- s06_4p5_A0A3P8DS39  ref=ref_PmTYR qtm=0.524  saved: CuB3 TRP286 4.03Å
load "6_4.5Årescue_A0A3P8DS39.cif", s06_4p5_A0A3P8DS39
set grid_slot, 6, s06_4p5_A0A3P8DS39
set grid_slot, 6, d06
hide everything, s06_4p5_A0A3P8DS39
show cartoon, s06_4p5_A0A3P8DS39
color gray90, s06_4p5_A0A3P8DS39
show spheres, s06_4p5_A0A3P8DS39 and elem Cu
color orange, s06_4p5_A0A3P8DS39 and elem Cu
select cx, s06_4p5_A0A3P8DS39 and resi 171+152+162+289+312
show sticks, cx and not name C+N+O
color cyan, cx
select gx, s06_4p5_A0A3P8DS39 and resi 286
show sticks, gx and not name C+N+O
color green, gx
label gx and name CA, "%s%s" % (resn, resi)
distance d06, gx and name CA, s06_4p5_A0A3P8DS39 and elem Cu

# --- s07_5p0_A0A8H6K3D9  ref=ref_2Y9W_Abisporus qtm=0.507  saved: CuA2 PRO266 4.53Å
load "7_5.0Årescue_A0A8H6K3D9.cif", s07_5p0_A0A8H6K3D9
set grid_slot, 7, s07_5p0_A0A8H6K3D9
set grid_slot, 7, d07
hide everything, s07_5p0_A0A8H6K3D9
show cartoon, s07_5p0_A0A8H6K3D9
color gray90, s07_5p0_A0A8H6K3D9
show spheres, s07_5p0_A0A8H6K3D9 and elem Cu
color orange, s07_5p0_A0A8H6K3D9 and elem Cu
select cx, s07_5p0_A0A8H6K3D9 and resi 283+274+473+443+447
show sticks, cx and not name C+N+O
color cyan, cx
select gx, s07_5p0_A0A8H6K3D9 and resi 266
show sticks, gx and not name C+N+O
color green, gx
label gx and name CA, "%s%s" % (resn, resi)
distance d07, gx and name CA, s07_5p0_A0A8H6K3D9 and elem Cu

# --- s08_5p0_A0A3P3Y6K2  ref=ref_1JS8_squid qtm=0.392  saved: CuA1 GLY65 4.6Å
load "8_5.0Årescue_A0A3P3Y6K2.cif", s08_5p0_A0A3P3Y6K2
set grid_slot, 8, s08_5p0_A0A3P3Y6K2
set grid_slot, 8, d08
hide everything, s08_5p0_A0A3P3Y6K2
show cartoon, s08_5p0_A0A3P3Y6K2
color gray90, s08_5p0_A0A3P3Y6K2
show spheres, s08_5p0_A0A3P3Y6K2 and elem Cu
color orange, s08_5p0_A0A3P3Y6K2 and elem Cu
select cx, s08_5p0_A0A3P3Y6K2 and resi 58+73+196+173+169
show sticks, cx and not name C+N+O
color cyan, cx
select gx, s08_5p0_A0A3P3Y6K2 and resi 65
show sticks, gx and not name C+N+O
color green, gx
label gx and name CA, "%s%s" % (resn, resi)
distance d08, gx and name CA, s08_5p0_A0A3P3Y6K2 and elem Cu

# --- s09_unsaved_A0A5C6SLG5  ref=ref_PmTYR qtm=0.578  saved: CuA3 ILE258 7.85Å
load "9_notsaved5Å_A0A5C6SLG5.cif", s09_unsaved_A0A5C6SLG5
set grid_slot, 9, s09_unsaved_A0A5C6SLG5
set grid_slot, 9, d09
hide everything, s09_unsaved_A0A5C6SLG5
show cartoon, s09_unsaved_A0A5C6SLG5
color gray90, s09_unsaved_A0A5C6SLG5
show spheres, s09_unsaved_A0A5C6SLG5 and elem Cu
color orange, s09_unsaved_A0A5C6SLG5 and elem Cu
select cx, s09_unsaved_A0A5C6SLG5 and resi 79+75+208+230+204
show sticks, cx and not name C+N+O
color cyan, cx
select gx, s09_unsaved_A0A5C6SLG5 and resi 258
show sticks, gx and not name C+N+O
color green, gx
label gx and name CA, "%s%s" % (resn, resi)
distance d09, gx and name CA, s09_unsaved_A0A5C6SLG5 and elem Cu

# --- s10_unsaved_A0A8H8SYG4  ref=ref_PmTYR qtm=0.639  saved: CuA2 VAL66 3.33Å; CuA3 LYS248 7.24Å
load "10_notsaved5Å_A0A8H8SYG4.cif", s10_unsaved_A0A8H8SYG4
set grid_slot, 10, s10_unsaved_A0A8H8SYG4
set grid_slot, 10, d10
hide everything, s10_unsaved_A0A8H8SYG4
show cartoon, s10_unsaved_A0A8H8SYG4
color gray90, s10_unsaved_A0A8H8SYG4
show spheres, s10_unsaved_A0A8H8SYG4 and elem Cu
color orange, s10_unsaved_A0A8H8SYG4 and elem Cu
select cx, s10_unsaved_A0A8H8SYG4 and resi 249+185+207+181
show sticks, cx and not name C+N+O
color cyan, cx
select gx, s10_unsaved_A0A8H8SYG4 and resi 66+248
show sticks, gx and not name C+N+O
color green, gx
label gx and name CA, "%s%s" % (resn, resi)
distance d10, gx and name CA, s10_unsaved_A0A8H8SYG4 and elem Cu

delete cx gx
deselect
set grid_mode, 1
orient
zoom all, 3
set ray_shadows, 0
