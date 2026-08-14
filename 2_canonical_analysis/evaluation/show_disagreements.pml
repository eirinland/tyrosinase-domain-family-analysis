reinitialize

# Load all structures
load 00_ref_B2ZB02_PmTYR.cif, ref
load 01_A0ABQ5PAJ6_taxID_2989713_model.cif, s01
load 02_A0A1Y1Z2V1_taxID_1314790_model.cif, s02
load 03_A0A1Y1XAL5_taxID_1314790_model.cif, s03
load 04_A0ABS7Q4E0_taxID_2873382_model.cif, s04
load 05_A0ABW0UH65_taxID_349910_model.cif, s05
load 06_A0A4P9XGI5_taxID_78915_model.cif, s06
load 07_A0A4P9Y249_taxID_1907219_model.cif, s07
load 08_A0A7S2FUY7_taxID_3111310_model.cif, s08
load 09_A0A7S2R8C7_taxID_1034831_model.cif, s09
load 10_F0Y8A8_taxID_44056_model.cif, s10

# Align all to reference
align s01, ref
align s02, ref
align s03, ref
align s04, ref
align s05, ref
align s06, ref
align s07, ref
align s08, ref
align s09, ref
align s10, ref

# Show everything as cartoon, hide all sticks/labels
hide everything
show cartoon, all
set cartoon_transparency, 0.7

# Color reference
color gray70, ref

# Colors per structure
color marine, s01
color teal, s02
color cyan, s03
color orange, s04
color tv_orange, s05
color salmon, s06
color raspberry, s07
color red, s08
color firebrick, s09
color chocolate, s10

# ---- Disagreement residues as sticks ----
# Kabsch = magenta sticks, Foldseek = green sticks

# s01 (A0ABQ5PAJ6, qTM 0.94): His204 and Ala221
select s01_kab, s01 and resi 183+199
select s01_fs, none
show sticks, s01_kab
color magenta, s01_kab

# s02 (A0A1Y1Z2V1, qTM 0.89): His204 and Ala221
select s02_kab, s02 and resi 164+180
select s02_fs, none
show sticks, s02_kab
color magenta, s02_kab

# s03 (A0A1Y1XAL5, qTM 0.88): His204 and Ala221
select s03_kab, s03 and resi 163+179
select s03_fs, none
show sticks, s03_kab
color magenta, s03_kab

# s04 (A0ABS7Q4E0, qTM 0.87): His42 and Gly46
# Kabsch: H38, F42  Foldseek: F42, G46
select s04_kab, s04 and resi 38
select s04_fs, s04 and resi 46
show sticks, s04_kab
show sticks, s04_fs
color magenta, s04_kab
color splitpea, s04_fs
# resi 42 is claimed by both (kabsch Gly46->F42, foldseek His42->F42)
select s04_shared, s04 and resi 42
show sticks, s04_shared
color yellow, s04_shared

# s05 (A0ABW0UH65, qTM 0.86): His42 and Gly46
# Kabsch: H39, Y43  Foldseek: Y43, G47
select s05_kab, s05 and resi 39
select s05_fs, s05 and resi 47
show sticks, s05_kab
show sticks, s05_fs
color magenta, s05_kab
color splitpea, s05_fs
select s05_shared, s05 and resi 43
show sticks, s05_shared
color yellow, s05_shared

# s06 (A0A4P9XGI5, qTM 0.85): Gly46, His60, Glu195, Ala221
# Kabsch: A42, H46, E165, S185  Foldseek: N43, -, L162, -
select s06_kab, s06 and resi 42+46+165+185
select s06_fs, s06 and resi 43+162
show sticks, s06_kab
show sticks, s06_fs
color magenta, s06_kab
color splitpea, s06_fs

# s07 (A0A4P9Y249, qTM 0.86): Glu195, His204
# Kabsch: T161, H164  Foldseek: E158, -
select s07_kab, s07 and resi 161+164
select s07_fs, s07 and resi 158
show sticks, s07_kab
show sticks, s07_fs
color magenta, s07_kab
color splitpea, s07_fs

# s08 (A0A7S2FUY7, qTM 0.41): 8 disagreements
# Kabsch: S31,H43,E186,H190,S191,H194,Q195,L342
# Foldseek: K27,V42,T289,E304,Y305,R327,I328,P341
select s08_kab, s08 and resi 31+43+186+190+191+194+195+342
select s08_fs, s08 and resi 27+42+289+304+305+327+328+341
show sticks, s08_kab
show sticks, s08_fs
color magenta, s08_kab
color splitpea, s08_fs

# s09 (A0A7S2R8C7, qTM 0.40): 8 disagreements
# Kabsch: A15,Y168,H173,G174,H177,M179,T295,A298
# Foldseek: D17,E252,D277,A278,I281,R282,S296,-
select s09_kab, s09 and resi 15+168+173+174+177+179+295+298
select s09_fs, s09 and resi 17+252+277+278+281+282+296
show sticks, s09_kab
show sticks, s09_fs
color magenta, s09_kab
color splitpea, s09_fs

# s10 (F0Y8A8, qTM 0.42): 8 disagreements
# Kabsch: A58,L212,H217,G218,H221,M223,T352,A355
# Foldseek: S60,E339,T347,-,-,-,S353,-
select s10_kab, s10 and resi 58+212+217+218+221+223+352+355
select s10_fs, s10 and resi 60+339+347+353
show sticks, s10_kab
show sticks, s10_fs
color magenta, s10_kab
color splitpea, s10_fs

# Show reference His + key positions as sticks for context
select ref_his, ref and resi 42+60+69+204+208+230+231 and name CA+CB+CG+ND1+CE1+NE2+CD2
show sticks, ref_his
color white, ref_his

# Show Cu in reference
show spheres, ref and name CU
color copper, ref and name CU
set sphere_scale, 0.5

# Show Cu in all query structures
show spheres, (s01 or s02 or s03 or s04 or s05 or s06 or s07 or s08 or s09 or s10) and name CU
color copper, (s01 or s02 or s03 or s04 or s05 or s06 or s07 or s08 or s09 or s10) and name CU

# Clean up
deselect
set stick_radius, 0.15
set label_size, 12
set label_color, black
bg_color white
set ray_opaque_background, 1
set ray_shadows, 0
zoom ref

# Disable all except ref initially -- enable one at a time to inspect
disable s01
disable s02
disable s03
disable s04
disable s05
disable s06
disable s07
disable s08
disable s09
disable s10

# Group by disagreement type for easy toggling
group high_qtm, s01 s02 s03
group shift_4res, s04 s05
group mixed_cub, s06 s07
group low_qtm, s08 s09 s10
