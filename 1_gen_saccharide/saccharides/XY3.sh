#!/bin/bash
# Generate a Rosetta-native xylotriose in RosettaCarbohydrate format.
# 3x beta-D-xylopyranose linked beta-(1->4)
# Run this ONCE on the cluster — takes ~10 seconds.
# Output: XY3.pdb (3x XYL on chain X, correct ->4) residue types)

/quobyte/jbsiegelgrp/software/Rosetta_314/rosetta/main/source/bin/pose_from_saccharide_sequence.static.linuxgccrelease \
  -database /quobyte/jbsiegelgrp/software/Rosetta_314/rosetta/main/database \
  -include_sugars \
  -alternate_3_letter_codes pdb_sugar \
  -carbohydrates:saccharide_sequence "b-D-Xylp-(1->4)-b-D-Xylp-(1->4)-b-D-Xylp" \
  -nstruct 1

mv out_saccharide.pdb XY3.pdb
echo "Output: XY3.pdb"
