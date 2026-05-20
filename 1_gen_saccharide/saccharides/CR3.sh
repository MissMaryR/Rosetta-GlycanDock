#!/bin/bash
# Generate a Rosetta-native beta-1,3-glucotrioside (curdlan/laminarin-type) in RosettaCarbohydrate format.
# 3x beta-D-glucopyranose linked beta-(1->3)
# Run this ONCE on the cluster — takes ~10 seconds.
# Output: CR3.pdb (3x GLC on chain X, correct ->3) residue types)

/quobyte/jbsiegelgrp/software/Rosetta_314/rosetta/main/source/bin/pose_from_saccharide_sequence.static.linuxgccrelease \
  -database /quobyte/jbsiegelgrp/software/Rosetta_314/rosetta/main/database \
  -include_sugars \
  -alternate_3_letter_codes pdb_sugar \
  -carbohydrates:saccharide_sequence "b-D-Glcp-(1->3)-b-D-Glcp-(1->3)-b-D-Glcp" \
  -nstruct 1

mv out_saccharide.pdb CR3.pdb
echo "Output: CR3.pdb"
