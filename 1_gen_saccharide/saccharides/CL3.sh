#!/bin/bash
# Generate a Rosetta-native cellotriose in RosettaCarbohydrate format.
# Run this ONCE on the cluster — takes ~10 seconds.
# Output: out_saccharide.pdb  (3x BGC on chain X, correct ->4) residue types)

/quobyte/jbsiegelgrp/software/Rosetta_314/rosetta/main/source/bin/pose_from_saccharide_sequence.static.linuxgccrelease \
  -database /quobyte/jbsiegelgrp/software/Rosetta_314/rosetta/main/database \
  -include_sugars \
  -alternate_3_letter_codes pdb_sugar \
  -carbohydrates:saccharide_sequence "b-D-Glcp-(1->4)-b-D-Glcp-(1->4)-b-D-Glcp" \
  -nstruct 1

mv out_saccharide.pdb CL3.pdb
echo "Output: CL3.pdb"
