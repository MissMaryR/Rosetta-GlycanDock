#!/bin/bash
# Generate a Rosetta-native H3B trisaccharide + xylose branch in RosettaCarbohydrate format.
# Structure: 3x beta-D-glucopyranose beta-(1->4), with alpha-D-xylopyranose alpha-(1->6)
#            branch on the NON-REDUCING-END glucose.
#
# Sequence reads left (non-reducing end) to right (reducing end).
# [a-D-Xylp-(1->6)] in brackets before the first b-D-Glcp means:
#   Xylp C1 connects to C6 of the non-reducing-end glucose.
#
# Run this ONCE on the cluster — takes ~10 seconds.
# Output: H3B.pdb

/quobyte/jbsiegelgrp/software/Rosetta_314/rosetta/main/source/bin/pose_from_saccharide_sequence.static.linuxgccrelease \
  -database /quobyte/jbsiegelgrp/software/Rosetta_314/rosetta/main/database \
  -include_sugars \
  -alternate_3_letter_codes pdb_sugar \
  -carbohydrates:saccharide_sequence "[a-D-Xylp-(1->6)]-b-D-Glcp-(1->4)-b-D-Glcp-(1->4)-b-D-Glcp" \
  -nstruct 1

mv out_saccharide.pdb H3B.pdb
echo "Output: H3B.pdb"
