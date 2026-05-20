#!/bin/bash
# Generate a Rosetta-native H2B disaccharide + xylose branch in RosettaCarbohydrate format.
# Structure: 2x beta-D-glucopyranose beta-(1->4), with alpha-D-xylopyranose alpha-(1->6)
#            branch on the REDUCING-END glucose.
#
# Sequence reads left (non-reducing end) to right (reducing end).
# [a-D-Xylp-(1->6)] in brackets before b-D-Glcp means:
#   Xylp C1 connects to C6 of the following Glcp (the reducing-end glucose).
#
# Run this ONCE on the cluster — takes ~10 seconds.
# Output: H2B.pdb

/quobyte/jbsiegelgrp/software/Rosetta_314/rosetta/main/source/bin/pose_from_saccharide_sequence.static.linuxgccrelease \
  -database /quobyte/jbsiegelgrp/software/Rosetta_314/rosetta/main/database \
  -include_sugars \
  -alternate_3_letter_codes pdb_sugar \
  -carbohydrates:saccharide_sequence "b-D-Glcp-(1->4)-[a-D-Xylp-(1->6)]-b-D-Glcp" \
  -nstruct 1

mv out_saccharide.pdb H2B.pdb
echo "Output: H2B.pdb"
