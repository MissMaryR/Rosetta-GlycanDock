#!/bin/bash
# Step 1: Generate a Rosetta-native saccharide in RosettaCarbohydrate format.
# Run this ONCE on the cluster — takes ~10 seconds.
# Output: out_saccharide.pdb  (chain X, correct ->4)-beta-D-Glcp residue types)
#
# Edit ROSETTA_BIN, ROSETTA_DB, and the saccharide_sequence for your system.

ROSETTA_BIN=/path/to/rosetta/main/source/bin
ROSETTA_DB=/path/to/rosetta/main/database

$ROSETTA_BIN/pose_from_saccharide_sequence.static.linuxgccrelease \
  -database $ROSETTA_DB \
  -include_sugars \
  -alternate_3_letter_codes pdb_sugar \
  -carbohydrates:saccharide_sequence "b-D-Glcp-(1->4)-b-D-Glcp-(1->4)-b-D-Glcp" \
  -nstruct 1

# Output is named out_saccharide.pdb by Rosetta
echo "Output: out_saccharide.pdb"
echo "Copy this file to your local machine for the PyMOL alignment step (2_align)."
