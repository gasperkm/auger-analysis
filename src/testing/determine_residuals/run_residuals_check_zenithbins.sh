#!/bin/bash

  datdir=/home/gkukec/Gasper/github/auger-analysis/results/analysis_new_20190121/zenith_angle

  make

  # Residuals from PAO data (zenith binning)
  # sec theta = 1-1.5 versus 1.5-2
  echo "Checking residuals for zenith binning (data)"
  ./determine_residuals $datdir/epos_relative_1-1.5/plots/lna_composition_results_TreeS6.root $datdir/epos_relative_1.5-2/plots/lna_composition_results_TreeS6.root > $datdir/epos_relative_1-1.5/plots/resid_rel-abs_printout.txt

  # Residuals from AugerMix
  # Relative versus absolute observables
  echo "Checking residuals for zenith binning (AugerMix)"
  ./determine_residuals $datdir/epos_relative_1-1.5/plots/lna_composition_results_TreeS5.root $datdir/epos_relative_1.5-2/plots/lna_composition_results_TreeS5.root > $datdir/epos_relative_1-1.5/plots/resid_rel-abs_printout.txt

  # Relative versus relative+theta observables
#  echo "Checking residuals for comparison relative vs. absolute observables"
#  ./determine_residuals $datdir/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS6.root $datdir/xmax_deltas38_deltarise_zenithFD/plots/lna_composition_results_TreeS6.root > $datdir/xmax_deltas38_deltarise/plots/resid_rel-rel-zenith_printout.txt

  # Relative versus relative+theta observables
#  echo "Checking residuals for comparison relative vs. absolute observables"
#  ./determine_residuals $datdir/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS5.root $datdir/xmax_deltas38_deltarise_zenithFD/plots/lna_composition_results_TreeS5.root > $datdir/xmax_deltas38_deltarise/plots/resid_rel-rel-zenith_AugerMix_printout.txt
