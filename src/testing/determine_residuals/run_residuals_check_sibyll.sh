#!/bin/bash

  datdir=/home/gkukec/Gasper/github/auger-analysis/results/analysis_new_20181128/v3r3p4-production/Auger_data_application/mva_analysis/sibyll

  make

  # Residuals from PAO data
  # Relative versus absolute observables
  echo "Checking residuals for comparison relative vs. absolute observables"
  ./determine_residuals $datdir/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS6.root $datdir/xmax_s1000_rise_zenithFD/plots/lna_composition_results_TreeS6.root > $datdir/xmax_deltas38_deltarise/plots/resid_rel-abs_printout.txt

  # Residuals from AugerMix
  # Relative versus absolute observables
  echo "Checking residuals for comparison relative vs. absolute observables"
  ./determine_residuals $datdir/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS5.root $datdir/xmax_s1000_rise_zenithFD/plots/lna_composition_results_TreeS5.root > $datdir/xmax_deltas38_deltarise/plots/resid_rel-abs_AugerMix_printout.txt

  # Relative versus relative+theta observables
#  echo "Checking residuals for comparison relative vs. absolute observables"
#  ./determine_residuals $datdir/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS6.root $datdir/xmax_deltas38_deltarise_zenithFD/plots/lna_composition_results_TreeS6.root > $datdir/xmax_deltas38_deltarise/plots/resid_rel-rel-zenith_printout.txt

  # Relative versus relative+theta observables
#  echo "Checking residuals for comparison relative vs. absolute observables"
#  ./determine_residuals $datdir/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS5.root $datdir/xmax_deltas38_deltarise_zenithFD/plots/lna_composition_results_TreeS5.root > $datdir/xmax_deltas38_deltarise/plots/resid_rel-rel-zenith_AugerMix_printout.txt
