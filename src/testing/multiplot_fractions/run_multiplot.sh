#!/bin/bash

  datdir=/home/gkukec/Gasper/github/auger-analysis/results/analysis_new_20181128/v3r3p4-production/Auger_data_application/mva_analysis

  make

  # Multiplot from PAO data (all three models, relative)
  echo "Creating multiplot for all three models (relative, PAO)"
  ./multiplot_fractions $datdir/epos/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS6.root $datdir/qgsjet/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS6.root $datdir/sibyll/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS6.root > $datdir/epos/xmax_deltas38_deltarise/plots/multiplot_relative_printout.txt

  # Multiplot from PAO data (all three models, absolute)
  echo "Creating multiplot for all three models (absolute, PAO)"
  ./multiplot_fractions $datdir/epos/xmax_s1000_rise_zenithFD/plots/lna_composition_results_TreeS6.root $datdir/qgsjet/xmax_s1000_rise_zenithFD/plots/lna_composition_results_TreeS6.root $datdir/sibyll/xmax_s1000_rise_zenithFD/plots/lna_composition_results_TreeS6.root > $datdir/epos/xmax_s1000_rise_zenithFD/plots/multiplot_absolute_printout.txt

  # Multiplot from AugerMix (all three models, relative)
  echo "Creating multiplot for all three models (relative, AugerMix)"
  ./multiplot_fractions $datdir/epos/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS5.root $datdir/qgsjet/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS5.root $datdir/sibyll/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS5.root > $datdir/epos/xmax_deltas38_deltarise/plots/multiplot_relative_AugerMix_printout.txt

  # Multiplot from AugerMix (all three models, absolute)
  echo "Creating multiplot for all three models (absolute, AugerMix)"
  ./multiplot_fractions $datdir/epos/xmax_s1000_rise_zenithFD/plots/lna_composition_results_TreeS5.root $datdir/qgsjet/xmax_s1000_rise_zenithFD/plots/lna_composition_results_TreeS5.root $datdir/sibyll/xmax_s1000_rise_zenithFD/plots/lna_composition_results_TreeS5.root > $datdir/epos/xmax_s1000_rise_zenithFD/plots/multiplot_absolute_AugerMix_printout.txt
