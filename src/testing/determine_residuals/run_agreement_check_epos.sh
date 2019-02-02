#!/bin/bash

  datdir=/home/gkukec/Gasper/github/auger-analysis/results/analysis_new_20181128/v3r3p4-production/Auger_data_application

  make

  # Agreement between PAO data and published data
  echo "Checking agreement between PAO data and published data (FD-only)"
  ./determine_agreement $datdir/xmax_only/epos/xmax-only/plots/lna_composition_results_TreeS6_xmax.root > $datdir/xmax_only/epos/xmax-only/plots/agreement_anal-pub_printout.txt

  # Agreement between PAO data and published data
  echo "Checking agreement between PAO data and published data (relative)"
  ./determine_agreement $datdir/mva_analysis/epos/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS6.root > $datdir/mva_analysis/epos/xmax_deltas38_deltarise/plots/agreement_anal-pub_printout.txt

  # Agreement between PAO data and published data
  echo "Checking agreement between PAO data and published data (absolute)"
  ./determine_agreement $datdir/mva_analysis/epos/xmax_s1000_rise_zenithFD/plots/lna_composition_results_TreeS6.root > $datdir/mva_analysis/epos/xmax_s1000_rise_zenithFD/plots/agreement_anal-pub_printout.txt

  # Agreement between AugerMix and published data
  echo "Checking agreement between AugerMix and published data (FD-only)"
  ./determine_agreement $datdir/xmax_only/epos/xmax-only/plots/lna_composition_results_TreeS5_xmax.root > $datdir/xmax_only/epos/xmax-only/plots/agreement_anal-pub_AugerMix_printout.txt

  # Agreement between AugerMix and published data
  echo "Checking agreement between AugerMix and published data (relative)"
  ./determine_agreement $datdir/mva_analysis/epos/xmax_deltas38_deltarise/plots/lna_composition_results_TreeS5.root > $datdir/mva_analysis/epos/xmax_deltas38_deltarise/plots/agreement_anal-pub_AugerMix_printout.txt

  # Agreement between AugerMix and published data
  echo "Checking agreement between AugerMix and published data (absolute)"
  ./determine_agreement $datdir/mva_analysis/epos/xmax_s1000_rise_zenithFD/plots/lna_composition_results_TreeS5.root > $datdir/mva_analysis/epos/xmax_s1000_rise_zenithFD/plots/agreement_anal-pub_AugerMix_printout.txt
