# Metacognition
Codes for "Macaque Prefrontal Cortex Integrates Multiple Components for Metacognitive Judgments of Working Memory"



The "BehavioralAnalyses" folder contains behavioral analyses corresponding to Figure 1 and pupil size analyses corresponding to Figure 5.

   * This file already contains some demo behaviral data to run the code.

   * run "test_temp_plot_metaEvidence.m" for reproducing behavioral analyses, related to Figure 1.

   * run "pupilSizeAnalysis_baseline_2p.m" for pupil size analyses, related to Figure 5.   
    



The "TwoPhotonDataAnalyses" folder includes my two-photon preprocessing pipeline (motion correction with NoRMCorre and source extraction with Suite2p) and comprehensive neuronal data analyses.

   * This folder does not contain the complete neuronal dataset required to run the code. The full dataset is available via a dataset-specific Zenodo link. This folder only includes data from a single FOV named as "demoData.mat" (from monkey D FOV 8, with two consecutive sessions merged). If the variable if_loadDemoData is set to 1, most of the results in "my2pScripts" can still be plotted.

   * NoRMCorre: my pipeline for motion correction with NoRMCorre

   * forSuite2p: my pipeline for source extraction with Suite2p

   * my2pScripts: comprehensive neuronal data analyses

   1. In my2pScripts, these files corresponds to example neurons, related to Figure 1G-1I:  
      "test_temp_plot_locTuning_newFigure_singleNeuron.m";  
      "test_temp_plot_metaTuning_newFigure_singleNeuron.m";

   2. In my2pScripts, these files corresponds to WM neural decoder, related to Figure 2:  
      "stepD3_train123_test123.m";

   3. In my2pScripts, these files corresponds to (baseline) meta-WM neural decoder, related to Figure 3 and Figure 5.  
      "stepDM1_test_decodingMeta.m";

   4. In my2pScripts, these files corresponds to decoding time courses of WM and meta-WM, related to Figure 4D:  
      "stepF5D_test_mismatchMechnism_meta_crossTime.m";  
      "stepF5E_test_mismatchMechnism_memoryPrecision_crossTime.m";

   5. In my2pScripts, these files corresponds to validation of trial proportions for low-strength mismatch in memory-error trials, related to Figure 4H:  
      "stepF5A_memoryMetaMismatch_twoDecoder.m";

   6. In my2pScripts, these files corresponds to meta-WM linear regression, related to Figure 5G-5I:  
      "test_temp_baselineMeta_meta_GLM.m"

   7. In my2pScripts, these files corresponds to neeronal spatial organization analyses, related to Figure 6F-6H:  
      "stepF4B_test_plot_FOV_wholeImage_complexTuning.m" (single-FOV);  
      "test_temp_plot_mulFOV_singleNeuron.m" (multi-FOVs);


   8. In my2pScripts, these files corresponds to conduct multi-FOVs analyses pipeline:  
      "stepZ1_summaryMultiFov.m";

   9. In my2pScripts, these files corresponds to plot multi-FOVs analyses results from above:  
      "test_temp_plot_mulFOV_singleNeuron.m";

   10. In my2pScripts, these files contain code to run monkey behavioral experiments by PsychtoolBoox in Matlab:  
      "trainDing_Step15_eye_for113Record_niMarker_56.m";  
      "trainZelku_Step15_eye_for113Record_niMarker_84.m";            
