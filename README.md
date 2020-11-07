# Intercellular-Communications-and-Multilayer-GRN-of-EMT

Walkthrough: Mouse skin squamous cell carcinoma (SCC dataset): This dataset of 382 cells on skin tumors contains FACS_isolated epithelial YFP+Epcam+ tumor cells, which are relatively homogeneous, and mesenchymal-like YFP+Epcam- tumor cells, which are more heterogeneous.

Usage: For matlab codes, download the source codes and unzip the MATLAB package. Change the current directory in MATLAB to the folder containing the scripts.

1) [QuanTC](https://github.com/yutongo/QuanTC/blob/master/Example/QuanTC_SCC.pdf): clustering and transition trajectory reconstruction. 
2) [PIDC](https://github.com/Tchanders/network_inference_tutorials): using partial information decomposition to identify GRN 
3) SoptSC: qualitatively characterizing cell-cell communications [code](https://github.com/yutongo/Intercellular-Communications-and-Multilayer-GRN-of-EMT/blob/main/apply_SoptSC.m)
4) Multilayer regulations of EMT [code](https://github.com/yutongo/Intercellular-Communications-and-Multilayer-GRN-of-EMT/blob/main/plot_multilayer_GRN.m)
5) [igraph](https://igraph.org/r/doc/closeness.html): measuring node centrality [code](https://github.com/yutongo/Intercellular-Communications-and-Multilayer-GRN-of-EMT/blob/main/apply_igraph.R)
