% Multilayer regulations of EMT

%load SCC datasets with gene-gene correlations computed by PIDC
addpath('data');
load('SCC_multi_GRN.mat')
multilayer_GRN(P_cluster_agg,No_cluster,cluster_order1,target_genes,...
    target_up_index,target_down_index,PIDC_target,PIDC_all_genes,...
    data,gene_name,marker_genes,C,marker_genes_index)