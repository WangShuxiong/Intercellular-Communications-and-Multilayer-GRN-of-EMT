% Multilayer regulations of EMT

%load SCC datasets with gene-gene correlations computed by PIDC
addpath('data');
load('SCC_multi_GRN.mat')
%Input: P_cluster_agg (cluster-cluster signaling probabilities)
%       No_cluster (# cluster)
%       target_genes (index of target genes)
%       target_up_index (index of up-regulated target genes)
%       target_down_index (index of down-regulated target genes)
%       cluster_order1 (order of inferred clusters)
%       data (log(gene expression matri+1))
%       PIDC_target (correlations matrics of target genes computed by PIDC)
%       PIDC_all_genes (correlations matrics of target&marker genes computed by PIDC)
%       gene_names (names of genes)
%       marker_genes (index of marker genes)
%       C (clustering result of cells)
%       marker_genes_index (cluster of each marker)

multilayer_GRN(P_cluster_agg,No_cluster,cluster_order1,target_genes,...
    target_up_index,target_down_index,PIDC_target,PIDC_all_genes,...
    data,gene_name,marker_genes,C,marker_genes_index)
