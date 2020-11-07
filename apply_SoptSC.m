% Using SoptSC to qualitatively characterize cell-cell communications
%Wang, Shuxiong, et al. "Cell lineage and communication network inference 
%via optimization for single-cell transcriptomics." Nucleic acids research 
%47.11 (2019): e66-e66.

%TGFB signaling pathway
%load SCC gene expression matrix with clustring result
addpath('data');
load('SCC_SoptSC.mat')
%ligands
Lig = {'Tgfb1','Tgfb2','Tgfb3','Tgfb1','Tgfb2','Tgfb3','Tgfb1','Tgfb2','Tgfb3','Tgfb1','Tgfb2','Tgfb3'};
%rescprtors
Rec = {'Tgfbr1','Tgfbr1','Tgfbr1','Tgfbr2','Tgfbr2','Tgfbr2','Acvr1','Acvr1','Acvr1','Acvr1b','Acvr1b','Acvr1b','Acvr1c','Acvr1c','Acvr1c'};
%targer genes
target_up = {'Fn1','Vtn','Cdh2','Col1a1','Col1a2','Mmp2','Mmp3','Mmp9'...
    ,'Twist1','Twist2','Ids','Zeb1','Zeb2'...
    ,'Sparc','Itga5','Itgb3','Ncam','Vim','Acta2','Plau','Dab2','Hic5','Tgfb1i1','Hmga2'};
target_down = {'Ocln','Cldn1','Cldn2','Cldn3','Cldn4','Cldn5','Cldn6','Cldn7','Cldn8'...
    ,'Cldn9','Cldn10','Cldn11','Cldn12','Cldn13','Cldn14','Cldn15','Cldn16'...
    ,'Cldn17','Cldn18','Cldn19','Cldn20','Cldn21','Cldn22','Cldn23','Cd34'...
    ,'Cdh1','Dsp','Pkp1','Pkp2','Pkp3','Crb3','Ck5','Ck14','Ck8','Ck18'...
    ,'Crb3','Esr1'};
%find the above genes in the entire gene set
%use Matlab code of SoptSC
addpath('Cell_Cell_Communication');
[Lig,Rec] = ligand_receptor_in_list(allgenes,Lig,Rec);
target_up = intersect(allgenes,target_up,'stable');
target_down = intersect(allgenes,target_down,'stable');

%compute cell-cell & cluster_cluster interactions
[P_cell_each_lr,P_cell_agg,P_cluster_each_lr,P_cluster_agg] = ...
    LR_Interaction(data_raw_log',allgenes,cell_clustering_index,Lig,Rec,target_up,target_down);
