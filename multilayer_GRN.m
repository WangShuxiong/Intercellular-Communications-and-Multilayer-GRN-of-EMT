function [] = multilayer_GRN(P_cluster_agg,No_cluster,cluster_order1,target_genes,...
    target_up_index,target_down_index,PIDC_target,PIDC_all_genes,...
    data,gene_name,marker_genes,C,marker_genes_index)

%Input: P_cluster_agg (cluster-clustet signaling probabilities)
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

% plot first cluster-cluster layer

threshold = 0.5; %only plot signling probability >0.5
zzz = get(gca,'colororder');
mycolor = zeros(11,3);
mycolor(1:7,:) = zzz;
mycolor(1,:) = 	[248, 118, 109]/255;
mycolor(2,:) = 	[124, 174, 0]/255;
mycolor(3,:) = 	[0, 191, 196]/255;
mycolor(4,:) = 	[199, 124, 255]/255;
mycolor(8,:) = [0 0 0];
mycolor(9,:) = [0 0 0.803922];
mycolor(10,:) = [1 0 1];
mycolor(11,:) = [0.5 1 0];
mycolor1 = mycolor(1:No_cluster,:);

adjacentM =P_cluster_agg;
adjacentM = adjacentM./max(adjacentM(:));
adjacentM(adjacentM < threshold) = 0;

trans_p= 0.85; %transparency
m_size = 7;
if max(adjacentM(:)) > 0
    bg = digraph(adjacentM);%,'omitselfloops');
    bg.Edges.LWidths = 2*bg.Edges.Weight/max(bg.Edges.Weight);
    
    % Set my edge color
    aa = bg.Edges.EndNodes;
    Gh = plot(bg,'Layout','layered','Sources',1:4, 'Marker','o','MarkerSize',30,...
        'NodeColor',mycolor1,'EdgeColor',[0  0  0],'LineWidth',...
        bg.Edges.LWidths,'ArrowSize',0,'EdgeAlpha',trans_p);
    Gh.NodeLabel=[];
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    axis off;
end

hold on
% layer
a=0.3;b=1;
h1x = [Gh.XData(1)-b, Gh.XData(end)+b,Gh.XData(end)+b,Gh.XData(1)-b];
h1y = [Gh.YData(1)-a,Gh.YData(1)-a,Gh.YData(1)+a,Gh.YData(1)+a ];
s = fill3(h1x,h1y,zeros(1,No_cluster),[0.4 0 0],'LineStyle','none');
alpha(s,.1)
hold on


% second layer target genes
% link between first & second layer
M_taget = zeros(No_cluster, length(target_genes));
for i = 1:No_cluster
    for j = 1:length(target_genes)
        M_taget(i,j) = mean(data(C==cluster_order1(i),target_genes(j)));
    end
end
M_taget = M_taget./max(max(M_taget));
%link >0.5
M_taget = select_top_links(M_taget,0.2);
[row1,col1] = find(M_taget>0);


%option 1 plot GRN of target genes within clusters

for i = 1:No_cluster
    PIDC_target1 = PIDC_target{i};
    G = graph(PIDC_target1,'OmitSelfLoops');
    h2 = plot(G,'EdgeColor',[0  0  0 ],'LineWidth',1.5,'MarkerSize',m_size,'EdgeAlpha',trans_p);
    hold on
    h2.NodeLabel=[];
    xd = get(h2, 'XData');
    yd = get(h2, 'YData');
    h2.XData = h2.XData/6.5+Gh.XData(i);
    h2.XData = h2.XData - mean(h2.XData) + Gh.XData(i);
    h2.YData = (h2.YData/7.5+Gh.YData(i))/2;
    h2.YData = h2.YData-mean(h2.YData)+Gh.YData(1);
    h2.ZData = -5*ones(1,size(h2.ZData,2));
    %text(h2.XData, h2.YData,h2.ZData, cellstr(gene_name(target_genes)), 'FontSize',10, 'FontWeight','bold', 'HorizontalAlignment','left', 'FontName','Arial','VerticalAlignment','top')

    highlight(h2,find(ismember(target_genes,target_up_index)),'NodeColor',mycolor(6,:)) %up
    highlight(h2,find(ismember(target_genes,target_down_index)),'NodeColor',mycolor(7,:)) %down
    

    %link between layers
    index = find(row1==i);
    for j = 1:length(index)
        p1=plot3([Gh.XData(row1(index(j))) h2.XData(col1(index(j)))],[Gh.YData(row1(index(j))) h2.YData(col1(index(j)))],[Gh.ZData(row1(index(j))) h2.ZData(col1(index(j)))],'Color',[0  0  0 ],'LineWidth',2*M_taget(row1(index(j)),col1(index(j))));
        p1.Color(4) = trans_p;
        hold on
    end
    


    s = fill3(h1x,h1y,-5*ones(size(h1x)),[0 0 0.4],'LineStyle','none');
    alpha(s,.1)
    hold on
    
% third layer includes marker&trans genes
    PIDC_all_genes1 = PIDC_all_genes{i};
    PIDC_mt = PIDC_all_genes1(1:length(marker_genes),1:length(marker_genes));
    G = graph(PIDC_mt,'OmitSelfLoops');
    h3 = plot(G,'EdgeColor',[0  0  0 ],'LineWidth',1.5,'MarkerSize',m_size,'EdgeAlpha',trans_p);
    hold on
    h3.NodeLabel=[];
    xd = get(h3, 'XData');
    yd = get(h3, 'YData');
    h3.XData = h3.XData/5.5+Gh.XData(i);
    h3.XData = h3.XData - mean(h3.XData) + Gh.XData(i);
    h3.YData = (h3.YData/7+Gh.YData(i))/2;
    h3.YData = h3.YData-mean(h3.YData)+Gh.YData(1);
    h3.ZData = -10*ones(1,size(h3.ZData,2));
    %text(h3.XData, h3.YData, h3.ZData,cellstr(gene_name(marker_genes)), 'FontSize',10, 'FontWeight','bold', 'HorizontalAlignment','left', 'FontName','Arial','VerticalAlignment','top')
    highlight(h3,find(ismember([marker_genes],marker_genes(marker_genes_index==1))),'NodeColor',mycolor(1,:)) %E
    highlight(h3,find(ismember([marker_genes],marker_genes(marker_genes_index==2))),'NodeColor',mycolor(2,:)) %I1
    highlight(h3,find(ismember([marker_genes],marker_genes(marker_genes_index==3))),'NodeColor',mycolor(3,:)) %I2
    highlight(h3,find(ismember([marker_genes],marker_genes(marker_genes_index==4))),'NodeColor',mycolor(4,:)) %M
    
    %link between second & third layer
    M_trans = PIDC_all_genes1(length(marker_genes)+1:end,1:length(marker_genes));
    M_trans = M_trans./max(max(M_trans));
    % only select top 1.5% of links
    M_trans = select_top_links(M_trans,0.015);
    [row,col] = find(M_trans>0);
    for j = 1:length(row)
        p1=plot3([h2.XData(row((j))) h3.XData(col((j)))],[h2.YData(row((j))) h3.YData(col((j)))],[h2.ZData(row((j))) h3.ZData(col((j)))],'Color',[0  0  0 ],'LineWidth',2*M_trans(row((j)),col((j))));
        p1.Color(4) = trans_p;
        hold on
    end
    
    s = fill3(h1x,h1y,-10*ones(size(h1x)),[0 0.4 0],'LineStyle','none');
    alpha(s,.1)
    hold on
end

end

