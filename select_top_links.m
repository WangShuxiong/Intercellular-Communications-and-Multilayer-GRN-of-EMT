function [PIDC_all_matrix] = select_top_links(PIDC_all_matrix,top_p)
a = sort(PIDC_all_matrix(:),'descend');
%top_p = 0.14;
top_cut = a(round(length(a)*top_p));
indices = find(PIDC_all_matrix<top_cut);
PIDC_all_matrix(indices)=0;
end

