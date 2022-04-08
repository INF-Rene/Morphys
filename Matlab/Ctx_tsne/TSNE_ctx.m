%%
t=readtable('/Users/elinemertens/Data/Projects/Hippocampus/Data/Excel/AP_property_TSNE/Ctx_morph_try1.csv')

mat_feat = t{:,:};

%% actual TSNE mapping, check these features and its explanation on
% mathworks help center TSNE
tm = tsne(mat_feat,'Distance', 'euclidean', 'Perplexity', 10,...
'Standardize', true);

%%
close all 
%%
% In feature you can change the feature you want to base color coding on
feature='SegmentLength';
c = mat_feat(:,strcmp(t.Properties.VariableNames, feature));

% here you choose how to do the color coding and add the legend
scatter(tm(:,1), tm(:,2), [], c, 'filled')
colorbar
colormap jet

%% here you can find the specific values of interest and see which one they
% are
find(tsne_mapping(:,1)<-100) %smaller than -100 on x axis
find(tsne_mapping(:,1)>80) %bigger than 100 

%% y axis cut off? to do  
find(tsne_mapping(:,1)<-100)