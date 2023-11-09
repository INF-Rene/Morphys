% Load in csv excel file 
fn = ('/Users/elinemertens/Data/Projects/Collaborations/EPFL/EPFL_conf_2022/Ctx_L5_ephys_features_tsne.csv');
opts = detectImportOptions(fn, "VariableNamingRule", "preserve");
charVars = contains(opts.VariableNames, 'char');
opts.VariableTypes(charVars) = {'double'};
t = readtable(fn, opts);

%%
t=readtable('C:\Users\femke\OneDrive\Documenten\Thesis\Clusteranalysis\human_cluster_analysis_mastersheet.csv')

%% My excel contains comma instead of dots, replace here 
% for morph necessary cell_feat = num2cell(t(:,7) ; 
mat_feat = t{:,:};
%cell_feat=cellfun(@(x) str2num(strrep(x, ',', '.')), cell_feat, 'UniformOutput', false);
%mat_feat = cell2mat(cell_feat);
patients = readcell('C:\Users\femke\OneDrive\Documenten\Thesis\Clusteranalysis\human_group_parameter_cell_type.csv');
% delete the second row of patients before continuing 

%% actual TSNE mapping, check these features and its explanation on
% mathworks help center TSNE
tm = tsne(mat_feat,'Distance', 'euclidean', 'Perplexity', 10,...
'Standardize', true);


%% based on group (patient, cell-type or layer)
% color the data based on parameters
% choose the feature you want to color the data 
feature='Threshold';
cdata = mat_feat(:,strcmp(t.Properties.VariableNames, feature));
c = jet(numel(cdata));
[sorted_data, idx] = sort(cdata);
c = c(idx, :);

%%

% give symbols based on unique groups
mkrs={'d', 'o', 's', 'p', '>', 'h'};  % 6 markers for cell-types
% mkrs={'d', 'o', 'h'};                 % 3 markers for layers
% mkrs={'d', 'o', 'h'};                   % 3 markers for species
unpats=unique(patients);
for i=1:numel(unpats)
    locs = strcmp(unpats{i}, patients);
    scatter(tm(locs,1),tm(locs,2), 40, c(locs,:), 'filled', mkrs{i})
    legend ('Exc L2 LAMP5 LTK', 'Exc L2-3 LINC00507 FREM3', 'Exc L2-4 LINC00507 GLP2R', ...
        'Exc L3-4 RORB CARM1P1','Exc L3-5 RORB COL22A1', 'noname')      % legend titles for cell-types
    % legend('L2/3', 'L4', 'L5/6')      % legend titles for layers
    % legend('human', 'mouse TeA', 'rat S1')     % legend titels for species
    hold on
end

% color bar
caxis([min(cdata), max(cdata)])
cmap = colormap(jet);
cbar = colorbar();

% gscatter(tm(:,1),tm(:,2),patients)
%%
close all 



%%
% In feature you can change the feature you want to base color coding on
     feature='InputResistance';
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