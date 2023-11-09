%%
clear, clc, close all
%% Load in csv excel file 
fn = ('/Users/elinemertens/Data/Projects/Ch3.Method section/Data/cluster analysis/T_L2L3.csv');
opts = detectImportOptions(fn, "VariableNamingRule", "preserve");
charVars = contains(opts.VariableNames, 'char');
opts.VariableTypes(charVars) = {'double'};
t = readtable(fn, opts);

%% My excel contains comma instead of dots, replace here 
% for morph necessary cell_feat = num2cell(t(:,7) ; 
mat_feat = t{:,:};
    %cell_feat=cellfun(@(x) str2num(strrep(x, ',', '.')), cell_feat, 'UniformOutput', false);
%% mat_feat = cell2mat(cell_feat);
patients = readcell('/Users/elinemertens/Data/Projects/Ch3.Method section/Data/cluster analysis/Patients.csv');
% delete the second row of patients before continuing 

%% alternatively color on hierarchical clusters isntead of groups
N_components= 6;          % set the number of principal components to include
n_clusters = 4;              % set the number of clusters you want 
N = length(mat_feat);        % number of observations from your file

[coeff, score, latent] = pca(normalize(mat_feat));
Z = linkage(normalize(score(:, 1:N_components)),'ward'); %not sure whether to re-normalize the scores here
cutoff = median([Z(end-n_clusters+1,3) Z(end-n_clusters+2, 3)]);
h = dendrogram(Z,N, 'ColorThreshold', cutoff); % h contains Line objects
%
hierclusters = cluster(Z,'Maxclust',n_clusters);
%
Explainedpercentage = latent/sum(latent) ;

%%
Summary_coeff = array2table(coeff) ;
writetable(Summary_coeff, 'summary_coeff.xlsx');

%% actual TSNE mapping, check these features and its explanation on
% mathworks help center TSNE

rng(2) % set random seed for reproducibility (3 works for tsne hipp)
tm = tsne(mat_feat,'Distance', 'euclidean', 'Perplexity', 10,...
'Standardize', true, 'NumPCAComponents', 7); % num pca isn't that 7???

% based on group (patient, cell-type or layer)
%c = [0 1 1; 1 0 0];
% c = [0 0 1; 1 0 0; 0 1 0];
 c = [1 0 1; 1 0 0; 0 1 1; 0 1 0 ; 1 1 0];

% colors based on hierarchical clusters:
figure()
gscatter(tm(:,1),tm(:,2),hierclusters,c, ".", 15, "on")

legend('Location','northeastoutside')

%% 
%c = [0 1 1; 1 0 0];
% c = [0 0 1; 1 0 0; 0 1 0];
 c = [1 0 1; 0 1 1; 0 1 0; 1 0 0];

mkrs={'s','d','o', 'h'};  % n=53 hipp 1

unpats=unique(patients);
for i=1:numel(unpats)
    locs = strcmp(unpats{i}, patients);
    clusternr = hierclusters(locs);
    scatter(tm(locs,1),tm(locs,2), 70, c(clusternr,:), 'filled', mkrs{i})
    ylim([-250 250]) ;
    legend('P117', 'P198', 'P209')
    
    legend('Location','northeastoutside')
    hold on
end
%% make a heatmap of coeff table
tbl = readtable('C:\Users\femke\OneDrive\Documenten\CNCR\Project Hippocampus\Clusteranalysis\heatmap parameters PCA.csv');

cdata = [-0.288	0.436	-0.103	-0.456	0.689	0.041	-0.180; 0.532	0.312	-0.268	0.045	0.204	-0.177	0.688; 
    0.332	0.365	-0.540	-0.217	-0.443	0.102	-0.460; 0.457	-0.239	0.364	-0.485	0.049	-0.558	-0.229; 
    -0.005	0.463	0.602	-0.326	-0.376	0.338	0.248; 0.248	0.452	0.353	0.633	0.195	-0.131	-0.392;
    0.504	-0.322	0.081	-0.040	0.325	0.717	-0.122];
xvalues = {'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7',};
yvalues = {'soma-SR', 'TDL total', '#nodes basal', '#obliques', '#bifurcations trunk', 'sag', 'downstroke'};

h = heatmap(xvalues,yvalues,cdata,"Colormap",autumn,"ColorScaling","scaledcolumns");

h.YLabel = 'parameter';







%% Heatplot
% In feature you can change the feature you want to base color coding on

i = 1;  % Index to create different plot positions
for featuresheatplot = {'soma-SR', 'TDL total', '#nodes basal', '#bifurcations trunk', 'sag', 'downstroke', '#obliques'}
c = mat_feat(:,strcmp(t.Properties.VariableNames, featuresheatplot));
figure('Position',[50*i,40+20*i,640,480])  % Create new figure

scatter(tm(:,1), tm(:,2), [], c, 'filled')
colorbar
colormap hsv
title(featuresheatplot)
end

%% make boxplot
i = 1;  % Index to create different plot positions

for featuresboxplot = {'soma-SR', 'TDL total', '#nodes basal', '#bifurcations trunk', 'sag', 'downstroke', '#obliques'}
    %{'resonance', '3db cut off', 'input resistance', 'sag ratio', 'threshold', 'amplitude', 'peak', 'halfwidth', 'upstroke', 'downstroke', 'up down ratio', 'max upstroke', 'max downstroke', 'rheobase', 'TDL basals', '#basals', '#nodes', 'TDL apical', 'TDL total', 'TDL obl', '#obliques', '#nodes obl', 'TDL trunk', '#trunks', '#bifurcations trunk', 'So-soma', 'soma-SLM_NEW', 'SR', 'SP', 'SP SR', 'so-soma/SP'}

    %{'sag ratio', '#basals', 'TDL apical', '#nodes obl', 'so-soma/SP'}
%'resonance', '3db cut off', 'input resistance', 'sag ratio', 'threshold', 'amplitude', 'peak', 'halfwidth', 'upstroke', 'downstroke', 'up down ratio', 'max upstroke', 'max downstroke', 'rheobase', 'TDL basals', '#basals', '#nodes', 'TDL apical', 'TDL total', 'TDL obl', '#obliques', '#nodes obl', 'TDL trunk', '#trunks', '#bifurcations trunk', 'TDL tuft', '#tufts', 'nodes', 'So-soma', 'soma-SLM_NEW', 'SR', 'SP', 'SP SR', 'fraction location in SP+SR', 'inverted location %' 

data = mat_feat(:,strcmp(t.Properties.VariableNames, featuresboxplot));
figure('Position',[50*i,40+20*i,250,400])  % Create new figure
boxplot(data, hierclusters, "BoxStyle","outline", "Colors","k", "Notch","on")
xlabel('Cluster number');
ylabel(featuresboxplot);
title (featuresboxplot)
axis tight

outliers = findobj( 'Tag', 'Outliers' );
delete(outliers)

%show datapoints in boxplot    
Cluster1 = data(hierclusters==1);
Cluster2 = data(hierclusters==2);
Cluster3 = data(hierclusters==3);
% Cluster4 = data(hierclusters==4);
% Cluster5 = data(hierclusters==5);
% Cluster6 = data(hierclusters==6);
% Cluster7 = data(hierclusters==7);
       
hold on

scatter(ones(size(Cluster1)).*(1+(rand(size(Cluster1))-0.5)/10),Cluster1, 20, 'blue','o')
scatter(ones(size(Cluster2)).*(2+(rand(size(Cluster2))-0.5)/10),Cluster2, 20, 'red','o')
scatter(ones(size(Cluster3)).*(3+(rand(size(Cluster3))-0.5)/10),Cluster3, 20, 'green','o')
% scatter(ones(size(Cluster4)).*(4+(rand(size(Cluster4))-0.5)/10),Cluster4, 20, 'red','o')
% scatter(ones(size(Cluster5)).*(5+(rand(size(Cluster5))-0.5)/10),Cluster5, 20, 'magenta','o')
% scatter(ones(size(Cluster6)).*(6+(rand(size(Cluster6))-0.5)/10),Cluster6, 20, 'yellow','o')
% scatter(ones(size(Cluster7)).*(7+(rand(size(Cluster7))-0.5)/10),Cluster7, 20, 'blue','o')


i = i + 1;
end

%% boxplots morph corresponding to ephys
fn = ('C:\Users\femke\OneDrive\Documenten\CNCR\Project Hippocampus\Clusteranalysis\Hipp_cluster_analysis_ephys_morph_1_n=31.csv');
opts = detectImportOptions(fn, "VariableNamingRule", "preserve");
charVars = contains(opts.VariableNames, 'char');
opts.VariableTypes(charVars) = {'double'};
t = readtable(fn, opts);

mat_feat = t{:,:};


i = 1;  % Index to create different plot positions
for featuresboxplot = { 'resonance', 'input resistance', 'sag', 'threshold', 'amplitude', 'peak', 'halfwidth', 'upstroke', 'downstroke', 'up down ratio', 'max upstroke', 'max downstroke', 'rheobase', 'TDL total', 'TDL basals', '#basals', '#nodes basal', 'nodes per branch basal', 'TDL apical', 'TDL obl', '#obliques', '#nodes obl', 'TDL trunk', '#bifurcations trunk', 'SO-soma', 'soma-SLM', 'soma-SR', 'SR', 'SP', 'SP SR', 'fraction location in SP+SR', 'inverted location %', 'fraction in SP'} 


data = mat_feat(:,strcmp(t.Properties.VariableNames, featuresboxplot));
clustdata = readmatrix('C:\Users\femke\OneDrive\Documenten\CNCR\Project Hippocampus\Clusteranalysis\hierclusters_n=31_ephys_morph_2.csv');
figure('Position',[50*i,40+20*i,250,400])  % Create new figure
boxplot(data, clustdata, "BoxStyle","outline", "Colors","k", "Notch","on")
xlabel('Cluster number');
ylabel(featuresboxplot);
title (featuresboxplot)

outliers = findobj( 'Tag', 'Outliers' );
delete( outliers );

Cluster1 = data(clustdata==1);
Cluster2 = data(clustdata==2);
Cluster3 = data(clustdata==3);
% Cluster4 = data(clustdata==4);
% Cluster5 = data(clustdata==5);
% Cluster6 = data(clustdata==6);
% Cluster7 = data(clustdata==7);
       
hold on

scatter(ones(size(Cluster1)).*(1+(rand(size(Cluster1))-0.5)/10),Cluster1, 20, 'green','o')
scatter(ones(size(Cluster2)).*(2+(rand(size(Cluster2))-0.5)/10),Cluster2, 20, 'red','o')
scatter(ones(size(Cluster3)).*(3+(rand(size(Cluster3))-0.5)/10),Cluster3, 20, 'blue','o')
% scatter(ones(size(Cluster4)).*(4+(rand(size(Cluster4))-0.5)/10),Cluster4, 20, 'yellow','o')
% scatter(ones(size(Cluster5)).*(5+(rand(size(Cluster5))-0.5)/10),Cluster5, 20, 'magenta','o')
% scatter(ones(size(Cluster6)).*(6+(rand(size(Cluster6))-0.5)/10),Cluster6, 20, 'green','o')
% scatter(ones(size(Cluster7)).*(7+(rand(size(Cluster7))-0.5)/10),Cluster7, 20, 'blue','o')


i = i + 1;
end
