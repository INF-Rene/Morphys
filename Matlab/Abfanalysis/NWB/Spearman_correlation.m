% Ephys spearman 

% voeg de namen toe van de parameters die je gebruikt
%Noem deze data_cell_table

% Load the table from the specified directory
input_directory = '/Users/elinemertens/Data/Projects/Ch3.Method section/Figures/Fig4/';
input_filename = 'data_cell_table.csv';
input_path = fullfile(input_directory, input_filename);
loaded_table = readtable(input_path);
%%
data_cell_table(data_cell_table == 0) = NaN;

%% With a new / adjusted table, make sure all of these are specified
%Specify column names
columnNames = {'RelativeLoc','Tau', 'Sagratio', 'Rinput', 'Rheobase' , 'Threshold', 'Amplitude', 'Halfwidth', 'Upstroke', 'Downstroke', 'updownratio', 'Adaptindex', 'ISIratio','f-Islope' 'Fres'};
% Create a new table with column titles
data_cell_table = array2table(data_cell_table, 'VariableNames', columnNames);

%% Replace NaN values with column means in a cell array
% no longer needed!!!
for col = 1:size(data_cell_temp, 2)
    colData = data_cell_temp(:, col);
    nanIndices = isnan(colData);
    colMean = mean(colData, 'omitnan');
    colData(nanIndices) = colMean;
    data_cell_temp(:, col) = colData;
end

%% Calculate Spearman correlation
rho = corr(data_cell_table{:,:}, 'Type', 'Spearman', 'rows', 'pairwise');
[rho, p] = corr(data_cell_table{:,:}, 'Type', 'Spearman', 'rows', 'pairwise');
%% multiple corrections, p value correlation
% Number of tests
num_tests = numel(p);
% Apply the Bonferroni correction
alpha_adjusted = 0.05 / ((num_tests / 2 ) -15);
% Identify significant correlations based on the adjusted significance level
significant_correlations = p < alpha_adjusted;
% Reshape the boolean vector back into a matrix of the same size as the original correlation matrix
smp = reshape(significant_correlations, size(p));
% Create a custom colormap for black (0) and white (1)
custom_colormap = [1, 1, 1; 0, 0, 0];

% Plot the significant matrix using imagesc
figure(1);
smatrixp = tril(smp);
smatrixp(logical(eye(size(smatrixp)))) = 0;
imagesc(smatrixp);

% Apply the custom colormap
colormap(custom_colormap);
title('Custom Colormap for P-Values Heatmap');
xticks(1:22);
yticks(1:22);
caxis([0, 0.5]); 
set(gca, 'TickDir', 'out')
xticklabels(columnNames);
yticklabels(columnNames);
xlabel('Variables');
ylabel('Variables');

%% Display or use the correlation matrix 'rho'
disp(rho); 
% rho = corr(data_cell, 'Type', 'Spearman');
% [rho, p] = corr(data_cell, 'Type', 'Spearman');

%% Plot the rho colormap with the custom colormap
%hier laad je de colormap in, custommie, staan bij ch3, fig4. deze bevatten
%de gradienten 
% Extract the lower triangular part of the matrix
figure(2)
lower_triangle = tril(rho);
lower_triangle(logical(eye(size(lower_triangle)))) = 0;
imagesc(lower_triangle);
colormap(custommie);  % Use the custom colormap directly
colorbar;
caxis([-1, 1]);  % Adjust the range as neededx
title('Custom Colormap for Rho Heatmap');
xticks(1:15);
yticks(1:15);
 set(gca, 'TickDir', 'out')
caxis([-1, 1]); 
xticklabels(columnNames);
yticklabels(columnNames);
xlabel('Variables');
ylabel('Variables');


%%
% Extract the data for the variables of interest
relativeloc = data_cell_table.RelativeLoc;
%nrAPsTS2 = data_cell_table.nrAPsTS2;
Rinput = data_cell_table.Rinput;
sagratio = data_cell_table.Sagratio;
halfwidth = data_cell_table.Halfw ; 
downstroke = data_cell_table.Downstroke ;
upstroke = data_cell_table.Upstroke ;
amplitude = data_cell_table.Amp ;
Fres = data_cell_table.Fres ; 
% [h, p, jbstat, cval] = jbtest(relativeloc);
% if h
%     disp('The data is not normally distributed.');
% else
%     disp('The data is normally distributed.');
% end

% Calculate Spearman correlation coefficients and p-values
%[rho_nrAPsTS2, p_nrAPsTS2] = corr(relativeloc, nrAPsTS2, 'Type', 'Spearman');
[rho_Rinput, p_Rinput] = corr(relativeloc, Rinput, 'Type', 'Spearman','rows', 'pairwise');
[rho_sagratio, p_sagratio] = corr(relativeloc, sagratio, 'Type', 'Spearman','rows', 'pairwise');
[rho_downstroke, p_downstroke] = corr(halfwidth, downstroke, 'Type', 'Spearman','rows', 'pairwise');
[rho_amplitude, p_amplitude] = corr(upstroke, amplitude, 'Type', 'Spearman','rows', 'pairwise');
[rho_Fres, p_Fres] = corr(sagratio, Fres, 'Type', 'Spearman','rows', 'pairwise');

% Perform linear regression for each correlation
%fit_nrAPsTS2 = polyfit(relativeloc, nrAPsTS2, 1);
fit_Rinput = polyfit(relativeloc, Rinput, 1);
fit_sagratio = polyfit(relativeloc, sagratio, 1);
fit_downstroke = polyfit(halfwidth,downstroke, 1); 
fit_amplitude = polyfit(upstroke, amplitude, 1);

% Create scatter plots and regression lines without confidence intervals
figure(2);

% subplot(2, 2, 1); % Create subplot for nrAPsTS2
% scatter(relativeloc, nrAPsTS2, 'filled');
% hold on;
% plot(relativeloc, polyval(fit_nrAPsTS2, relativeloc), 'r', 'LineWidth', 2);
% xlabel('Relative Loc');
% ylabel('nrAPsTS2');
% title(['Spearman Correlation = ' num2str(rho_nrAPsTS2), ', p = ', num2str(p_nrAPsTS2)]);
% legend('Data', 'Regression Line');

subplot(1, 4, 1); % Create subplot for Rinput
scatter(relativeloc, Rinput, 'filled');
hold on;
set(gca, 'TickDir', 'out')
%plot(relativeloc, polyval(fit_Rinput, relativeloc), 'r', 'LineWidth', 2);
xlabel('RelativeLoc');
ylabel('Rinput (mOhm)');
title(['Spearman Correlation = ' num2str(rho_Rinput), ', p = ', num2str(p_Rinput)]);
%legend('Data', 'Regression Line');

subplot(1, 4, 2); % Create subplot for sagratio
scatter(relativeloc, sagratio, 'filled');
hold on;
set(gca, 'TickDir', 'out')
%plot(relativeloc, polyval(fit_sagratio, relativeloc), 'r', 'LineWidth', 2);
xlabel('RelativeLoc');
ylabel('Sagratio');
title(['Spearman Correlation = ' num2str(rho_sagratio), ', p = ', num2str(p_sagratio)]);
%legend('Data', 'Regression Line');

subplot(1, 4, 3); % Create subplot for halfwidth vs downstroke
scatter(halfwidth, downstroke, 'filled');
hold on;
set(gca, 'TickDir', 'out')
%plot(halfwidth, polyval(fit_downstroke, halfwidth), 'r', 'LineWidth', 2);
xlabel('Halfwidth (ms)');
ylabel('Downstroke (mV/ms)');
title(['Spearman Correlation = ' num2str(rho_downstroke), ', p = ', num2str(p_downstroke)]);
%legend('Data', 'Regression Line');

subplot(1, 4, 4) ; %subplot for fres vs sag
scatter(Fres, sagratio, 'filled')
 hold on;
 set(gca, 'TickDir', 'out')
xlabel('Fres (mV/ms)');
 ylabel('sag ratio (mV)');
 title(['Spearman Correlation = ' num2str(rho_Fres), ', p = ', num2str(p_Fres)]);



% subplot(1, 4, 4); % Create subplot for amplitude
% scatter(upstroke, amplitude, 'filled');
% hold on;
% set(gca, 'TickDir', 'out')
% %plot(upstroke, polyval(fit_amplitude, upstroke), 'r', 'LineWidth', 2);
% xlabel('Upstroke (mV/ms)');
% ylabel('Amplitude (mV)');
% title(['Spearman Correlation = ' num2str(rho_amplitude), ', p = ', num2str(p_amplitude)]);
% %legend('Data', 'Regression Line');

% Check if the correlations are significant
alpha = alpha_adjusted ;% Significance level

% if p_nrAPsTS2 < alpha
%     disp('nrAPsTS2 is significantly correlated with Relative Loc.');
% else
%     disp('nrAPsTS2 is not significantly correlated with Relative Loc.');
% end

if p_Rinput < alpha
    disp('Rinput is significantly correlated with Relative Loc.');
else
    disp('Rinput is not significantly correlated with Relative Loc.');
end

if p_sagratio < alpha
    disp('Sag Ratio is significantly correlated with Relative Loc.');
else
    disp('Sag Ratio is not significantly correlated with Relative Loc.');
end

if p_downstroke < alpha
    disp('Downstroke is significantly correlated with Halfwidth.');
else
    disp('Downstroke is not significantly correlated with halfwidth.');
end

if p_amplitude < alpha
    disp('Amplitude is significantly correlated with upstroke.');
else
    disp('Amplitude is not significantly correlated with upstroke.');
end

%% if you want to extract the n number for each of the correlations
% Assuming your table is named 'data_cell_table'
num_values_per_variable = sum(~isnan(table2array(data_cell_table)));

disp('Number of non-NaN values per variable:');
disp(num_values_per_variable);

% in a matrix it would work like this
data_matrix = table2array(data_cell_table);

n = size(data_matrix, 2);

rho_matrix = zeros(n);
p_matrix = zeros(n);
num_cells_matrix = zeros(n);

for i = 1:n
    for j = 1:n
        if i ~= j
            non_nan_idx = ~isnan(data_matrix(:, i)) & ~isnan(data_matrix(:, j));

            [rho_matrix(i, j), p_matrix(i, j)] = corr(data_matrix(non_nan_idx, i), data_matrix(non_nan_idx, j), 'Type', 'Spearman');
            num_cells_matrix(i, j) = sum(non_nan_idx);
        end
    end
end

disp('Spearman correlation matrix:');
disp(rho_matrix);

disp('P-value matrix:');
disp(p_matrix);

disp('Number of cells matrix:');
disp(num_cells_matrix);


%%
% Save the table to the specified directory
output_directory = '/Users/elinemertens/Data/Projects/Ch3.Method section/Figures/Fig4/';
output_filename = 'data_cell_table.csv';
output_path = fullfile(output_directory, output_filename);
writetable(data_table, output_path);


















