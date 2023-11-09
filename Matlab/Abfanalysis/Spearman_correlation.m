%je maakt een variable aan met de 10 kolommen, met alle waardes voor alle
%cellen die in je database zitten dus 10 in kolommen en in rijen (n = alle
%cellen)

% voeg de namen toe van de parameters die je gebruikt
%Noem deze data_cell

% Load the table from the specified directory
input_directory = '/Users/elinemertens/Data/Projects/Ch3.Method section/Figures/Fig4/';
input_filename = 'data_cell_table.csv';
input_path = fullfile(input_directory, input_filename);
loaded_table = readtable(input_path);

%% With a new / adjusted table, make sure all of these are specified
%Specify column names
columnNames = {'Halfwidth','Downstroke', 'Tau', 'Threshold', 'nrAPsTS2', 'Rinput', 'AHPFAP', 'Upstroke', 'Amplitude', 'sagratio', 'Rheobase', 'ISIratio','relativeloc'};
% Create a new table with column titles
data_cell_table = array2table(data_cell, 'VariableNames', columnNames);

%% Replace NaN values with column means in a cell array
for col = 1:size(data_cell_temp, 2)
    colData = data_cell_temp(:, col);
    nanIndices = isnan(colData);
    colMean = mean(colData, 'omitnan');
    colData(nanIndices) = colMean;
    data_cell_temp(:, col) = colData;
end

%% maak een correlatie plot van 10x10 voor alle parameters met elkaar

rho = corr(data_cell, 'Type', 'Spearman');
[rho, p] = corr(data_cell, 'Type', 'Spearman');

%% % Define the number of colormap steps (you can adjust this value)
nSteps = 64;

% Create a custom colormap for blue
blue_colormap = [linspace(0, 0, nSteps)', linspace(0, 0, nSteps)', linspace(1, 1, nSteps)'];  % Blue

% Create a custom colormap for red
red_colormap = [linspace(1, 1, nSteps)', linspace(0, 0, nSteps)', linspace(0, 0, nSteps)'];  % Red

% Determine the index where the white color starts and ends
white_start_index = round(nSteps / 3);  % Adjust this value for the desired width of the white area
white_end_index = nSteps - white_start_index;

% Create a gradient colormap from blue to white (-1 to -0.2)
gradient_blue_to_white = [linspace(0, 1, white_end_index-white_start_index)', linspace(0, 1, white_end_index-white_start_index)', linspace(1, 1, white_end_index-white_start_index)'];  % Gradient from Blue to White

% Create a gradient colormap from white to red (0.2 to 1)
gradient_white_to_red = [linspace(1, 1, white_end_index-white_start_index)', linspace(0, 0, white_end_index-white_start_index)', linspace(0, 0, white_end_index-white_start_index)'];  % Gradient from White to Red

% Combine the two colormaps for blue and red, and the gradient in the middle
%custom_colormap = [blue_colormap(1:white_start_index, :); gradient_blue_to_white; repmat([1, 1, 1], white_end_index-white_start_index, 1); gradient_white_to_red; red_colormap(1:white_start_index, :)];
%custom_colormap = [blue_colormap(1:white_start_index, :); gradient_blue_to_white; repmat([1, 1, 1], white_end_index-white_start_index, 1); gradient_white_to_red; red_colormap(white_end_index:end, :)];
custom_colormap = [blue_colormap(1:white_start_index, :); gradient_blue_to_white; repmat([1, 1, 1], white_end_index-white_start_index, 1); gradient_white_to_red; red_colormap(white_end_index:end, :)];


% Plot the rho colormap with the custom colormap
imagesc(rho);
colormap(custom_colormap);  % Use the custom colormap directly
colorbar;
caxis([-1, 1]);  % Adjust the range as needed
title('Custom Colormap for Rho Heatmap');
% Add xticks, yticks, labels, and other settings as needed


%% 
figure;
figure(1)
imagesc(rho);
colormap(redblue(45));  % You can choose a different colormap if needed
colorbar;
title('Spearman Rank Correlation Heatmap');
axis square;
xticks(1:13);
yticks(1:13);
 set(gca, 'TickDir', 'out')
caxis([-1, 1]); 
xticklabels(variable_names);
yticklabels(variable_names);
xlabel('Variables');
ylabel('Variables');
%% make a binary colormap
% Define the number of colormap steps (you can adjust this value)
nSteps = 16;

% Threshold for switching to white
threshold = 0.1;  % Corrected threshold value

% Create a custom colormap with white (for values < 0.05) and black (for values >= 0.05)
white_colormap = [linspace(1, 1, nSteps)', linspace(1, 1, nSteps)', linspace(1, 1, nSteps)'];  % White

% Determine the index where the threshold falls
threshold_index = round(threshold * nSteps);

% Create a black colormap for values >= 0.1
black_colormap = [linspace(0, 0, nSteps-threshold_index)', linspace(0, 0, nSteps-threshold_index)', linspace(0, 0, nSteps-threshold_index)'];  % Black

% Combine the two colormaps
custom_colormap = [white_colormap(1:threshold_index, :); black_colormap];

% Plot the p-values colormap with the custom colormap
figure(2)
imagesc(p);
colormap(custom_colormap);  % Use the custom colormap directly
colorbar;
caxis([0, 0.5]); 
title('Custom Colormap for P-Values Heatmap');
xticks(1:13);
yticks(1:13);
caxis([0, 0.5]); 
set(gca, 'TickDir', 'out')
xticklabels(variable_names);
yticklabels(variable_names);
xlabel('Variables');
ylabel('Variables');
%%
% % Create a custom colormap with bright red (for values < 0.05) and a gradient from red to white (for values >= 0.05)
% red_colormap = [linspace(1, 1, nSteps)', linspace(0, 0, nSteps)', linspace(0, 0, nSteps)'];  % Bright Red
% 
% % Determine the index where the threshold falls
% threshold_index = round(threshold * nSteps);
% 
% % Create a gradient from red to white for values >= 0.05
% gradient_colormap = [linspace(1, 1, nSteps-threshold_index)', linspace(0, 1, nSteps-threshold_index)', linspace(0, 1, nSteps-threshold_index)'];  % Gradient from Red to White
% 
% % Combine the two colormaps
% custom_colormap = [red_colormap(1:threshold_index, :); gradient_colormap];

% Plot the p-values colormap with the custom colormap
figure(2)
imagesc(p);
colormap(custom_colormap);  % Use the custom colormap directly
colorbar;
caxis([0, 0.3]); 
title('Custom Colormap for P-Values Heatmap');
xticks(1:13);
yticks(1:13);
caxis([0, 0.5]); 
 set(gca, 'TickDir', 'out')
xticklabels(variable_names);
yticklabels(variable_names);
xlabel('Variables');
ylabel('Variables');


%%
figure(2)
imagesc(p);
colormap(custom_colormap);  % You can choose a different colormap if needed
colorbar;
title('Spearman Rank p-value Heatmap');
axis square;
xticks(1:13);
yticks(1:13);
caxis([0, 0.5]); 
xticklabels(variable_names);
yticklabels(variable_names);
xlabel('Variables');
ylabel('Variables');

%%
% Extract the data for the variables of interest
relativeloc = data_cell_table.relativeloc;
nrAPsTS2 = data_cell_table.nrAPsTS2;
Rinput = data_cell_table.Rinput;
sagratio = data_cell_table.sagratio;
halfwidth = data_cell_table.Halfwidth ; 
downstroke = data_cell_table.Downstroke ;
upstroke = data_cell_table.Upstroke ;
amplitude = data_cell_table.Amplitude ;
[h, p, jbstat, cval] = jbtest(relativeloc);
if h
    disp('The data is not normally distributed.');
else
    disp('The data is normally distributed.');
end

% Calculate Spearman correlation coefficients and p-values
[rho_nrAPsTS2, p_nrAPsTS2] = corr(relativeloc, nrAPsTS2, 'Type', 'Spearman');
[rho_Rinput, p_Rinput] = corr(relativeloc, Rinput, 'Type', 'Spearman');
[rho_sagratio, p_sagratio] = corr(relativeloc, sagratio, 'Type', 'Spearman');
[rho_downstroke, p_downstroke] = corr(halfwidth, downstroke, 'Type', 'Spearman');
[rho_amplitude, p_amplitude] = corr(upstroke, amplitude, 'Type', 'Spearman');


% Perform linear regression for each correlation
fit_nrAPsTS2 = polyfit(relativeloc, nrAPsTS2, 1);
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

subplot(2, 2, 1); % Create subplot for Rinput
scatter(relativeloc, Rinput, 'filled');
hold on;
plot(relativeloc, polyval(fit_Rinput, relativeloc), 'r', 'LineWidth', 2);
xlabel('Relative Loc');
ylabel('Rinput');
title(['Spearman Correlation = ' num2str(rho_Rinput), ', p = ', num2str(p_Rinput)]);
legend('Data', 'Regression Line');

subplot(2, 2, 2); % Create subplot for sagratio
scatter(relativeloc, sagratio, 'filled');
hold on;
plot(relativeloc, polyval(fit_sagratio, relativeloc), 'r', 'LineWidth', 2);
xlabel('Relative Loc');
ylabel('Sag Ratio');
title(['Spearman Correlation = ' num2str(rho_sagratio), ', p = ', num2str(p_sagratio)]);
legend('Data', 'Regression Line');

subplot(2, 2, 3); % Create subplot for halfwidth vs downstroke
scatter(halfwidth, downstroke, 'filled');
hold on;
plot(halfwidth, polyval(fit_downstroke, halfwidth), 'r', 'LineWidth', 2);
xlabel('halfwidth');
ylabel('downstroke');
title(['Spearman Correlation = ' num2str(rho_downstroke), ', p = ', num2str(p_downstroke)]);
legend('Data', 'Regression Line');

subplot(2, 2, 4); % Create subplot for amplitude
scatter(upstroke, amplitude, 'filled');
hold on;
plot(upstroke, polyval(fit_amplitude, upstroke), 'r', 'LineWidth', 2);
xlabel('upstroke');
ylabel('amplitude');
title(['Spearman Correlation = ' num2str(rho_amplitude), ', p = ', num2str(p_amplitude)]);
legend('Data', 'Regression Line');

% Check if the correlations are significant
alpha = 0.05; % Significance level

if p_nrAPsTS2 < alpha
    disp('nrAPsTS2 is significantly correlated with Relative Loc.');
else
    disp('nrAPsTS2 is not significantly correlated with Relative Loc.');
end

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

%%
% Save the table to the specified directory
output_directory = '/Users/elinemertens/Data/Projects/Ch3.Method section/Figures/Fig4/';
output_filename = 'data_cell_table.csv';
output_path = fullfile(output_directory, output_filename);
writetable(data_table, output_path);


















