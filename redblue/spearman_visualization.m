% Define the number of colormap steps (you can adjust this value)
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
