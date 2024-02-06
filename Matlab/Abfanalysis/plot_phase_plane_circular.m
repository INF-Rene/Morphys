function plot_smoothed_phase_plane(ApsTimeArray, ApsDataArray, DerivativesArray)
    % Plot the smoothed phase plane of action potentials.

    % Apply a moving average filter to smooth the data
    smoothed_ApsDataArray = smoothdata(ApsDataArray, 'movmean', 5);
    smoothed_ApsDerivativesArray = smoothdata(DerivativesArray, 'movmean', 5);

    % Plot the original and smoothed data
    figure;
    subplot(2,1,1);
    plot(ApsDataArray, DerivativesArray, 'b', 'LineWidth', 2);
    title('Original Phase Plane');
    xlabel('Voltage');
    ylabel('dV/dt');

    subplot(2,1,2);
    plot(smoothed_ApsDataArray, smoothed_ApsDerivativesArray, 'r', 'LineWidth', 2);
    title('Smoothed Phase Plane');
    xlabel('Voltage');
    ylabel('dV/dt');
end
