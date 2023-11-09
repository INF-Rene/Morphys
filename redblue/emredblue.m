function c = emredblue(m)
    % Adjusted Shades of Red and Blue Colormap with Extended White Area
    %   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
    %   The colors begin with bright blue, range through shades of
    %   blue to white (extended), and then through shades of red to bright red.
    %   REDBLUE, by itself, is the same length as the current figure's
    %   colormap. If no figure exists, MATLAB creates one.
    %
    %   For example, to reset the colormap of the current figure:
    %
    %             colormap(redblue)
    
    if nargin < 1
        m = size(get(gcf, 'colormap'), 1);
    end

    if (mod(m, 2) == 0)
        % Extend the white area around zero
        whiteExtent = 0.25;  % Adjust the extent as needed
        
        m1 = m * 0.5;
        
        % Initialize arrays
        r = zeros(m1, 1);
        g = zeros(m1, 1);
        b = zeros(m1, 1);
        
        % White area
        whiteLength = round(whiteExtent * m1);
        r(1:whiteLength) = 1;
        g(1:whiteLength) = 1;
        b(1:whiteLength) = 1;
        
        % Transition from white to red
        redLength = m1 - whiteLength;
        r(whiteLength + 1:end) = (0:redLength-1) / redLength;
        g(whiteLength + 1:end) = 0;
        b(whiteLength + 1:end) = 0;
    else
        % For odd 'm', keep the original logic
        m1 = floor(m * 0.5);
        r = (0:m1-1) / max(m1, 1);
        g = r;
        r = [r; ones(m1+1, 1)];
        g = [g; 1; flipud(g)];
        b = flipud(r);
    end

    c = [r g b];
end
