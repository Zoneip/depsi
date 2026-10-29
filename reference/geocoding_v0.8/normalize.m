function [line_normal,pixel_normal] = normalize(line_original, ...
                                                pixel_original, ...
                                                line_total,pixel_total)

% $$$ NORMALIZES PARSED VALUES TO [-2,2] range
    
% normalize lines
line_normal = (line_original - 0.5  * (1 + line_total));
line_normal = line_normal / (0.25 * (line_total - 1));

% normalize pixels
pixel_normal = (pixel_original - 0.5  * (1 + pixel_total));
pixel_normal = pixel_normal / (0.25 * (pixel_total - 1));
