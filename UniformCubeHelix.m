function colour_map = UniformCubeHelix(n_lev,H_start,L_start,L_stop,rots,colourful,Gamma)
%colour_map = UniformCubeHelix(n_lev,H_start,L_start,L_stop,rots,colourful,Gamma)
%Create color map colour_map which is monotonically varying in brightness
%and rotates uniformly in hue. n_lev - number of levels; H_start - starting
%hue (in radians: 0 ~ red, pi/2 ~ yellow, pi ~ blue-green, 3pi/2 ~ 
%purple-blue); L_start - starting lightness (0 - black, 100 - white);
%L_stop - ending lightness; rots - number of rotations (0 for constant
%hue); colourful - scaling of chroma (0 - greyscale, 1 - maximum without
%clipping); Gamma - stretching factor for intensity. colour_map consists of
%RGB values in the range [0 1]. Suggested default values
%(256,0,0,100,1.2,1,1) for black-white (coloured); (256,0,0,100,0,0,1) 
%for black-white (greyscale); (256,0,0,64.5,0,1,1) for black-red (constant 
%hue).

    % Define levels
    i_lev = 1:n_lev;
    fract = (i_lev - 1) / (n_lev - 1);
    
    % Rotation in hue
    H_angle = H_start + (2 * pi) * (rots * fract);
    
    % Luminance with gamma scaling
    L_mag = L_start + (L_stop - L_start) * fract.^Gamma;
    
    % Parameters for chroma function
    C_start = 11.0;
    C_peak = 58.0;
    L_break = 64.5;
    L_zero = 99.5;
    
    % Non-clipping chroma
    C_amp = colourful * ((C_peak - C_start) * L_mag / L_break + C_start) .* (L_mag <= L_break) + (C_peak * (L_break - L_mag)/(L_zero - L_break) + C_peak) .* ((L_mag > L_break) .* (L_mag <= L_zero));
    
    [colour_map, ~] = convert_colour_space(C_amp, H_angle, L_mag, 3, 1);
        
end