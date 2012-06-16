function colour_map = UniformCubeHelix(n_lev,H_start,rots,colourful,Gamma)

    %colour_map = zeros(n_lev,3);

    i_lev = 1:n_lev;
    fract = (i_lev - 1) / (n_lev - 1);
    H_angle = 2*pi*(rots * fract) + H_start;
    fract = fract.^Gamma;
    C_amp = colourful * 0.5 * fract .* (1 - fract) * 100;
    [colour_map, ~] = convert_colour_space(C_amp, H_angle, 100 * fract);
        
end