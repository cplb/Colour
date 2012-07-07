function [V_RGB, n_clip] = convert_colour_space(C_Lab, H_Lab, L_Lab, white, rgb_spec)
%[V_RGB, n_clip] = convert_colour_space(C_Lab, H_Lab, L_Lab, white, rgb_spec)
%C_Lab - Croma, H_Lab - Hue, L_Lab - Lightness. white - reference white
%point, rgb_spec - system for RGB conversion. H_Lab is given in radians,
%L_Lab should be between 0 (black) and 100 (white). C_Lab can be set to 0 
%for a greyscale. V_RGB is a vector containing red, green and blue values 
%between 0 and 1. n_clip is a list of levels which have been clipped: the
%RGB values have been readjusted to lie within the range [0 1] because the
%input LCH coordinates are outside of the permitted gamut. Suggested 
%default values are white = 3 (standard illuminat C) and rgb_spec = 1 
%(NTSC).
    
    % Define white point
	% ASTM E 308-2006
    if white == 1 % Illuminant A
        X_n = 1.0985;
        Y_n = 1.0000;
        Z_n = 0.3558;
        x_n = 0.4476;
        y_n = 0.4074;
    elseif white == 2 % B
        X_n = 0.9909;
        Y_n = 1.0000;
        Z_n = 0.8531;
        x_n = 0.3484;
        y_n = 0.3516;
    elseif white == 3 % C
        X_n = 0.98074;
        Y_n = 1.00000;
        Z_n = 1.18232;
        x_n = 0.3101;
        y_n = 0.3162;
    elseif white == 4 % D50
        X_n = 0.9641;
        Y_n = 1.0000;
        Z_n = 0.8250;
        x_n = 0.3457;
        y_n = 0.3585;
    elseif white == 5 % D55
        X_n = 0.9568;
        Y_n = 1.0000;
        Z_n = 0.9214;
        x_n = 0.3324;
        y_n = 0.3474;
    elseif white == 6 % D65
        X_n = 0.9504;
        Y_n = 1.0000;
        Z_n = 1.0888;
        x_n = 0.3127;
        y_n = 0.3290;
    elseif white == 7 % D75
        X_n = 0.9497;
        Y_n = 1.0000;
        Z_n = 1.2257;
        x_n = 0.2991;
        y_n = 0.3149;
    elseif white == 8 % ID50
        X_n = 0.9528;
        Y_n = 1.0000;
        Z_n = 0.8233;
        x_n = 0.3433;
        y_n = 0.3602;
    elseif white == 9 % ID65
        X_n = 0.9395;
        Y_n = 1.0000;
        Z_n = 1.0846;
        x_n = 0.3107;
        y_n = 0.3307;
    end

    % Convert LCH(ab) to Lab
    a_Lab = C_Lab .* cos(H_Lab);
    b_Lab = C_Lab .* sin(H_Lab);
    
    % Convert from Lab to XYZ
    % Define conversion function
    delta = 6/29;
    %f_t = (@(t)((t > delta^3) .* (t.^(1 / 3)) + (t < delta^3) .* (t / delta^2 + 2 * delta) / 3)); % Forward conversion function
    inf_t = (@(t)((t > delta^3) .* (t.^3) + (t < delta^3) .* (3 * delta^2 * (t - 2 * delta / 3)))); % Inverse conversion function
    
    Y_XYZ = Y_n * inf_t((L_Lab + 16) / 116 + zeros(size(a_Lab)));
    X_XYZ = X_n * inf_t((L_Lab + 16) / 116 + a_Lab / 500);
    Z_XYZ = Z_n * inf_t((L_Lab + 16) / 116 - b_Lab / 200);                                      
    
    % Convert from XYZ to rgb
    if rgb_spec == 1 % NTSC        
        x_rgb = [0.6700, 0.2100, 0.1400];
        y_rgb = [0.3300, 0.7100, 0.0800];
    elseif rgb_spec == 2 % sRGB
        x_rgb = [0.6400, 0.3000, 0.1500];
        y_rgb = [0.3300, 0.6000, 0.0600];
    end
    
    % Select transfromation method
    method = 1;
    if method == 0 % Inverse
        X_rgb = x_rgb ./ y_rgb;
        Y_rgb = [1, 1, 1];
        Z_rgb = (1 - x_rgb - y_rgb) ./ y_rgb;

        S_rgb = [X_rgb; Y_rgb; Z_rgb] \ [X_n; Y_n; Z_n];

        M_rgb = [S_rgb(1) * X_rgb(1), S_rgb(2) * X_rgb(2), S_rgb(3) * X_rgb(3)
                S_rgb(1) * Y_rgb(1), S_rgb(2) * Y_rgb(2), S_rgb(3) * Y_rgb(3)
                S_rgb(1) * Z_rgb(1), S_rgb(2) * Z_rgb(2), S_rgb(3) * Z_rgb(3)];

        v_rgb = M_rgb \ [X_XYZ; Y_XYZ; Z_XYZ];
    else % Forward
       % See Table 11.1 of 'Measuring Colour', Hunt & Pointer, 2011
       a_mat = [x_rgb', y_rgb', (1 - x_rgb - y_rgb)'];
       
       b_mat = det(a_mat)*inv(a_mat);
       
       %Consistency check: A_mag should be constant
       %A_mag = sum(b_mat,2);
       
       J_w = [x_n / y_n; 1; (1 - x_n - y_n) / y_n];
       
       k_w = 1 ./ (b_mat' * J_w);
       
       M_XYZ = [b_mat(:,1) * k_w(1), b_mat(:,2) * k_w(2), b_mat(:,3) * k_w(3)]';
       
       v_rgb = M_XYZ * [X_XYZ; Y_XYZ; Z_XYZ];
    end
    
    %Y_check = [0.298839, 0.586811, 0.114350] * v_rgb;
    
    % Remap into permitted range
    v_clip = max(min(v_rgb, 1), 0);
    
    % Find levels which have been clipped
    n_clip = unique([find((v_rgb(1,:) > 1)), find((v_rgb(1,:) < 0)), find((v_rgb(2,:) > 1)), find((v_rgb(2,:) < 0)), find((v_rgb(3,:) > 1)), find((v_rgb(3,:) < 0))]);
    
    % Convert from rgb to RGB
    if rgb_spec == 1 % NTSC
        gamma = 2.2; 
        V_RGB = (v_clip.^(1/gamma))';
    elseif rgb_spec == 2 % sRGB
        c_break = 3.130668442500634e-03;
        a_offset = 0.055;
        gamma = 2.4; % Effective gamma for entire range is ~2.2
        V_RGB = (((1 + a_offset) * v_clip.^(1/gamma) - a_offset) .* (v_clip > c_break) + (12.92 * v_clip) .* (v_clip <= c_break))';
    end
    
end