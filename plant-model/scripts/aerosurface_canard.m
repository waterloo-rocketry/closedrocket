function [canard_pos_x, canard_pos_r_mean, canard_area] = aerosurface_canard(canard_chord_root, canard_height, canard_sweep_angle, canard_spine_width, canard_spine_effectiveness, canard_radial_offset, canard_pos_x_roottip, rocket_diameter)

    %% ASSUME COP FOR RECTANGLE IS AT CENTRE SUPERSONIC

    area_delta = canard_chord_root * canard_height / 2;
    area_spine = canard_spine_effectiveness * canard_spine_width * canard_chord_root;
    canard_area = area_delta + area_spine; 
    
    % center of pressure radial
    height_aero_mean_delta = canard_height / 3; %mean aerodynamic chord distance, radially from root
    height_aero_mean_spine = - canard_spine_width / 2;
    height_aero_mean = (area_delta * height_aero_mean_delta + area_spine * height_aero_mean_spine) / canard_area;
    canard_pos_r_mean = rocket_diameter/2 + height_aero_mean + canard_radial_offset; % mean aerodynamic chord distance with the radius added

    % center of pressure axial
    chord_aero_mean_delta = 2/3 * canard_chord_root;
    pos_x_aero_mean_delta = height_aero_mean * tan(canard_sweep_angle);
    delta_cp_loc = - pos_x_aero_mean_delta - 1/4 * chord_aero_mean_delta;
    spine_cp_loc = - 0.5 * canard_chord_root;

    canard_cp_loc = (area_delta * delta_cp_loc + area_spine * spine_cp_loc) / canard_area;
    canard_pos_x = canard_pos_x_roottip + canard_cp_loc;
end