function [coordinate_scale_to_mm, detected_coordinate_units] = ...
    detect_coordinate_scale_to_mm(nodes, Centerline)

    spatial_bounds = [
        min(nodes, [], 1);
        max(nodes, [], 1);
        min(Centerline, [], 1);
        max(Centerline, [], 1)];
    spatial_extent = max(max(spatial_bounds) - min(spatial_bounds));

    if spatial_extent < 1
        coordinate_scale_to_mm = 1000;
        detected_coordinate_units = 'm';
    else
        coordinate_scale_to_mm = 1;
        detected_coordinate_units = 'mm';
    end
end
