function geometry = build_section_geometry(data, config)
%BUILD_SECTION_GEOMETRY Geometria de los cortes transversales.
%   Distribuye N cortes equidistantes sobre la centerline, calcula el centro
%   y la normal (tangente de la centerline) de cada uno, ancla las posiciones
%   anatomicas I y IV (Urbina/Sotelo) por longitud de arco fija desde la raiz
%   y asigna a cada corte su nivel de Laplace. Devuelve la estructura
%   'geometry'.

    nodes = data.nodes;
    Centerline = data.Centerline;
    Laplace = data.Laplace;

    number_of_aortic_sections = config.number_of_aortic_sections;
    position_I_arc_from_root_mm = config.position_I_arc_from_root_mm;
    position_IV_arc_from_root_mm = config.position_IV_arc_from_root_mm;

    % ---- Posiciones equidistantes sobre la centerline ----
    centerline_segment_lengths_mm = sqrt(sum(diff(Centerline, 1, 1).^2, 2));
    unique_centerline_point_mask = [true; centerline_segment_lengths_mm > eps];
    centerline_points_mm = Centerline(unique_centerline_point_mask, :);

    centerline_segment_lengths_mm = ...
        sqrt(sum(diff(centerline_points_mm, 1, 1).^2, 2));
    centerline_arc_length_mm = [0; cumsum(centerline_segment_lengths_mm)];
    centerline_total_length_mm = centerline_arc_length_mm(end);

    section_arc_positions_mm = ...
        linspace(0, centerline_total_length_mm, number_of_aortic_sections + 2);
    section_arc_positions_mm = section_arc_positions_mm(2:end-1);
    section_arc_positions_m = section_arc_positions_mm / 1000;

    % Se eligen los planos I y IV como los cortes cuyo arco desde la raiz queda
    % mas cerca de las posiciones fijas configuradas. Asi I y IV caen en el
    % mismo lugar anatomico en todos los datasets (dentro de +/- media
    % separacion) y la distancia arc(IV)-arc(I) deja de depender de cuanto se
    % segmento la aorta. La distancia correcta es esa resta (la funcion del par
    % ya la calcula); NO se mueve el origen al inflow, porque el delay se
    % obtiene correlacionando las ondas de I y IV y la distancia debe empezar
    % en I.
    [~, section_I_urbina] = ...
        min(abs(section_arc_positions_mm - position_I_arc_from_root_mm));
    [~, section_IV_urbina] = ...
        min(abs(section_arc_positions_mm - position_IV_arc_from_root_mm));

    urbina_section_arc_positions_m = section_arc_positions_m(:);

    section_centers_mm = interp1( ...
        centerline_arc_length_mm, ...
        centerline_points_mm, ...
        section_arc_positions_mm, ...
        'linear');

    centerline_tangent_x = gradient( ...
        centerline_points_mm(:,1), centerline_arc_length_mm);
    centerline_tangent_y = gradient( ...
        centerline_points_mm(:,2), centerline_arc_length_mm);
    centerline_tangent_z = gradient( ...
        centerline_points_mm(:,3), centerline_arc_length_mm);

    centerline_tangent_vectors = ...
        [centerline_tangent_x, centerline_tangent_y, centerline_tangent_z];

    section_normal_vectors = interp1( ...
        centerline_arc_length_mm, ...
        centerline_tangent_vectors, ...
        section_arc_positions_mm, ...
        'linear');
    section_normal_vectors = ...
        section_normal_vectors ./ vecnorm(section_normal_vectors, 2, 2);

    fprintf('Numero de cortes: %d\n', number_of_aortic_sections)
    fprintf('Longitud de la centerline: %.1f mm\n', centerline_total_length_mm)
    fprintf('Separacion entre cortes: %.1f mm\n', ...
        mean(diff(section_arc_positions_mm)))
    fprintf(['Anclaje I-IV por arco fijo: I=%.1f mm (plano %d), ', ...
        'IV=%.1f mm (plano %d), distancia I-IV=%.1f mm\n'], ...
        section_arc_positions_mm(section_I_urbina), section_I_urbina, ...
        section_arc_positions_mm(section_IV_urbina), section_IV_urbina, ...
        section_arc_positions_mm(section_IV_urbina) - ...
        section_arc_positions_mm(section_I_urbina))

    % ---- Niveles de Laplace asociados a cada corte ----
    nearest_node_to_section_center = zeros(number_of_aortic_sections, 1);
    laplace_cut_levels = zeros(number_of_aortic_sections, 1);

    for section_index = 1:number_of_aortic_sections
        distance_from_nodes_to_center = sqrt(sum( ...
            (nodes - section_centers_mm(section_index,:)).^2, 2));

        [~, nearest_node_to_section_center(section_index)] = ...
            min(distance_from_nodes_to_center);

        laplace_cut_levels(section_index) = ...
            Laplace(nearest_node_to_section_center(section_index));
    end

    if numel(unique(laplace_cut_levels)) < number_of_aortic_sections
        warning('Dos o mas cortes usan el mismo nivel de Laplace.')
    end

    geometry = struct();
    geometry.number_of_aortic_sections = number_of_aortic_sections;
    geometry.centerline_points_mm = centerline_points_mm;
    geometry.centerline_arc_length_mm = centerline_arc_length_mm;
    geometry.centerline_total_length_mm = centerline_total_length_mm;
    geometry.section_arc_positions_mm = section_arc_positions_mm;
    geometry.section_arc_positions_m = section_arc_positions_m;
    geometry.urbina_section_arc_positions_m = urbina_section_arc_positions_m;
    geometry.section_I_urbina = section_I_urbina;
    geometry.section_IV_urbina = section_IV_urbina;
    geometry.section_centers_mm = section_centers_mm;
    geometry.section_normal_vectors = section_normal_vectors;
    geometry.laplace_cut_levels = laplace_cut_levels;
end
