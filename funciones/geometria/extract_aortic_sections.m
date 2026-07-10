function sections = extract_aortic_sections(data, geometry)
%EXTRACT_AORTIC_SECTIONS Cortes oblicuos por interseccion malla-Laplace.
%   Cada corte se obtiene como la interseccion entre la malla tetraedrica y el
%   nivel de Laplace seleccionado (qmeshcut). El area transversal usada para
%   el flujo se toma de Area.mat, promediando los nodos de superficie externa
%   asociados a los tetraedros intersectados. El area geometrica de qmeshcut se
%   conserva como control metodologico y fallback. Devuelve la estructura
%   'sections' con el arreglo de cortes y los vectores de area por corte.

    elem = data.elem;
    nodes = data.nodes;
    Laplace = data.Laplace;
    Area = data.Area;
    surface_node_ids = data.surface_node_ids;
    number_of_time_phases = data.number_of_time_phases;

    number_of_aortic_sections = geometry.number_of_aortic_sections;
    section_centers_mm = geometry.section_centers_mm;
    section_normal_vectors = geometry.section_normal_vectors;
    section_arc_positions_mm = geometry.section_arc_positions_mm;
    laplace_cut_levels = geometry.laplace_cut_levels;

    empty_aortic_section = struct( ...
        'center_mm', [], ...
        'normal_vector', [], ...
        'arc_position_mm', [], ...
        'laplace_level', [], ...
        'vertices_mm', [], ...
        'faces', [], ...
        'element_ids', [], ...
        'face_area_mm2', [], ...
        'geometric_area_m2', [], ...
        'area_from_area_mat_m2', [], ...
        'cross_section_area_m2', [], ...
        'area_source', [], ...
        'cut_surface_node_ids', [], ...
        'normal_velocity_mps', []);

    aortic_sections = repmat(empty_aortic_section, number_of_aortic_sections, 1);
    cross_section_area_m2 = nan(number_of_aortic_sections, 1);
    geometric_area_m2 = nan(number_of_aortic_sections, 1);
    area_source = strings(number_of_aortic_sections, 1);

    for section_index = 1:number_of_aortic_sections
        [section_vertices_mm, ~, section_faces, intersected_element_ids] = ...
            qmeshcut(elem, nodes, Laplace, laplace_cut_levels(section_index));

        if isempty(section_vertices_mm) || isempty(section_faces)
            warning('El corte %d no intersecta la malla.', section_index)
            continue
        end

        [section_vertices_mm, section_faces] = removedupnodes( ...
            section_vertices_mm, section_faces, 1e-10);

        section_face_area_mm2 = calculate_cross_section_face_areas( ...
            section_vertices_mm, section_faces);
        section_geometric_area_m2 = ...
            sum(section_face_area_mm2, 'omitnan') * 1e-6;

        element_node_ids = unique(elem(intersected_element_ids,:));
        cut_surface_node_ids = intersect(element_node_ids, surface_node_ids);

        if isempty(cut_surface_node_ids)
            warning(['Corte %d sin nodos de superficie asociados. ', ...
                'Se usa area geometrica directa como fallback.'], section_index)

            section_area_from_area_mat_m2 = NaN;
            section_area_m2 = section_geometric_area_m2;
            section_area_source = "geometric_fallback";
        else
            section_area_cm2 = mean(Area(cut_surface_node_ids), 'omitnan');
            section_area_from_area_mat_m2 = section_area_cm2 * 1e-4;
            section_area_m2 = section_area_from_area_mat_m2;
            section_area_source = "Area.mat";
        end

        if ~isfinite(section_area_m2) || section_area_m2 <= 0
            warning(['Area.mat entrego un area invalida para el corte %d. ', ...
                'Se usa area geometrica directa como fallback.'], section_index)

            section_area_m2 = section_geometric_area_m2;
            section_area_source = "geometric_fallback_invalid_AreaMat";
        end

        aortic_sections(section_index).center_mm = ...
            section_centers_mm(section_index,:);
        aortic_sections(section_index).normal_vector = ...
            section_normal_vectors(section_index,:);
        aortic_sections(section_index).arc_position_mm = ...
            section_arc_positions_mm(section_index);
        aortic_sections(section_index).laplace_level = ...
            laplace_cut_levels(section_index);
        aortic_sections(section_index).vertices_mm = section_vertices_mm;
        aortic_sections(section_index).faces = section_faces;
        aortic_sections(section_index).element_ids = intersected_element_ids;
        aortic_sections(section_index).face_area_mm2 = section_face_area_mm2;
        aortic_sections(section_index).geometric_area_m2 = ...
            section_geometric_area_m2;
        aortic_sections(section_index).area_from_area_mat_m2 = ...
            section_area_from_area_mat_m2;
        aortic_sections(section_index).cross_section_area_m2 = section_area_m2;
        aortic_sections(section_index).area_source = section_area_source;
        aortic_sections(section_index).cut_surface_node_ids = cut_surface_node_ids;
        aortic_sections(section_index).normal_velocity_mps = ...
            nan(size(section_vertices_mm,1), number_of_time_phases);

        cross_section_area_m2(section_index) = section_area_m2;
        geometric_area_m2(section_index) = section_geometric_area_m2;
        area_source(section_index) = section_area_source;
    end

    sections = struct();
    sections.list = aortic_sections;
    sections.cross_section_area_m2 = cross_section_area_m2;
    sections.geometric_area_m2 = geometric_area_m2;
    sections.area_source = area_source;
end
