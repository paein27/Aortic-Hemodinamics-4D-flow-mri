function [sections, flow] = project_section_velocities(data, geometry, sections)
%PROJECT_SECTION_VELOCITIES Proyeccion de velocidades y calculo de flujo.
%   Para cada fase temporal se actualizan tres interpoladores espaciales
%   dispersos (uno por componente de velocidad). La velocidad de cada vertice
%   de corte se proyecta en la normal local; el flujo del corte es la media de
%   la velocidad normal ponderada por area de cara, por el area transversal.
%   Actualiza sections.list(:).normal_velocity_mps y devuelve la estructura
%   'flow' con velocidad media/peak y flujo por corte y fase.

    nodes = data.nodes;
    velocity_nodes_mps = data.velocity_nodes_mps;
    number_of_time_phases = data.number_of_time_phases;

    number_of_aortic_sections = geometry.number_of_aortic_sections;

    aortic_sections = sections.list;
    cross_section_area_m2 = sections.cross_section_area_m2;

    velocity_interpolant_x = scatteredInterpolant( ...
        nodes(:,1), nodes(:,2), nodes(:,3), zeros(size(nodes,1),1), ...
        'linear', 'none');
    velocity_interpolant_y = scatteredInterpolant( ...
        nodes(:,1), nodes(:,2), nodes(:,3), zeros(size(nodes,1),1), ...
        'linear', 'none');
    velocity_interpolant_z = scatteredInterpolant( ...
        nodes(:,1), nodes(:,2), nodes(:,3), zeros(size(nodes,1),1), ...
        'linear', 'none');

    mean_normal_velocity_mps = ...
        nan(number_of_aortic_sections, number_of_time_phases);
    peak_normal_velocity_mps = ...
        nan(number_of_aortic_sections, number_of_time_phases);
    flow_rate_m3_per_s = nan(number_of_aortic_sections, number_of_time_phases);
    flow_rate_mL_per_s = nan(number_of_aortic_sections, number_of_time_phases);

    for time_phase_index = 1:number_of_time_phases
        velocity_interpolant_x.Values = ...
            velocity_nodes_mps(:,1,time_phase_index);
        velocity_interpolant_y.Values = ...
            velocity_nodes_mps(:,2,time_phase_index);
        velocity_interpolant_z.Values = ...
            velocity_nodes_mps(:,3,time_phase_index);

        for section_index = 1:number_of_aortic_sections
            section_vertices_mm = aortic_sections(section_index).vertices_mm;

            if isempty(section_vertices_mm) || ...
                    ~isfinite(cross_section_area_m2(section_index))
                continue
            end

            section_normal_vector = ...
                aortic_sections(section_index).normal_vector;
            section_faces = aortic_sections(section_index).faces;
            section_face_area_mm2 = ...
                aortic_sections(section_index).face_area_mm2;

            velocity_x_mps = velocity_interpolant_x(section_vertices_mm);
            velocity_y_mps = velocity_interpolant_y(section_vertices_mm);
            velocity_z_mps = velocity_interpolant_z(section_vertices_mm);

            section_normal_velocity_mps = ...
                velocity_x_mps*section_normal_vector(1) + ...
                velocity_y_mps*section_normal_vector(2) + ...
                velocity_z_mps*section_normal_vector(3);

            aortic_sections(section_index).normal_velocity_mps( ...
                :, time_phase_index) = section_normal_velocity_mps;

            face_normal_velocity_mps = calculate_face_mean_vertex_values( ...
                section_normal_velocity_mps, section_faces);

            valid_face_mask = ...
                isfinite(face_normal_velocity_mps) & ...
                isfinite(section_face_area_mm2) & ...
                section_face_area_mm2 > 0;

            if ~any(valid_face_mask)
                continue
            end

            mean_normal_velocity_mps(section_index,time_phase_index) = sum( ...
                face_normal_velocity_mps(valid_face_mask) .* ...
                section_face_area_mm2(valid_face_mask)) / ...
                sum(section_face_area_mm2(valid_face_mask));

            valid_vertex_velocity_mps = ...
                section_normal_velocity_mps(isfinite(section_normal_velocity_mps));

            if ~isempty(valid_vertex_velocity_mps)
                peak_normal_velocity_mps(section_index,time_phase_index) = ...
                    max(valid_vertex_velocity_mps);
            end

            flow_rate_m3_per_s(section_index,time_phase_index) = ...
                mean_normal_velocity_mps(section_index,time_phase_index) * ...
                cross_section_area_m2(section_index);
            flow_rate_mL_per_s(section_index,time_phase_index) = ...
                flow_rate_m3_per_s(section_index,time_phase_index) * 1e6;
        end
    end

    sections.list = aortic_sections;

    flow = struct();
    flow.mean_normal_velocity_mps = mean_normal_velocity_mps;
    flow.peak_normal_velocity_mps = peak_normal_velocity_mps;
    flow.flow_rate_m3_per_s = flow_rate_m3_per_s;
    flow.flow_rate_mL_per_s = flow_rate_mL_per_s;
end
