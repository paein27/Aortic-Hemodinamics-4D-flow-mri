clc
close all
clear

%% Carga de datos procesados
normal_rest_data_folder = 'C:\Users\ignac\OneDrive\Desktop\USM\Memoria\Data\Datos_Procesados\13 mm rest\MATLAB FILES';

addpath(fullfile('C:\Users\ignac\OneDrive\Desktop\USM\Memoria\4D-Flow-Matlab-Toolbox-main\iso2mesh'))

load(fullfile(normal_rest_data_folder, 'FE Mesh', 'elem.mat'))
load(fullfile(normal_rest_data_folder, 'FE Mesh', 'nodes.mat'))
load(fullfile(normal_rest_data_folder, 'FE Mesh', 'faces.mat'))
load(fullfile(normal_rest_data_folder, 'FE Laplace', 'Laplace.mat'))
load(fullfile(normal_rest_data_folder, 'FE Velocity', 'VEL.mat'))
load(fullfile(normal_rest_data_folder, 'FE Centerline', 'Centerline.mat'))
load(fullfile(normal_rest_data_folder, 'FE Area', 'Area.mat'))

[coordinate_scale_to_mm, detected_coordinate_units] = ...
    detect_coordinate_scale_to_mm(nodes, Centerline);
nodes = nodes * coordinate_scale_to_mm;
Centerline = Centerline * coordinate_scale_to_mm;

fprintf('Unidad de coordenadas detectada: %s\n', detected_coordinate_units)

surface_node_ids = unique(faces(:));

time_file = fullfile(normal_rest_data_folder, '2D Flow', 'Time.mat');

%% Configuracion general
number_of_aortic_sections = 30;
section_I_urbina = 3;
section_IV_urbina = 30;
number_of_time_phases = size(VEL, 3);

[time_vector_seconds, time_vector_source] = get_time_vector_for_pwv(time_file, number_of_time_phases);

fprintf('Numero de fases temporales: %d\n', number_of_time_phases)

%% Deteccion y conversion de unidades de velocidad
velocity_magnitude_samples = abs(VEL(:));
velocity_magnitude_samples = ...
    velocity_magnitude_samples(isfinite(velocity_magnitude_samples));
velocity_magnitude_samples = sort(velocity_magnitude_samples);

percentile_99_index = max(1, round(0.99*numel(velocity_magnitude_samples)));
velocity_percentile_99 = velocity_magnitude_samples(percentile_99_index);

if velocity_percentile_99 > 20
    velocity_scale_to_mps = 0.01;
    detected_velocity_units = 'cm/s';
else
    velocity_scale_to_mps = 1.0;
    detected_velocity_units = 'm/s';
end

velocity_nodes_mps = VEL * velocity_scale_to_mps;

fprintf('Unidad de VEL detectada: %s\n', detected_velocity_units)

%--------------------- Posiciones equidistantes sobre la centerline -----------------------%

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

%----------------------- Niveles de Laplace asociados a cada corte -----------------------%

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

%----------------------- Generacion de cortes oblicuos con qmeshcut -----------------------%

% Cada corte se obtiene como la interseccion entre la malla tetraedrica y el
% nivel de Laplace seleccionado. El area transversal usada para el flujo se
% obtiene desde Area.mat, promediando los nodos de superficie externa
% asociados a los tetraedros intersectados. El area geometrica de qmeshcut se
% conserva como control metodologico y fallback.

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

%% Interpolacion y proyeccion de velocidades sobre cada corte
% Para cada fase temporal se actualizan tres interpoladores espaciales
% dispersos, uno por componente de velocidad. La velocidad de cada vertice de
% corte se proyecta en la normal local

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

%% Visualizacion 3D de malla, centerline y cortes

average_flow_per_time_phase_mL_per_s = ...
    mean(flow_rate_mL_per_s, 1, 'omitnan');

if all(~isfinite(average_flow_per_time_phase_mL_per_s))
    error('No fue posible calcular flujo valido en ningun corte.')
elseif any(average_flow_per_time_phase_mL_per_s > 0)
    [~, visualization_phase_index] = ...
        max(average_flow_per_time_phase_mL_per_s);
else
    [~, visualization_phase_index] = ...
        max(abs(average_flow_per_time_phase_mL_per_s));
    warning(['No se encontro una fase con flujo medio positivo. ', ...
        'Revise la orientacion de la centerline.'])
end

figure('Name', 'Normal rest: malla, 30 cortes y velocidad normal')
patch( ...
    'Faces', faces, ...
    'Vertices', nodes, ...
    'EdgeColor', 'none', ...
    'FaceColor', [0.75 0.75 0.75], ...
    'FaceAlpha', 0.15)
hold on
plot3( ...
    centerline_points_mm(:,1), ...
    centerline_points_mm(:,2), ...
    centerline_points_mm(:,3), ...
    'r', 'LineWidth', 1.5)

velocity_vector_scale_mm_per_mps = 8;

for section_index = 1:number_of_aortic_sections
    section_vertices_mm = aortic_sections(section_index).vertices_mm;
    section_faces = aortic_sections(section_index).faces;

    if isempty(section_vertices_mm)
        continue
    end

    section_normal_velocity_mps = aortic_sections(section_index). ...
        normal_velocity_mps(:, visualization_phase_index);
    section_normal_vector = aortic_sections(section_index).normal_vector;

    patch( ...
        'Faces', section_faces, ...
        'Vertices', section_vertices_mm, ...
        'FaceVertexCData', section_normal_velocity_mps, ...
        'FaceColor', 'interp', ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.85)

    quiver_step = max(1, ceil(size(section_vertices_mm,1)/100));
    quiver_vertex_indices = 1:quiver_step:size(section_vertices_mm,1);
    quiver_velocity_values_mps = ...
        section_normal_velocity_mps(quiver_vertex_indices);
    valid_quiver_mask = isfinite(quiver_velocity_values_mps);

    quiver_vertex_indices = quiver_vertex_indices(valid_quiver_mask);
    quiver_velocity_values_mps = ...
        quiver_velocity_values_mps(valid_quiver_mask);

    quiver3( ...
        section_vertices_mm(quiver_vertex_indices,1), ...
        section_vertices_mm(quiver_vertex_indices,2), ...
        section_vertices_mm(quiver_vertex_indices,3), ...
        section_normal_vector(1) * ...
            quiver_velocity_values_mps * velocity_vector_scale_mm_per_mps, ...
        section_normal_vector(2) * ...
            quiver_velocity_values_mps * velocity_vector_scale_mm_per_mps, ...
        section_normal_vector(3) * ...
            quiver_velocity_values_mps * velocity_vector_scale_mm_per_mps, ...
        0, 'k', 'LineWidth', 0.6)
end

axis equal
grid on
view(3)
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Z [mm]')
title(sprintf( ...
    '30 cortes qmeshcut y velocidad normal | fase %d', ...
    visualization_phase_index))
colorbar
colormap turbo
hold off

%% Resumen de area y flujo por corte

peak_flow_per_section_mL_per_s = nan(number_of_aortic_sections, 1);
mean_flow_per_section_mL_per_s = nan(number_of_aortic_sections, 1);

for section_index = 1:number_of_aortic_sections
    section_flow_curve_mL_per_s = flow_rate_mL_per_s(section_index,:);
    section_flow_curve_mL_per_s = ...
        section_flow_curve_mL_per_s(isfinite(section_flow_curve_mL_per_s));

    if ~isempty(section_flow_curve_mL_per_s)
        peak_flow_per_section_mL_per_s(section_index) = ...
            max(section_flow_curve_mL_per_s);
        mean_flow_per_section_mL_per_s(section_index) = ...
            mean(section_flow_curve_mL_per_s);
    end
end

FlowSummary = table( ...
    (1:number_of_aortic_sections)', ...
    section_arc_positions_m(:), ...
    laplace_cut_levels(:), ...
    cross_section_area_m2, ...
    cross_section_area_m2*1e6, ...
    geometric_area_m2, ...
    geometric_area_m2*1e6, ...
    cross_section_area_m2 ./ geometric_area_m2, ...
    area_source, ...
    peak_flow_per_section_mL_per_s, ...
    mean_flow_per_section_mL_per_s, ...
    'VariableNames', { ...
        'Plane', ...
        's_m', ...
        'LaplaceLevel', ...
        'AreaMatUsed_m2', ...
        'AreaMatUsed_mm2', ...
        'GeometricArea_m2', ...
        'GeometricArea_mm2', ...
        'AreaMat_to_Geometric_Ratio', ...
        'AreaSource', ...
        'PeakFlow_mLs', ...
        'MeanFlow_mLs'});

disp(FlowSummary)

%% Graficar curvas de flujo por corte

figure('Name', 'Normal rest: curvas de flujo por corte')
plot(time_vector_seconds, flow_rate_mL_per_s', 'LineWidth', 1.2)
xlabel('Tiempo [s]')
ylabel('Flujo [mL/s]')
title('Flujo = velocidad normal media por area transversal')
legend(compose('Corte %d', 1:number_of_aortic_sections), ...
    'Location', 'eastoutside')
grid on

%% Calculo de PWV con TTP, TTF y correlacion cruzada

%   TTP  = Time To Peak, usa el tiempo del maximo de flujo.
%   TTF  = Time To Foot, estima el pie de onda desde la pendiente inicial.
%   XCor = Cross-correlation global tipo Markl.
%
% La distancia usada en el ajuste lineal es la posicion acumulada sobre la
% centerline, en metros. PWV se calcula como 1/pendiente del ajuste tiempo
% vs distancia.

[PWV_TTP_mps, TTP_s, fit_TTP, R2_TTP] = estimate_pwv_ttp( ...
    section_arc_positions_m(:), time_vector_seconds, flow_rate_mL_per_s);

[PWV_TTF_mps, TTF_s, fit_TTF, R2_TTF] = estimate_pwv_ttf( ...
    section_arc_positions_m(:), time_vector_seconds, flow_rate_mL_per_s);

xcor_reference_section_index = section_I_urbina;

[PWV_XCor_mps, Delay_XCor_s, fit_XCor, R2_XCor] = ...
    estimate_pwv_xcor( ...
        section_arc_positions_m(:), ...
        time_vector_seconds, ...
        flow_rate_mL_per_s, ...
        xcor_reference_section_index);

valid_ttp_section_count = nnz(isfinite(TTP_s));
valid_ttf_section_count = nnz(isfinite(TTF_s));
valid_xcor_section_count = nnz(isfinite(Delay_XCor_s));

PWVSummary = table( ...
    ["TTP"; "TTF"; "XCor"], ...
    [PWV_TTP_mps; PWV_TTF_mps; PWV_XCor_mps], ...
    [fit_TTP(1); fit_TTF(1); fit_XCor(1)], ...
    [R2_TTP; R2_TTF; R2_XCor], ...
    [valid_ttp_section_count; ...
        valid_ttf_section_count; ...
        valid_xcor_section_count], ...
    repmat(number_of_aortic_sections, 3, 1), ...
    ["TTP all sections"; ...
        "TTF baseline-corrected all sections"; ...
        "XCor: all sections vs reference section"], ...
    'VariableNames', { ...
        'Method', ...
        'PWV_mps', ...
        'Slope_s_per_m', ...
        'R2', ...
        'ValidPlanes', ...
        'TotalPlanes', ...
        'Notes'});

PWVPlaneSummary = table( ...
    (1:number_of_aortic_sections)', ...
    section_arc_positions_m(:), ...
    TTP_s(:), ...
    TTF_s(:), ...
    Delay_XCor_s(:), ...
    'VariableNames', {'Plane', 's_m', 'TTP_s', 'TTF_s', 'XCor'});

fprintf('\n============================================================\n')
fprintf('RESULTADOS NORMAL REST\n')
fprintf('============================================================\n')
fprintf('Numero de fases temporales: %d\n', numel(time_vector_seconds))
fprintf('Tiempo inicial: %.6f s\n', time_vector_seconds(1))
fprintf('Tiempo final: %.6f s\n', time_vector_seconds(end))
fprintf('dt medio: %.6f s\n', mean(diff(time_vector_seconds)))
fprintf('\nNumero de cortes: %d\n', number_of_aortic_sections)
fprintf('XCor %d\n', ...
    xcor_reference_section_index)
fprintf(['PWV XCor global = %.4f m/s | slope = %.6f s/m | ', ...
    'R2 = %.4f | planos validos = %d\n'], ...
    PWV_XCor_mps, fit_XCor(1), R2_XCor, valid_xcor_section_count)
fprintf('Separacion media: %.1f mm\n', mean(diff(section_arc_positions_mm)))
fprintf('Fase visualizada: %d\n', visualization_phase_index)
fprintf('\nResumen de PWV:\n')
disp(PWVSummary)
fprintf('\nTiempos caracteristicos por corte:\n')
disp(PWVPlaneSummary)

plot_pwv_fit( ...
    section_arc_positions_m, TTP_s, fit_TTP, PWV_TTP_mps, R2_TTP, ...
    'PWV por TTP', 'TTP [s]');
plot_pwv_fit( ...
    section_arc_positions_m, TTF_s, fit_TTF, PWV_TTF_mps, R2_TTF, ...
    'PWV por TTF', 'TTF [s]');

relative_distance_xcor_m = ...
    section_arc_positions_m(:) - ...
    section_arc_positions_m(xcor_reference_section_index);
valid_xcor_plot_mask = ...
    isfinite(relative_distance_xcor_m) & isfinite(Delay_XCor_s);
x_fit = linspace( ...
    min(relative_distance_xcor_m(valid_xcor_plot_mask)), ...
    max(relative_distance_xcor_m(valid_xcor_plot_mask)), ...
    200);
y_fit = polyval(fit_XCor, x_fit);

figure('Name', 'PWV por XCor')
plot( ...
    relative_distance_xcor_m(valid_xcor_plot_mask), ...
    Delay_XCor_s(valid_xcor_plot_mask), ...
    'o', 'MarkerSize', 7, 'LineWidth', 1.5)
hold on
plot(x_fit, y_fit, '-', 'LineWidth', 1.8)

xlabel('Distancia sobre centerline [m]')
ylabel('Delay relativo [s]')
title(sprintf( ...
    'PWV por XCor global = %.3f m/s | R2 = %.3f', ...
    PWV_XCor_mps, R2_XCor))
legend({'Delay por corte', 'Ajuste lineal'}, 'Location', 'best')
grid on
hold off

%% Funciones auxiliares
% Las funciones locales mantienen aislados los calculos repetitivos:
%   - area de caras triangulares/cuadrilaterales del corte,
%   - promedios por cara desde valores de vertice,
%   - lectura robusta del eje temporal,
%   - estimadores de PWV usados tambien como referencia en main.m.

function face_area_mm2 = calculate_cross_section_face_areas( ...
    vertices_mm, face_node_indices)

    % qmeshcut puede devolver caras triangulares o cuadrilateras. Para
    % homogeneizar, cada cara se divide en dos triangulos: (1,2,3) y
    % (1,3,4). Si la cara es triangular, qmeshcut repite el indice 3 en la
    % columna 4, por lo que el segundo triangulo tiene area cero.

    first_vertex_mm = vertices_mm(face_node_indices(:,1),:);
    second_vertex_mm = vertices_mm(face_node_indices(:,2),:);
    third_vertex_mm = vertices_mm(face_node_indices(:,3),:);
    fourth_vertex_mm = vertices_mm(face_node_indices(:,4),:);

    first_triangle_area_mm2 = 0.5 * vecnorm( ...
        cross(second_vertex_mm - first_vertex_mm, ...
        third_vertex_mm - first_vertex_mm, 2), 2, 2);
    second_triangle_area_mm2 = 0.5 * vecnorm( ...
        cross(third_vertex_mm - first_vertex_mm, ...
        fourth_vertex_mm - first_vertex_mm, 2), 2, 2);

    face_area_mm2 = first_triangle_area_mm2 + second_triangle_area_mm2;
end


function face_mean_value = calculate_face_mean_vertex_values( ...
    vertex_values, face_node_indices)

    % Promedia los valores nodales de cada cara. Las caras triangulares se
    % tratan con solo sus tres vertices reales para evitar contar dos veces
    % el vertice repetido por qmeshcut.

    values_per_face = vertex_values(face_node_indices);
    face_mean_value = mean(values_per_face, 2);

    triangular_face_mask = face_node_indices(:,3) == face_node_indices(:,4);
    face_mean_value(triangular_face_mask) = ...
        mean(values_per_face(triangular_face_mask,1:3), 2);
end


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


function [time_vector_seconds, time_vector_source] = get_time_vector_for_pwv( ...
    time_file, number_of_time_phases)

    % Busca dentro de Time.mat la primera variable numerica compatible con
    % el numero de fases de VEL. Si el eje temporal parece estar en ms, lo
    % convierte a segundos y lo desplaza para iniciar en t = 0.

    time_vector_seconds = [];
    time_vector_source = '';

    if exist(time_file, 'file') == 2
        time_data = load(time_file);
        variable_names = fieldnames(time_data);

        for variable_index = 1:numel(variable_names)
            candidate_time_vector = time_data.(variable_names{variable_index});

            if isnumeric(candidate_time_vector) && ~isempty(candidate_time_vector)
                candidate_time_vector = squeeze(candidate_time_vector);

                if isvector(candidate_time_vector) && ...
                        numel(candidate_time_vector) >= number_of_time_phases

                    time_vector_seconds = double(candidate_time_vector(:));
                    time_vector_seconds = ...
                        time_vector_seconds(1:number_of_time_phases);
                    time_vector_source = ...
                        ['Time.mat -> variable: ', variable_names{variable_index}];
                    break
                end
            end
        end
    end

    if isempty(time_vector_seconds)
        error('No se encontro un eje temporal valido en %s.', time_file)
    end

    if max(abs(time_vector_seconds)) > 10
        time_vector_seconds = time_vector_seconds / 1000;
    end

    time_vector_seconds = time_vector_seconds - time_vector_seconds(1);

    if numel(time_vector_seconds) ~= number_of_time_phases
        error( ...
            'El eje temporal tiene %d puntos, pero VEL tiene %d fases.', ...
            numel(time_vector_seconds), number_of_time_phases)
    end

    if any(~isfinite(time_vector_seconds)) || ...
            any(diff(time_vector_seconds) <= 0)
        error('El eje temporal debe ser finito y estrictamente creciente.')
    end
end


function [PWV_mps, TTP_s, fit_coefficients, R2] = estimate_pwv_ttp( ...
    section_distance_m, time_vector_seconds, flow_rate_mL_per_s)

    % Time To Peak: para cada corte identifica el tiempo del maximo de flujo.
    % Si el maximo no esta en los extremos, ajusta una parabola local de tres
    % puntos para refinar el instante del peak con resolucion sub-fase.

    number_of_sections = size(flow_rate_mL_per_s, 1);
    TTP_s = nan(number_of_sections, 1);
    mean_time_step_seconds = mean(diff(time_vector_seconds));

    for section_index = 1:number_of_sections
        section_flow_curve_mL_per_s = flow_rate_mL_per_s(section_index,:)';

        if all(~isfinite(section_flow_curve_mL_per_s))
            continue
        end

        section_flow_curve_mL_per_s = fillmissing( ...
            section_flow_curve_mL_per_s, 'linear', 'EndValues', 'nearest');
        [~, peak_sample_index] = max(section_flow_curve_mL_per_s);
        floating_peak_sample_index = peak_sample_index - 1;

        if peak_sample_index > 1 && ...
                peak_sample_index < numel(section_flow_curve_mL_per_s)

            left_flow_value = section_flow_curve_mL_per_s(peak_sample_index-1);
            center_flow_value = section_flow_curve_mL_per_s(peak_sample_index);
            right_flow_value = section_flow_curve_mL_per_s(peak_sample_index+1);
            parabola_denominator = ...
                left_flow_value - 2*center_flow_value + right_flow_value;

            if isfinite(parabola_denominator) && ...
                    abs(parabola_denominator) > eps

                peak_sample_offset = ...
                    0.5*(left_flow_value - right_flow_value) / ...
                    parabola_denominator;
                peak_sample_offset = max(min(peak_sample_offset, 1), -1);
                floating_peak_sample_index = ...
                    peak_sample_index - 1 + peak_sample_offset;
            end
        end

        TTP_s(section_index) = time_vector_seconds(1) + ...
            floating_peak_sample_index*mean_time_step_seconds;
    end

    [PWV_mps, fit_coefficients, R2] = fit_pwv_linear( ...
        section_distance_m, TTP_s);
end


function [PWV_mps, TTF_s, fit_coefficients, R2] = estimate_pwv_ttf( ...
    section_distance_m, time_vector_seconds, flow_rate_mL_per_s)

    % Time To Foot: usa el tramo ascendente entre 20% y 80% del peak de
    % flujo. Una recta sobre ese tramo se extrapola hasta flujo cero para
    % estimar el pie de onda.

    number_of_sections = size(flow_rate_mL_per_s, 1);
    TTF_s = nan(number_of_sections, 1);

    for section_index = 1:number_of_sections
        section_flow_curve_mL_per_s = flow_rate_mL_per_s(section_index,:)';

        if all(~isfinite(section_flow_curve_mL_per_s))
            continue
        end

        section_flow_curve_mL_per_s = fillmissing( ...
            section_flow_curve_mL_per_s, 'linear', 'EndValues', 'nearest');
        baseline_sample_count = min(3, numel(section_flow_curve_mL_per_s));
        baseline_flow_mL_per_s = median( ...
            section_flow_curve_mL_per_s(1:baseline_sample_count), ...
            'omitnan');
        flow_relative_mL_per_s = ...
            section_flow_curve_mL_per_s - baseline_flow_mL_per_s;
        [peak_relative_flow_mL_per_s, peak_sample_index] = ...
            max(flow_relative_mL_per_s);

        if peak_sample_index < 3 || ...
                ~isfinite(peak_relative_flow_mL_per_s) || ...
                peak_relative_flow_mL_per_s <= 0
            continue
        end

        threshold_20_percent = 0.20 * peak_relative_flow_mL_per_s;
        threshold_80_percent = 0.80 * peak_relative_flow_mL_per_s;
        threshold_20_index = find( ...
            flow_relative_mL_per_s(1:peak_sample_index) >= ...
            threshold_20_percent, 1, 'first');
        threshold_80_index = find( ...
            flow_relative_mL_per_s(1:peak_sample_index) >= ...
            threshold_80_percent, 1, 'first');

        if isempty(threshold_20_index) || ...
                isempty(threshold_80_index) || ...
                threshold_80_index <= threshold_20_index
            continue
        end

        upslope_sample_indices = ( ...
            max(1, threshold_20_index-1): ...
            min(peak_sample_index, threshold_80_index+1))';

        if numel(upslope_sample_indices) < 2
            continue
        end

        upslope_fit_coefficients = polyfit( ...
            time_vector_seconds(upslope_sample_indices), ...
            flow_relative_mL_per_s(upslope_sample_indices), ...
            1);

        if ~isfinite(upslope_fit_coefficients(1)) || ...
                abs(upslope_fit_coefficients(1)) <= eps
            continue
        end

        TTF_s(section_index) = ...
            -upslope_fit_coefficients(2) / upslope_fit_coefficients(1);
    end

    [PWV_mps, fit_coefficients, R2] = fit_pwv_linear( ...
        section_distance_m, TTF_s);
end


function [PWV_mps, Delay_s, fit_coefficients, R2] = estimate_pwv_xcor( ...
    section_distance_m, time_vector_seconds, flow_rate_mL_per_s, ...
    reference_section_index)

    % Cross-correlation: compara cada curva de flujo contra un corte de
    % referencia. El lag de maxima correlacion se transforma en retraso
    % temporal, y luego se ajusta retraso vs distancia relativa.

    number_of_sections = size(flow_rate_mL_per_s, 1);
    number_of_time_phases = size(flow_rate_mL_per_s, 2);
    Delay_s = nan(number_of_sections, 1);
    mean_time_step_seconds = mean(diff(time_vector_seconds));

    reference_flow_curve_mL_per_s = ...
        flow_rate_mL_per_s(reference_section_index,:)';

    if all(~isfinite(reference_flow_curve_mL_per_s))
        PWV_mps = NaN;
        fit_coefficients = [NaN NaN];
        R2 = NaN;
        return
    end

    reference_flow_curve_mL_per_s = fillmissing( ...
        reference_flow_curve_mL_per_s, 'linear', 'EndValues', 'nearest');
    reference_flow_curve_mL_per_s = ...
        reference_flow_curve_mL_per_s - ...
        mean(reference_flow_curve_mL_per_s, 'omitnan');

    maximum_lag_seconds = 0.30;
    maximum_lag_samples = min( ...
        floor(number_of_time_phases/2), ...
        max(1, round(maximum_lag_seconds/mean_time_step_seconds)));
    lag_sample_vector = (-maximum_lag_samples:maximum_lag_samples)';

    for section_index = 1:number_of_sections
        section_flow_curve_mL_per_s = flow_rate_mL_per_s(section_index,:)';

        if all(~isfinite(section_flow_curve_mL_per_s))
            continue
        end

        section_flow_curve_mL_per_s = fillmissing( ...
            section_flow_curve_mL_per_s, 'linear', 'EndValues', 'nearest');
        section_flow_curve_mL_per_s = ...
            section_flow_curve_mL_per_s - ...
            mean(section_flow_curve_mL_per_s, 'omitnan');
        correlation_by_lag = nan(numel(lag_sample_vector), 1);

        for lag_index = 1:numel(lag_sample_vector)
            current_lag_samples = lag_sample_vector(lag_index);

            if current_lag_samples >= 0
                reference_segment = ...
                    reference_flow_curve_mL_per_s(1:end-current_lag_samples);
                section_segment = ...
                    section_flow_curve_mL_per_s(1+current_lag_samples:end);
            else
                reference_segment = ...
                    reference_flow_curve_mL_per_s(1-current_lag_samples:end);
                section_segment = ...
                    section_flow_curve_mL_per_s(1:end+current_lag_samples);
            end

            valid_sample_mask = ...
                isfinite(reference_segment) & isfinite(section_segment);

            if nnz(valid_sample_mask) < 3
                continue
            end

            reference_segment = reference_segment(valid_sample_mask);
            section_segment = section_segment(valid_sample_mask);
            correlation_denominator = sqrt( ...
                sum(reference_segment.^2) * sum(section_segment.^2));

            if correlation_denominator > 0
                correlation_by_lag(lag_index) = ...
                    sum(reference_segment .* section_segment) / ...
                    correlation_denominator;
            end
        end

        if all(~isfinite(correlation_by_lag))
            continue
        end

        [~, maximum_correlation_index] = max(correlation_by_lag);
        refined_lag_samples = lag_sample_vector(maximum_correlation_index);

        if maximum_correlation_index > 1 && ...
                maximum_correlation_index < numel(correlation_by_lag)

            left_correlation = correlation_by_lag(maximum_correlation_index-1);
            center_correlation = correlation_by_lag(maximum_correlation_index);
            right_correlation = correlation_by_lag(maximum_correlation_index+1);
            parabola_denominator = ...
                left_correlation - 2*center_correlation + right_correlation;

            if isfinite(parabola_denominator) && ...
                    abs(parabola_denominator) > eps

                lag_sample_offset = ...
                    0.5*(left_correlation - right_correlation) / ...
                    parabola_denominator;
                lag_sample_offset = max(min(lag_sample_offset, 1), -1);
                refined_lag_samples = refined_lag_samples + lag_sample_offset;
            end
        end

        Delay_s(section_index) = ...
            refined_lag_samples*mean_time_step_seconds;
    end

    Delay_s = Delay_s - Delay_s(reference_section_index);
    Delay_s(reference_section_index) = 0;
    relative_section_distance_m = ...
        section_distance_m(:) - section_distance_m(reference_section_index);

    [PWV_mps, fit_coefficients, R2] = fit_pwv_linear( ...
        relative_section_distance_m, Delay_s);
end


function [PWV_mps, Delay_s, fit_coefficients, R2, ...
    pair_delay_s, pair_distance_m, maximum_correlation] = ...
    estimate_pwv_xcor_urbina_pair( ...
        section_distance_m, time_vector_seconds, flow_rate_mL_per_s, ...
        section_I, section_IV)

    number_of_sections = size(flow_rate_mL_per_s, 1);
    Delay_s = nan(number_of_sections, 1);
    fit_coefficients = [NaN NaN];
    R2 = NaN;
    PWV_mps = NaN;
    pair_delay_s = NaN;
    pair_distance_m = NaN;
    maximum_correlation = NaN;

    if section_I < 1 || section_I > number_of_sections || ...
            section_IV < 1 || section_IV > number_of_sections
        warning('Las secciones XCor I-IV estan fuera del rango de cortes.')
        return
    end

    flow_I = flow_rate_mL_per_s(section_I,:)';
    flow_IV = flow_rate_mL_per_s(section_IV,:)';

    if all(~isfinite(flow_I)) || all(~isfinite(flow_IV))
        warning('XCor I-IV no tiene curvas de flujo validas.')
        return
    end

    flow_I = fillmissing(flow_I, 'linear', 'EndValues', 'nearest');
    flow_IV = fillmissing(flow_IV, 'linear', 'EndValues', 'nearest');

    flow_I = flow_I - mean(flow_I, 'omitnan');
    flow_IV = flow_IV - mean(flow_IV, 'omitnan');

    mean_time_step_seconds = mean(diff(time_vector_seconds));
    maximum_lag_seconds = 0.30;
    maximum_lag_samples = min( ...
        floor(numel(time_vector_seconds)/2), ...
        max(1, round(maximum_lag_seconds/mean_time_step_seconds)));
    lag_sample_vector = (-maximum_lag_samples:maximum_lag_samples)';
    correlation_by_lag = nan(numel(lag_sample_vector), 1);

    for lag_index = 1:numel(lag_sample_vector)
        current_lag_samples = lag_sample_vector(lag_index);

        if current_lag_samples >= 0
            reference_segment = flow_I(1:end-current_lag_samples);
            section_segment = flow_IV(1+current_lag_samples:end);
        else
            reference_segment = flow_I(1-current_lag_samples:end);
            section_segment = flow_IV(1:end+current_lag_samples);
        end

        valid_sample_mask = ...
            isfinite(reference_segment) & isfinite(section_segment);

        if nnz(valid_sample_mask) < 3
            continue
        end

        reference_segment = reference_segment(valid_sample_mask);
        section_segment = section_segment(valid_sample_mask);
        correlation_denominator = sqrt( ...
            sum(reference_segment.^2) * sum(section_segment.^2));

        if correlation_denominator > 0
            correlation_by_lag(lag_index) = ...
                sum(reference_segment .* section_segment) / ...
                correlation_denominator;
        end
    end

    if all(~isfinite(correlation_by_lag))
        warning('XCor I-IV no pudo calcular correlaciones validas.')
        return
    end

    [maximum_correlation, maximum_correlation_index] = ...
        max(correlation_by_lag);
    refined_lag_samples = lag_sample_vector(maximum_correlation_index);

    if maximum_correlation_index > 1 && ...
            maximum_correlation_index < numel(correlation_by_lag)

        left_correlation = correlation_by_lag(maximum_correlation_index-1);
        center_correlation = correlation_by_lag(maximum_correlation_index);
        right_correlation = correlation_by_lag(maximum_correlation_index+1);
        parabola_denominator = ...
            left_correlation - 2*center_correlation + right_correlation;

        if isfinite(parabola_denominator) && ...
                abs(parabola_denominator) > eps

            lag_sample_offset = ...
                0.5*(left_correlation - right_correlation) / ...
                parabola_denominator;
            lag_sample_offset = max(min(lag_sample_offset, 1), -1);
            refined_lag_samples = refined_lag_samples + lag_sample_offset;
        end
    end

    pair_delay_s = refined_lag_samples * mean_time_step_seconds;
    pair_distance_m = ...
        section_distance_m(section_IV) - section_distance_m(section_I);

    if ~isfinite(pair_delay_s) || pair_delay_s <= 0 || ...
            ~isfinite(pair_distance_m) || pair_distance_m <= 0
        warning('XCor I-IV entrego delay o distancia invalida.')
        return
    end

    PWV_mps = pair_distance_m / pair_delay_s;

    Delay_s(section_I) = 0;
    Delay_s(section_IV) = pair_delay_s;
end


function [PWV_mps, fit_coefficients, R2] = fit_pwv_linear( ...
    section_distance_m, characteristic_time_seconds)

    % Ajuste lineal comun a los tres metodos de PWV. Si t = a*s + b,
    % entonces PWV = 1/a. Se exige un minimo de tres cortes validos para
    % evitar ajustes degenerados.

    valid_sample_mask = ...
        isfinite(section_distance_m) & isfinite(characteristic_time_seconds);

    if nnz(valid_sample_mask) < 3
        PWV_mps = NaN;
        fit_coefficients = [NaN NaN];
        R2 = NaN;
        return
    end

    valid_section_distance_m = section_distance_m(valid_sample_mask);
    valid_characteristic_time_seconds = ...
        characteristic_time_seconds(valid_sample_mask);
    fit_coefficients = polyfit( ...
        valid_section_distance_m, valid_characteristic_time_seconds, 1);
    slope_seconds_per_meter = fit_coefficients(1);

    if ~isfinite(slope_seconds_per_meter) || ...
            abs(slope_seconds_per_meter) <= eps
        PWV_mps = NaN;
    else
        PWV_mps = 1/slope_seconds_per_meter;

        if PWV_mps < 0
            warning( ...
                'PWV negativa (%.6f m/s). Revise la orientacion.', ...
                PWV_mps)
        end
    end

    fitted_time_seconds = polyval(fit_coefficients, valid_section_distance_m);
    residual_sum_of_squares = sum( ...
        (valid_characteristic_time_seconds - fitted_time_seconds).^2);
    total_sum_of_squares = sum( ...
        (valid_characteristic_time_seconds - ...
        mean(valid_characteristic_time_seconds)).^2);

    if total_sum_of_squares > 0
        R2 = 1 - residual_sum_of_squares/total_sum_of_squares;
    else
        R2 = NaN;
    end
end


function plot_pwv_fit( ...
    section_distance_m, characteristic_time_seconds, fit_coefficients, ...
    PWV_mps, R2, figure_name, y_axis_label)

    % Grafica los tiempos caracteristicos validos y su ajuste lineal. Es una
    % herramienta de control visual para revisar dispersion, outliers y signo
    % de la pendiente antes de interpretar la PWV.

    valid_sample_mask = ...
        isfinite(section_distance_m(:)) & ...
        isfinite(characteristic_time_seconds(:));

    if nnz(valid_sample_mask) < 3
        warning('No hay suficientes cortes validos para graficar %s.', ...
            figure_name)
        return
    end

    valid_section_distance_m = section_distance_m(valid_sample_mask);
    valid_characteristic_time_seconds = ...
        characteristic_time_seconds(valid_sample_mask);
    fitted_distance_m = linspace( ...
        min(valid_section_distance_m), max(valid_section_distance_m), 200);
    fitted_time_seconds = polyval(fit_coefficients, fitted_distance_m);

    figure('Name', figure_name)
    plot(valid_section_distance_m, valid_characteristic_time_seconds, ...
        'o', 'MarkerSize', 7, 'LineWidth', 1.5)
    hold on
    plot(fitted_distance_m, fitted_time_seconds, '-', 'LineWidth', 1.8)
    xlabel('Distancia sobre centerline [m]')
    ylabel(y_axis_label)
    title(sprintf('%s = %.3f m/s | R^2 = %.3f', ...
        figure_name, PWV_mps, R2))
    legend({'Tiempo por corte', 'Ajuste lineal'}, 'Location', 'best')
    grid on
    hold off
end
