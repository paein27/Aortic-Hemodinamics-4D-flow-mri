function data = load_pwv_dataset(config)
%LOAD_PWV_DATASET Carga un dataset 4D Flow ya procesado y normaliza unidades.
%   Lee malla FE, campo de Laplace, velocidad, centerline, area y eje
%   temporal desde la carpeta config.data_folder. Detecta y convierte las
%   coordenadas a mm y la velocidad a m/s. Devuelve todo en la estructura
%   'data', que alimenta al resto del pipeline.

    normal_rest_data_folder = config.data_folder;
    dataset_label = erase(string(normal_rest_data_folder), "\MATLAB FILES");
    [~, dataset_label] = fileparts(dataset_label);

    % Carga explicita por struct: dentro de una funcion garantiza que estos
    % nombres sean variables (no funciones del path) y evita colisiones.
    mesh_elem = load(fullfile(normal_rest_data_folder, 'FE Mesh', 'elem.mat'));
    mesh_nodes = load(fullfile(normal_rest_data_folder, 'FE Mesh', 'nodes.mat'));
    mesh_faces = load(fullfile(normal_rest_data_folder, 'FE Mesh', 'faces.mat'));
    laplace_file = load(fullfile(normal_rest_data_folder, 'FE Laplace', 'Laplace.mat'));
    velocity_file = load(fullfile(normal_rest_data_folder, 'FE Velocity', 'VEL.mat'));
    centerline_file = load(fullfile(normal_rest_data_folder, 'FE Centerline', 'Centerline.mat'));
    area_file = load(fullfile(normal_rest_data_folder, 'FE Area', 'Area.mat'));

    elem = mesh_elem.elem;
    nodes = mesh_nodes.nodes;
    faces = mesh_faces.faces;
    Laplace = laplace_file.Laplace;
    VEL = velocity_file.VEL;
    Centerline = centerline_file.Centerline;
    Area = area_file.Area;

    [coordinate_scale_to_mm, detected_coordinate_units] = ...
        detect_coordinate_scale_to_mm(nodes, Centerline);
    nodes = nodes * coordinate_scale_to_mm;
    Centerline = Centerline * coordinate_scale_to_mm;

    fprintf('Unidad de coordenadas detectada: %s\n', detected_coordinate_units)

    surface_node_ids = unique(faces(:));

    time_file = fullfile(normal_rest_data_folder, '2D Flow', 'Time.mat');

    number_of_time_phases = size(VEL, 3);

    [time_vector_seconds, time_vector_source] = ...
        get_time_vector_for_pwv(time_file, number_of_time_phases);

    fprintf('Numero de fases temporales: %d\n', number_of_time_phases)

    % ---- Deteccion y conversion de unidades de velocidad ----
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

    data = struct();
    data.dataset_label = dataset_label;
    data.elem = elem;
    data.nodes = nodes;
    data.faces = faces;
    data.Laplace = Laplace;
    data.Centerline = Centerline;
    data.Area = Area;
    data.surface_node_ids = surface_node_ids;
    data.velocity_nodes_mps = velocity_nodes_mps;
    data.number_of_time_phases = number_of_time_phases;
    data.time_vector_seconds = time_vector_seconds;
    data.time_vector_source = time_vector_source;
    data.detected_coordinate_units = detected_coordinate_units;
    data.detected_velocity_units = detected_velocity_units;
end
