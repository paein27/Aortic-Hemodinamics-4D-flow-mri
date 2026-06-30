clc
clear

%%%%%%%%%%%%%%%%%%%%%%%% CARGA DE DATOS PROCESADOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile('Data','Datos_Procesados','LRFFE','MATLAB FILES','Segmentation ROI','SEG.mat'))
load(fullfile('Data','Datos_Procesados','LRFFE','MATLAB FILES','FE Centerline','Centerline.mat'))
load(fullfile('Data','Datos_Procesados','LRFFE','MATLAB FILES','Segmentation ROI','VoxelSize.mat'))
load(fullfile('Data','Datos_Procesados','LRFFE','MATLAB FILES','FE Mesh','nodes.mat'))
load(fullfile('Data','Datos_Procesados','LRFFE','MATLAB FILES','FE Velocity','VEL.mat'))
time_file = fullfile('Data','Datos_Procesados','LRFFE','MATLAB FILES','2D Flow','Time.mat');


%%%%%%%%%%%%%%% CREACION DE GRILLA EN COORDENADAS REALES [m] %%%%%%%%%%%%%%%%%%%%%%
[grid_idx_x, grid_idx_y, grid_idx_z] = meshgrid(0:size(SEG,2)-1, 0:size(SEG,1)-1, 0:size(SEG,3)-1);
grid_X_m = grid_idx_x * VoxelSize(2)/1000;
grid_Y_m = grid_idx_y * VoxelSize(1)/1000;
grid_Z_m = grid_idx_z * VoxelSize(3)/1000;

%% VISUALIZACION 3D

fig_aorta = figure('Name','Aorta, planos y velocidades');

aorta_patch = patch(isosurface(grid_X_m, grid_Y_m, grid_Z_m, SEG, 0.5));

set(aorta_patch, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none', 'FaceAlpha',0.25)

camlight
lighting gouraud

hold on

plot3(Centerline(:,2), Centerline(:,1), Centerline(:,3), 'r', 'LineWidth',2)

axis equal
grid on
view(3)

%%%%%%%%%%%%%%%%%%%%%% GENERAR 30 PLANOS EQUIDISTANTES + INTERPOLAR VEL + PROYECTAR NORMAL %%%%%%%%%%%%%%%%%%%%%%

n_planes = 30;

% Tamaño del plano: 0.020 m hacia cada lado => plano total de 4 cm x 4 cm
plane_half_size_m = 0.020; % [m]

% Resolución del plano: grilla de 41x41 puntos en cada plano
n_grid = 41;

% Paso entre puntos de la grilla
grid_step_m = (2*plane_half_size_m) / n_grid;

% Puntos tipo centro de celda (no en los bordes)
grid_offsets_m = linspace(-plane_half_size_m + grid_step_m/2, plane_half_size_m - grid_step_m/2, n_grid);

% Matrices de offset 2D para construir los planos
[offset_a, offset_b] = meshgrid(grid_offsets_m, grid_offsets_m);

% Área asociada a cada punto de muestreo del plano
cell_area_m2 = grid_step_m^2; % [m^2]

% Reordena las coordenadas del centerline y nodos a orden X Y Z
centerline_xyz = [Centerline(:,2), Centerline(:,1), Centerline(:,3)];
nodes_xyz      = [nodes(:,2),      nodes(:,1),      nodes(:,3)];


n_phases = size(VEL,3);

% -------------------------------------------------------------------------
% Eje temporal para cálculo de PWV
% -------------------------------------------------------------------------

[time_s, time_source] = get_time_vector_for_pwv(time_file, n_phases);

fprintf('Número de fases temporales: %d\n', numel(time_s))
fprintf('Tiempo inicial: %.6f s\n', time_s(1))
fprintf('Tiempo final: %.6f s\n', time_s(end))
fprintf('dt medio: %.6f s\n', mean(diff(time_s)))

% Detectar unidades de VEL:
% Si el percentil 99 de la magnitud supera 20, están en cm/s; si no, en m/s.
vel_magnitude        = abs(VEL(:));
vel_magnitude        = vel_magnitude(isfinite(vel_magnitude));
vel_magnitude_sorted = sort(vel_magnitude);
idx_p99              = max(1, round(0.99*numel(vel_magnitude_sorted)));
vel_p99_val          = vel_magnitude_sorted(idx_p99);

if vel_p99_val > 20
    vel_scale_mps = 0.01;   % cm/s -> m/s
    vel_units     = 'cm/s';
else
    vel_scale_mps = 1.0;    % ya está en m/s
    vel_units     = 'm/s';
end

fprintf('\nUnidad de VEL detectada: %s\n', vel_units)
fprintf('Factor de conversión a m/s: %.4f\n', vel_scale_mps)

% Velocidad en los nodos FE en orden X Y Z [m/s]
Vx_nodes = squeeze(VEL(:,2,:)) * vel_scale_mps;
Vy_nodes = squeeze(VEL(:,1,:)) * vel_scale_mps;
Vz_nodes = squeeze(VEL(:,3,:)) * vel_scale_mps;

% -------------------------------------------------------------------------
% Calcular posiciones equidistantes en distancia física sobre la centerline
% -------------------------------------------------------------------------

% Distancia entre puntos consecutivos del centerline
delta_cl_m = sqrt(sum(diff(centerline_xyz,1,1).^2, 2));

% Eliminar puntos repetidos o muy cercanos
valid_cl_pts  = [true; delta_cl_m > eps];
centerline_xyz = centerline_xyz(valid_cl_pts, :);

% Recalcular distancia acumulada (longitud de arco)
delta_cl_m   = sqrt(sum(diff(centerline_xyz,1,1).^2, 2)); % [m]
arc_length_m = [0; cumsum(delta_cl_m)];                   % [m]

L_total = arc_length_m(end);

% Crear 30 planos equidistantes, evitando los extremos del centerline,
% donde la tangente local es menos estable.
s_planes = linspace(0, L_total, n_planes+2);
s_planes = s_planes(2:end-1);

fprintf('Número de planos: %d\n', n_planes)
fprintf('Longitud total de la centerline: %.1f mm\n', L_total*1000)
fprintf('Separación entre planos: %.1f mm\n', mean(diff(s_planes))*1000)

% Centros de los planos interpolados sobre el centerline
plane_centers = interp1(arc_length_m, centerline_xyz, s_planes, 'linear');

% Tangente local del centerline en función de la longitud de arco
Tan_x = gradient(centerline_xyz(:,1), arc_length_m);
Tan_y = gradient(centerline_xyz(:,2), arc_length_m);
Tan_z = gradient(centerline_xyz(:,3), arc_length_m);

tangent_cl = [Tan_x, Tan_y, Tan_z];

% Normales de los planos = tangente del centerline en cada posición
plane_normals = interp1(arc_length_m, tangent_cl, s_planes, 'linear');
plane_normals = plane_normals ./ vecnorm(plane_normals, 2, 2);

% -------------------------------------------------------------------------
% Crear interpoladores de velocidad desde nodos FE hacia puntos arbitrarios
% -------------------------------------------------------------------------

interp_Vx = scatteredInterpolant(nodes_xyz(:,1), nodes_xyz(:,2), nodes_xyz(:,3), zeros(size(nodes_xyz,1),1), 'linear','none');
interp_Vy = scatteredInterpolant(nodes_xyz(:,1), nodes_xyz(:,2), nodes_xyz(:,3), zeros(size(nodes_xyz,1),1), 'linear','none');
interp_Vz = scatteredInterpolant(nodes_xyz(:,1), nodes_xyz(:,2), nodes_xyz(:,3), zeros(size(nodes_xyz,1),1), 'linear','none');

% -------------------------------------------------------------------------
% Prealocar outputs
% -------------------------------------------------------------------------

Vx_planes  = nan(n_grid, n_grid, n_phases, n_planes); % [m/s]
Vy_planes  = nan(n_grid, n_grid, n_phases, n_planes); % [m/s]
Vz_planes  = nan(n_grid, n_grid, n_phases, n_planes); % [m/s]
Vn_planes  = nan(n_grid, n_grid, n_phases, n_planes); % [m/s] componente normal

Lumen_masks  = false(n_grid, n_grid, n_planes);

Mean_Vn    = nan(n_planes, n_phases); % [m/s]
Peak_Vn    = nan(n_planes, n_phases); % [m/s]

Area_planes = nan(n_planes, 1);       % [m^2]
Flow_m3s    = nan(n_planes, n_phases); % [m^3/s]
Flow_mLs    = nan(n_planes, n_phases); % [mL/s]

Planes = struct();

% -------------------------------------------------------------------------
% Loop de procesamiento sobre planos equidistantes
% -------------------------------------------------------------------------

for k = 1:n_planes

    %% Centro y normal del plano k

    plane_center = plane_centers(k,:);
    plane_normal = plane_normals(k,:);
    plane_normal = plane_normal / norm(plane_normal);

    %% Vector auxiliar para construir base ortonormal del plano

    aux_vec = [0 0 1];

    if abs(dot(plane_normal, aux_vec)) > 0.9
        aux_vec = [0 1 0];
    end

    %% Base ortonormal del plano (dos vectores perpendiculares a la normal)

    plane_u1 = cross(plane_normal, aux_vec);
    plane_u1 = plane_u1 / norm(plane_u1);

    plane_u2 = cross(plane_normal, plane_u1);
    plane_u2 = plane_u2 / norm(plane_u2);

    %% Coordenadas 3D de los puntos de la grilla del plano

    plane_X = plane_center(1) + plane_u1(1)*offset_a + plane_u2(1)*offset_b;
    plane_Y = plane_center(2) + plane_u1(2)*offset_a + plane_u2(2)*offset_b;
    plane_Z = plane_center(3) + plane_u1(3)*offset_a + plane_u2(3)*offset_b;

    %% Máscara del lumen en este plano usando SEG

    seg_on_plane = interp3(grid_X_m, grid_Y_m, grid_Z_m, double(SEG), plane_X, plane_Y, plane_Z, 'nearest', 0);
    lumen_mask   = seg_on_plane > 0.5;

    Lumen_masks(:,:,k) = lumen_mask;

    % Área aproximada del lumen en este plano
    Area_planes(k) = nnz(lumen_mask) * cell_area_m2; % [m^2]

    %% Dibujar plano

    surf(plane_X, plane_Y, plane_Z, ...
        'FaceAlpha', 0.20, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'g')

    %% Dibujar normal del plano

    quiver3(plane_center(1), plane_center(2), plane_center(3), ...
            plane_normal(1)*0.01, ...
            plane_normal(2)*0.01, ...
            plane_normal(3)*0.01, ...
            'b', 'LineWidth', 2)

    %% Interpolar velocidades al plano y proyectar en dirección normal

    for t = 1:n_phases

        % Actualizar valores de velocidad en el interpolador para la fase t
        interp_Vx.Values = Vx_nodes(:,t);
        interp_Vy.Values = Vy_nodes(:,t);
        interp_Vz.Values = Vz_nodes(:,t);

        % Interpolar velocidades desde nodos FE hacia los puntos del plano
        Vx_plane = interp_Vx(plane_X, plane_Y, plane_Z);
        Vy_plane = interp_Vy(plane_X, plane_Y, plane_Z);
        Vz_plane = interp_Vz(plane_X, plane_Y, plane_Z);

        % Componente normal de la velocidad: Vn = V · normal_plano
        Vn_plane = Vx_plane*plane_normal(1) + Vy_plane*plane_normal(2) + Vz_plane*plane_normal(3);

        % Enmascarar puntos fuera del lumen
        Vx_plane(~lumen_mask) = NaN;
        Vy_plane(~lumen_mask) = NaN;
        Vz_plane(~lumen_mask) = NaN;
        Vn_plane(~lumen_mask) = NaN;

        % Guardar resultados
        Vx_planes(:,:,t,k) = Vx_plane;
        Vy_planes(:,:,t,k) = Vy_plane;
        Vz_planes(:,:,t,k) = Vz_plane;
        Vn_planes(:,:,t,k) = Vn_plane;

        % -------------------------------------------------------------
        % Velocidad normal media y flujo volumétrico en el plano
        % -------------------------------------------------------------

        Vn_in_lumen = Vn_plane(lumen_mask);
        Vn_in_lumen = Vn_in_lumen(isfinite(Vn_in_lumen));

        if ~isempty(Vn_in_lumen)

            Mean_Vn(k,t) = mean(Vn_in_lumen);       % velocidad media [m/s]
            Peak_Vn(k,t) = max(Vn_in_lumen);         % velocidad pico [m/s]

            % Flujo: Q = integral(Vn dA) ≈ sum(Vn_i * cell_area)
            Flow_m3s(k,t) = sum(Vn_in_lumen) * cell_area_m2; % [m^3/s]
            Flow_mLs(k,t) = Flow_m3s(k,t) * 1e6;             % [mL/s]

        end

    end

    %% Guardar geometría del plano

    Planes(k).center = plane_center;
    Planes(k).normal = plane_normal;
    Planes(k).u1     = plane_u1;
    Planes(k).u2     = plane_u2;
    Planes(k).X      = plane_X;
    Planes(k).Y      = plane_Y;
    Planes(k).Z      = plane_Z;
    Planes(k).s      = s_planes(k);
    Planes(k).mask   = lumen_mask;

end


%% Visualizar velocidades interpoladas en los planos usando quiver3

% Usar la fase cardíaca con mayor flujo medio positivo para la visualización
mean_flow_per_phase = mean(Flow_mLs, 1, 'omitnan');

if any(mean_flow_per_phase > 0)
    [~, phase_vis] = max(mean_flow_per_phase);
else
    [~, phase_vis] = max(abs(mean_flow_per_phase));
    warning(['No se encontró una fase con flujo medio positivo. ', ...
             'La visualización puede mostrar flujo reverso o una orientación ', ...
             'invertida del centerline/normales.'])
end

fprintf('Flujo medio de la fase visualizada: %.6f mL/s\n', mean_flow_per_phase(phase_vis))

figure(fig_aorta)
hold on

quiver_step  = 4;     % muestra 1 de cada 4 vectores para no saturar
quiver_scale = 0.008; % [m por cada m/s]

for k = 1:n_planes

    plane_X    = Planes(k).X;
    plane_Y    = Planes(k).Y;
    plane_Z    = Planes(k).Z;
    lumen_mask = Planes(k).mask;

    Vn_plane    = Vn_planes(:,:,phase_vis,k);
    plane_normal = Planes(k).normal;
    plane_normal = plane_normal / norm(plane_normal);

    % Descomponer velocidad normal en componentes vectoriales
    Vn_vec_x = Vn_plane * plane_normal(1);
    Vn_vec_y = Vn_plane * plane_normal(2);
    Vn_vec_z = Vn_plane * plane_normal(3);

    % Submuestreo espacial para no saturar la figura
    [row_idx, col_idx] = ndgrid(1:size(lumen_mask,1), 1:size(lumen_mask,2));

    quiver_mask = lumen_mask & ...
                  mod(row_idx-1, quiver_step) == 0 & ...
                  mod(col_idx-1, quiver_step) == 0 & ...
                  isfinite(Vn_plane);

    quiver3( ...
        plane_X(quiver_mask), plane_Y(quiver_mask), plane_Z(quiver_mask), ...
        Vn_vec_x(quiver_mask)*quiver_scale, ...
        Vn_vec_y(quiver_mask)*quiver_scale, ...
        Vn_vec_z(quiver_mask)*quiver_scale, ...
        0, 'k', 'LineWidth', 0.8)

end

title(sprintf('Aorta, planos y velocidades proyectadas | fase %d', phase_vis))

axis equal
grid on
view(3)
hold off

%% Resumen de área y flujo por plano

PeakFlow_mLs = nan(n_planes,1);
MeanFlow_mLs = nan(n_planes,1);

for k = 1:n_planes
    flow_curve_k = Flow_mLs(k,:);
    flow_curve_k = flow_curve_k(isfinite(flow_curve_k));

    if ~isempty(flow_curve_k)
        PeakFlow_mLs(k) = max(flow_curve_k);
        MeanFlow_mLs(k) = mean(flow_curve_k);
    end
end

FlowSummary = table( ...
    (1:n_planes)', ...
    s_planes(:), ...
    Area_planes, ...
    Area_planes*1e6, ...
    PeakFlow_mLs, ...
    MeanFlow_mLs, ...
    'VariableNames',{'Plane','s_m','Area_m2','Area_mm2','PeakFlow_mLs','MeanFlow_mLs'} ...
);

disp(FlowSummary)

%% Graficar curvas de flujo por plano

figure('Name','Curvas de flujo por plano')

plot(time_s, Flow_mLs', 'LineWidth',1.5)

xlabel('Tiempo [s]')
ylabel('Flujo [mL/s]')
title('Curvas de flujo calculadas en planos perpendiculares a la aorta')
legend(compose('Plano %d', 1:n_planes), 'Location','eastoutside')
grid on

%% Calcular PWV con métodos de Markl: TTP, TTF y XCor

flow_for_pwv_mLs = Flow_mLs;

% -------------------------------------------------------------------------
% Método 1: TTP = Time To Peak
% -------------------------------------------------------------------------

[PWV_TTP_mps, TTP_s, fit_TTP, R2_TTP] = estimate_pwv_ttp( ...
    s_planes(:), time_s, flow_for_pwv_mLs);

% -------------------------------------------------------------------------
% Método 2: TTF = Time To Foot
% -------------------------------------------------------------------------

[PWV_TTF_mps, TTF_s, fit_TTF, R2_TTF] = estimate_pwv_ttf( ...
    s_planes(:), time_s, flow_for_pwv_mLs);

% -------------------------------------------------------------------------
% Método 3: XCor = Cross-correlation
% -------------------------------------------------------------------------

ref_plane = 1;

[PWV_XCor_mps, Delay_XCor_s, fit_XCor, R2_XCor] = estimate_pwv_xcor( ...
    s_planes(:), time_s, flow_for_pwv_mLs, ref_plane);

% -------------------------------------------------------------------------
% Conteo de planos válidos por método
% -------------------------------------------------------------------------

n_valid_TTP  = nnz(isfinite(TTP_s));
n_valid_TTF  = nnz(isfinite(TTF_s));
n_valid_XCor = nnz(isfinite(Delay_XCor_s));

% -------------------------------------------------------------------------
% Tabla resumen de PWV
% -------------------------------------------------------------------------

PWVSummary = table( ...
    ["TTP";"TTF";"XCor"], ...
    [PWV_TTP_mps; PWV_TTF_mps; PWV_XCor_mps], ...
    [fit_TTP(1); fit_TTF(1); fit_XCor(1)], ...
    [R2_TTP; R2_TTF; R2_XCor], ...
    [n_valid_TTP; n_valid_TTF; n_valid_XCor], ...
    repmat(n_planes,3,1), ...
    'VariableNames',{'Method','PWV_mps','Slope_s_per_m','R2','ValidPlanes','TotalPlanes'} ...
);

% -------------------------------------------------------------------------
% Tabla de tiempos característicos por plano
% -------------------------------------------------------------------------

PWVPlaneSummary = table( ...
    (1:n_planes)', ...
    s_planes(:), ...
    TTP_s(:), ...
    TTF_s(:), ...
    Delay_XCor_s(:), ...
    'VariableNames',{'Plane','s_m','TTP_s','TTF_s','XCorDelay_s'} ...
);

% -------------------------------------------------------------------------
% Impresión por consola
% -------------------------------------------------------------------------

fprintf('\n============================================================\n')
fprintf('RESULTADOS PWV - METODOS DE MARKL\n')
fprintf('============================================================\n')

fprintf('Fuente temporal: %s\n', time_source)
fprintf('Numero de fases temporales: %d\n', numel(time_s))
fprintf('Tiempo inicial: %.6f s\n', time_s(1))
fprintf('Tiempo final: %.6f s\n', time_s(end))
fprintf('dt medio: %.6f s\n', mean(diff(time_s)))

fprintf('\nNumero de planos: %d\n', n_planes)
fprintf('Plano de referencia XCor: %d\n', ref_plane)
fprintf('Distancia total centerline: %.1f mm\n', L_total*1000)
fprintf('Separacion entre planos: %.1f mm\n', mean(diff(s_planes))*1000)

fprintf('\nResumen global:\n')
fprintf('Metodo     PWV [m/s]        Pendiente [s/m]       R2          Planos validos\n')
fprintf('------     ---------        ----------------      -------     --------------\n')

fprintf('TTP        %12.6f     %16.8f     %8.4f     %d/%d\n', ...
    PWV_TTP_mps,  fit_TTP(1),  R2_TTP,  n_valid_TTP,  n_planes)

fprintf('TTF        %12.6f     %16.8f     %8.4f     %d/%d\n', ...
    PWV_TTF_mps,  fit_TTF(1),  R2_TTF,  n_valid_TTF,  n_planes)

fprintf('XCor       %12.6f     %16.8f     %8.4f     %d/%d\n', ...
    PWV_XCor_mps, fit_XCor(1), R2_XCor, n_valid_XCor, n_planes)

fprintf('\nTabla resumen PWV:\n')
disp(PWVSummary)

fprintf('\nTiempos caracteristicos por plano:\n')
disp(PWVPlaneSummary)

%% Visualización del ajuste lineal — TTP

valid_planes_ttp = isfinite(s_planes(:)) & isfinite(TTP_s(:));

if nnz(valid_planes_ttp) < 3
    warning('No hay suficientes planos válidos para graficar el ajuste TTP.')
else
    dist_valid_ttp = s_planes(valid_planes_ttp);
    ttp_values     = TTP_s(valid_planes_ttp);

    dist_fit    = linspace(min(dist_valid_ttp), max(dist_valid_ttp), 200);
    ttp_fit_line = polyval(fit_TTP, dist_fit);

    figure('Name','PWV por TTP interpolado')
    plot(dist_valid_ttp, ttp_values, 'o', 'MarkerSize',7, 'LineWidth',1.5)
    hold on
    plot(dist_fit, ttp_fit_line, '-', 'LineWidth',1.8)
    xlabel('Distancia sobre centerline [m]')
    ylabel('TTP interpolado [s]')
    title(sprintf('PWV por TTP = %.3f m/s | R^2 = %.3f', PWV_TTP_mps, R2_TTP))
    legend({'TTP por plano','Ajuste lineal'}, 'Location','best')
    grid on
    hold off
end

%% Visualización del ajuste lineal — TTF

valid_planes_ttf = isfinite(s_planes(:)) & isfinite(TTF_s(:));

if nnz(valid_planes_ttf) < 3
    warning('No hay suficientes planos válidos para graficar el ajuste TTF.')
else
    dist_valid_ttf = s_planes(valid_planes_ttf);
    ttf_values     = TTF_s(valid_planes_ttf);

    dist_fit     = linspace(min(dist_valid_ttf), max(dist_valid_ttf), 200);
    ttf_fit_line = polyval(fit_TTF, dist_fit);

    figure('Name','PWV por TTF')
    plot(dist_valid_ttf, ttf_values, 'o', 'MarkerSize',7, 'LineWidth',1.5)
    hold on
    plot(dist_fit, ttf_fit_line, '-', 'LineWidth',1.8)
    xlabel('Distancia sobre centerline [m]')
    ylabel('TTF [s]')
    title(sprintf('PWV por TTF = %.3f m/s | R^2 = %.3f', PWV_TTF_mps, R2_TTF))
    legend({'TTF por plano','Ajuste lineal'}, 'Location','best')
    grid on
    hold off
end

%% Visualización del ajuste lineal — XCor

valid_planes_xcor = isfinite(s_planes(:)) & isfinite(Delay_XCor_s(:));

if nnz(valid_planes_xcor) < 3
    warning('No hay suficientes planos válidos para graficar el ajuste XCor.')
else
    dist_valid_xcor = s_planes(valid_planes_xcor);
    xcor_values     = Delay_XCor_s(valid_planes_xcor);

    dist_fit      = linspace(min(dist_valid_xcor), max(dist_valid_xcor), 200);
    xcor_fit_line = polyval(fit_XCor, dist_fit);

    figure('Name','PWV por XCor')
    plot(dist_valid_xcor, xcor_values, 'o', 'MarkerSize',7, 'LineWidth',1.5)
    hold on
    plot(dist_fit, xcor_fit_line, '-', 'LineWidth',1.8)
    xlabel('Distancia sobre centerline [m]')
    ylabel('Delay XCor relativo [s]')
    title(sprintf('PWV por XCor = %.3f m/s | R^2 = %.3f', PWV_XCor_mps, R2_XCor))
    legend({'Delay XCor por plano','Ajuste lineal'}, 'Location','best')
    grid on
    hold off
end

fprintf('============================================================\n')

%% Control rápido de equidistancia

fprintf('\nDistancia total centerline: %.1f mm\n', L_total*1000)
fprintf('Separacion entre planos: %.1f mm\n', mean(diff(s_planes))*1000)
fprintf('Numero de planos generados: %d\n', n_planes)

fprintf('size(Vn_planes) = %s\n', mat2str(size(Vn_planes)))
fprintf('size(Mean_Vn)   = %s\n', mat2str(size(Mean_Vn)))

% =========================================================================
%                           FUNCIONES AUXILIARES
% =========================================================================

function [time_s, time_source] = get_time_vector_for_pwv(time_file, n_phases)

    time_s      = [];
    time_source = '';

    if exist(time_file,'file') == 2

        S     = load(time_file);
        names = fieldnames(S);

        for i = 1:numel(names)

            x = S.(names{i});

            if isnumeric(x) && ~isempty(x)

                x = squeeze(x);

                if isvector(x) && numel(x) >= n_phases
                    time_s      = double(x(:));
                    time_s      = time_s(1:n_phases);
                    time_source = ['Time.mat -> variable: ', names{i}];
                    break
                end

            end

        end

    end

    if isempty(time_s)

        default_heart_rate_bpm = 68;

        heart_rate_input = input(sprintf( ...
            '\nNo se pudo usar Time.mat. Ingrese frecuencia cardíaca [bpm] [default = %g]: ', ...
            default_heart_rate_bpm));

        if isempty(heart_rate_input)
            heart_rate_bpm = default_heart_rate_bpm;
        elseif isnumeric(heart_rate_input) && ...
               isscalar(heart_rate_input) && ...
               isfinite(heart_rate_input) && ...
               heart_rate_input > 0
            heart_rate_bpm = heart_rate_input;
        else
            error('La frecuencia cardíaca debe ser un número positivo.')
        end

        T_cycle_s = 60 / heart_rate_bpm;
        time_s    = linspace(0, T_cycle_s, n_phases+1)';
        time_s    = time_s(1:end-1);

        time_source = sprintf('Estimado desde heart_rate = %.2f bpm', heart_rate_bpm);

    else

        time_s = time_s(:);

        % Convertir a segundos si parece estar en milisegundos
        if max(abs(time_s)) > 10
            time_s = time_s / 1000;
        end

        % Llevar el inicio a cero
        time_s = time_s - time_s(1);

    end

    if numel(time_s) ~= n_phases
        error('El eje temporal tiene %d puntos, pero VEL tiene %d fases.', numel(time_s), n_phases)
    end

    if any(~isfinite(time_s)) || any(diff(time_s) <= 0)
        error('El eje temporal para PWV debe ser finito y estrictamente creciente.')
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%% FUNCION CALCULO TTP %%%%%%%%%%%%%%%%%%%%%%%%%
% estimate_pwv_ttp  Estima la PWV con el método Time-To-Peak (TTP).
%
%   Entradas:
%     s_m     [n_planes x 1]           Posición de cada plano [m]
%     time_s  [n_phases x 1]           Eje temporal del ciclo cardíaco [s]
%     Q_mLs   [n_planes x n_phases]    Curvas de flujo volumétrico [mL/s]
%
%   Salidas:
%     PWV_mps   escalar          PWV estimada [m/s]  (1/pendiente_ajuste)
%     TTP_s     [n_planes x 1]  TTP de cada plano [s]; NaN si no válido
%     fit_coeff [1 x 2]          Coeficientes del ajuste lineal
%     R2        escalar          Coeficiente de determinación del ajuste

function [PWV_mps, TTP_s, fit_coeff, R2] = estimate_pwv_ttp(s_m, time_s, Q_mLs)

    % --- Inicialización ---------------------------------------------------
    n_planes = size(Q_mLs, 1);   % número de planos de corte
    TTP_s    = nan(n_planes, 1); % preasignar TTP como NaN para cada plano

    % Paso de tiempo medio [s]; se asume muestreo aproximadamente uniforme.
    % Se usa para convertir el índice de peak (entero o fraccionario tras
    % refinamiento parabólico) a tiempo absoluto en segundos.
    dt_s = mean(diff(time_s));

    % --- Bucle sobre planos -----------------------------------------------
    for k = 1:n_planes

        % Extraer la curva de flujo del plano k como vector columna
        flow_curve = Q_mLs(k,:)';

        % Descartar planos sin ningún dato válido (todo NaN o Inf).
        % Si se continuara con all-NaN, max() devolvería NaN y los cálculos
        % posteriores serían inválidos.
        if all(~isfinite(flow_curve))
            continue
        end

        % Rellenar posibles NaN internos con interpolación lineal.
        % Los extremos se rellenan con el valor más cercano válido ('nearest').
        % Esto garantiza que max() y las operaciones aritméticas funcionen
        % sin sesgos por huecos en la señal.
        flow_curve = fillmissing(flow_curve, 'linear', 'EndValues','nearest');

        % Localizar el índice entero del máximo de flujo sistólico
        [~, idx_peak] = max(flow_curve);

        % Convertir a índice base-0 para la aritmética temporal:
        idx_peak_float = idx_peak - 1;

        % --- Refinamiento parabólico del peak ------------------------------
        % La curva de flujo es discreta; el verdadero peak puede caer entre muestras.
        if idx_peak > 1 && idx_peak < numel(flow_curve)

            q_left   = flow_curve(idx_peak-1); % muestra anterior al peak
            q_center = flow_curve(idx_peak);   % muestra del peak discreto
            q_right  = flow_curve(idx_peak+1); % muestra posterior al peak

            % Denominador de la fórmula del vértice parabólico
            parab_denom = q_left - 2*q_center + q_right;

            % Solo aplicar refinamiento si el denominador es numéricamente
            % significativo (denominador < 0 garantiza concavidad hacia abajo)
            if isfinite(parab_denom) && abs(parab_denom) > eps
                % Desplazamiento fraccionario respecto al índice del peak
                delta = 0.5*(q_left - q_right) / parab_denom;

                % Limitar delta a [-1, 1] para no saltar más de un paso temporal
                delta = max(min(delta,1), -1);

                % Actualizar el índice fraccionario incluyendo el ajuste
                idx_peak_float = idx_peak - 1 + delta;
            end

        end

        % Convertir el índice fraccionario a tiempo absoluto [s]
        TTP_s(k) = time_s(1) + idx_peak_float * dt_s;

    end

    % --- Ajuste lineal TTP vs distancia → PWV ----------------------------
    % Llama a fit_pwv_linear que hace polyfit(s, TTP, 1) y devuelve
    % PWV = 1/pendiente, los coeficientes del polinomio y el R2.
    [PWV_mps, fit_coeff, R2] = fit_pwv_linear(s_m, TTP_s);

end


%%%%%%%%%%%%%%%%%%%%%%%%% FUNCION CALCULO TTF %%%%%%%%%%%%%%%%%%%%%%%%%

function [PWV_mps, TTF_s, fit_coeff, R2] = estimate_pwv_ttf(s_m, time_s, Q_mLs)

    n_planes = size(Q_mLs, 1);
    TTF_s    = nan(n_planes, 1);

    % --- Bucle sobre planos -----------------------------------------------
    for k = 1:n_planes

        % Extraer curva de flujo del plano k como vector columna
        flow_curve = Q_mLs(k,:)';

        % Descartar plano si no tiene ningún dato válido
        if all(~isfinite(flow_curve))
            continue
        end

        % Rellenar NaN internos con interpolación lineal (igual que en TTP)
        flow_curve = fillmissing(flow_curve, 'linear', 'EndValues','nearest');

        % --- Localizar el peak de flujo -----------------------------------
        [Qpeak, idx_peak] = max(flow_curve);

        % Validaciones previas antes de buscar el pie de la onda:
        %   · idx_peak < 3  -> el peak está demasiado al inicio del ciclo;
        %                      no hay suficiente rampa de subida para ajustar.
        %   · Qpeak <= 0    -> no hay flujo sistólico positivo detectado;
        %                      la curva podría ser toda diastólica o ruidosa.
        if idx_peak < 3 || ~isfinite(Qpeak) || Qpeak <= 0
            continue
        end

        % --- Definir umbrales del tramo de subida (20 %-80 % del peak) ---
        thr20 = 0.20 * Qpeak; % umbral inferior: 20 % del flujo pico [mL/s]
        thr80 = 0.80 * Qpeak; % umbral superior: 80 % del flujo pico [mL/s]

        % Primer índice en la rampa ascendente donde la curva supera el 20 %
        idx20 = find(flow_curve(1:idx_peak) >= thr20, 1, 'first');

        % Primer índice en la rampa ascendente donde la curva supera el 80 %
        idx80 = find(flow_curve(1:idx_peak) >= thr80, 1, 'first');

        % Validar que ambos umbrales existen y que idx80 > idx20
        if isempty(idx20) || isempty(idx80) || idx80 <= idx20
            continue
        end

        % --- Definir el tramo de subida para el ajuste lineal -------------
        % Se expande levemente el rango [idx20-1, idx80+1] para incluir un
        % punto extra en cada extremo, mejorando la estabilidad del polyfit.
        % Se clampea al rango válido [1, idx_peak] para no salirse del array.
        idx_upslope = (max(1, idx20-1) : min(idx_peak, idx80+1))';

        % Necesitar al menos 2 puntos para ajustar una recta (polyfit grado 1)
        if numel(idx_upslope) < 2
            continue
        end

        % --- Ajuste lineal sobre el tramo de subida -----------------------
        slope_fit = polyfit(time_s(idx_upslope), flow_curve(idx_upslope), 1);

        % Validar que la pendiente es positiva y numéricamente significativa.
        if ~isfinite(slope_fit(1)) || abs(slope_fit(1)) <= eps
            continue
        end

        % --- Calcular el TTF por extrapolación de la recta ----------------
        % El "pie" de la onda es el tiempo en que la recta ajustada cruza Q=0:
        TTF_s(k) = -slope_fit(2) / slope_fit(1);

    end

    % --- Ajuste lineal TTF vs distancia → PWV ----------------------------
    [PWV_mps, fit_coeff, R2] = fit_pwv_linear(s_m, TTF_s);

end


%%%%%%%%%%%%%%%%%%%%%%%%% FUNCION CALCULO XCOR %%%%%%%%%%%%%%%%%%%%%%%%%
% estimate_pwv_xcor  Estima la PWV con correlación cruzada (XCor).
%
%   Entradas:
%     s_m       [n_planes x 1]           Posición de cada plano [m]
%     time_s    [n_phases x 1]           Eje temporal [s]
%     Q_mLs     [n_planes x n_phases]    Curvas de flujo [mL/s]
%     ref_plane escalar                  Índice del plano de referencia
%                                         (normalmente el más proximal = 1)
%
%   Salidas:
%     PWV_mps    escalar          PWV estimada [m/s]
%     Delay_s    [n_planes x 1]  Desfase de cada plano respecto a ref [s]
%     fit_coeff  [1 x 2]          Coeficientes del ajuste lineal
%     R2         escalar          R2 del ajuste

function [PWV_mps, Delay_s, fit_coeff, R2] = estimate_pwv_xcor(s_m, time_s, Q_mLs, ref_plane)

    n_planes = size(Q_mLs, 1);
    n_phases = size(Q_mLs, 2);
    Delay_s  = nan(n_planes, 1); % desfase temporal de cada plano [s]
    dt_s     = mean(diff(time_s)); % paso de tiempo medio [s]

    % --- Preparar señal de referencia -------------------------------------
    % Extraer la curva de flujo del plano de referencia (el más proximal)
    ref_flow = Q_mLs(ref_plane,:)';

    % Si la referencia no tiene ningún dato válido, la correlación cruzada
    % no puede calcularse para ningún plano → retornar NaN en todo
    if all(~isfinite(ref_flow))
        PWV_mps   = NaN;
        fit_coeff = [NaN NaN];
        R2        = NaN;
        return
    end

    % Rellenar NaN en la referencia por interpolación lineal
    ref_flow = fillmissing(ref_flow, 'linear', 'EndValues','nearest');

    % Restar la media para hacer la correlación insensible a offsets DC.
    % La correlación cruzada normalizada mide similitud de FORMA, no de nivel.
    % Si no se resta la media, una señal con offset alto domina el numerador
    % independientemente del desfase temporal.
    ref_flow = ref_flow - mean(ref_flow, 'omitnan');

    % --- Definir el rango de lags a explorar ------------------------------
    max_lag_s  = 0.30;

    % Convertir el límite temporal a número de muestras enteras.
    max_lag    = min(floor(n_phases/2), max(1, round(max_lag_s/dt_s)));

    % Vector de lags enteros a evaluar: desde -max_lag hasta +max_lag.
    lag_vector = (-max_lag:max_lag)';

    % --- Bucle sobre planos -----------------------------------------------
    for k = 1:n_planes

        % Extraer curva de flujo del plano k
        flow_curve = Q_mLs(k,:)';

        % Descartar planos sin datos válidos
        if all(~isfinite(flow_curve))
            continue
        end

        % Rellenar NaN y restar la media
        flow_curve = fillmissing(flow_curve, 'linear', 'EndValues','nearest');
        flow_curve = flow_curve - mean(flow_curve, 'omitnan');

        % Preasignar vector de valores de correlación para cada lag evaluado
        corr_vals = nan(numel(lag_vector), 1);

        % --- Bucle interno: correlación cruzada para cada lag discreto ----
        for ii = 1:numel(lag_vector)

            lag = lag_vector(ii);

            % Alinear los segmentos según el lag temporal:
            if lag >= 0
                seg_ref  = ref_flow(1:end-lag);
                seg_flow = flow_curve(1+lag:end);
            else
                seg_ref  = ref_flow(1-lag:end);
                seg_flow = flow_curve(1:end+lag);
            end

            % Máscara de posiciones donde ambas señales son finitas
            valid_pts = isfinite(seg_ref) & isfinite(seg_flow);

            % Se necesitan al menos 3 puntos para que la correlación sea
            % estadísticamente significativa
            if nnz(valid_pts) < 3
                continue
            end

            % Aplicar la máscara para trabajar solo con puntos válidos
            seg_ref  = seg_ref(valid_pts);
            seg_flow = seg_flow(valid_pts);

            % --- Correlación cruzada normalizada  -----------
            corr_denom = sqrt(sum(seg_ref.^2) * sum(seg_flow.^2));

            if corr_denom > 0
                corr_vals(ii) = sum(seg_ref .* seg_flow) / corr_denom;
            end

        end

        % Si no se obtuvo ningún valor de correlación válido, saltar este plano
        if all(~isfinite(corr_vals))
            continue
        end

        % Encontrar el lag entero con correlación máxima
        [~, idx_max_corr] = max(corr_vals);
        lag_refined = lag_vector(idx_max_corr); % lag óptimo entero (muestras)

        % --- Refinamiento parabólico del lag óptimo ----------------------
        % Al igual que en TTP, el verdadero máximo de la función de
        % correlación puede caer entre muestras enteras de lag.
        % Ajustar una parábola a los tres valores de correlación alrededor
        % del máximo discreto permite obtener el lag sub-muestra, con
        % resolución temporal mejor que dt_s.
        if idx_max_corr > 1 && idx_max_corr < numel(corr_vals)

            c_left   = corr_vals(idx_max_corr-1); % correlación un lag antes
            c_center = corr_vals(idx_max_corr);   % correlación en el lag óptimo
            c_right  = corr_vals(idx_max_corr+1); % correlación un lag después

            % Denominador de la fórmula parabólica del vértice
            parab_denom = c_left - 2*c_center + c_right;

            if isfinite(parab_denom) && abs(parab_denom) > eps
                % Desplazamiento fraccionario (en muestras) respecto al lag máximo
                delta = 0.5*(c_left - c_right) / parab_denom;

                % Limitar a ±1 muestra para no saltar fuera del entorno local
                delta = max(min(delta,1), -1);

                % Lag refinado (puede ser no entero, en unidades de muestras)
                lag_refined = lag_refined + delta;
            end

        end

        % Convertir el lag fraccionario de muestras a tiempo [s]
        Delay_s(k) = lag_refined * dt_s;

    end

    % --- Referenciar delays al plano proximal ----------------------------
    Delay_s = Delay_s - Delay_s(ref_plane);
    Delay_s(ref_plane) = 0; % garantizar exactamente cero en el plano ref

    % --- Ajuste lineal sobre distancia relativa al plano ref → PWV -------
    % Se trabaja con distancia relativa al plano de referencia para que el
    % ajuste pase por el origen (delay = 0 en distancia = 0), lo que es
    % físicamente correcto: la onda llega en t=0 a la posición de referencia.
    s_relative = s_m(:) - s_m(ref_plane);

    [PWV_mps, fit_coeff, R2] = fit_pwv_linear(s_relative, Delay_s);

end


function [PWV_mps, fit_coeff, R2] = fit_pwv_linear(s_m, t_s)

    valid_pts = isfinite(s_m) & isfinite(t_s);

    if nnz(valid_pts) < 3
        PWV_mps   = NaN;
        fit_coeff = [NaN NaN];
        R2        = NaN;
        return
    end

    dist_valid = s_m(valid_pts);
    time_valid = t_s(valid_pts);

    % Ajuste lineal: tiempo = pendiente * distancia + intercepto
    fit_coeff = polyfit(dist_valid, time_valid, 1);
    slope     = fit_coeff(1);

    if ~isfinite(slope) || abs(slope) <= eps
        PWV_mps = NaN;
    else
        PWV_mps = 1 / slope;

        if PWV_mps < 0
            warning('PWV calculada negativa (%.6f m/s). Revisar orientación del centerline.', PWV_mps)
        end
    end

    time_fitted = polyval(fit_coeff, dist_valid);
    SS_res      = sum((time_valid - time_fitted).^2);
    SS_tot      = sum((time_valid - mean(time_valid)).^2);

    if SS_tot > 0
        R2 = 1 - SS_res/SS_tot;
    else
        R2 = NaN;
    end

end
