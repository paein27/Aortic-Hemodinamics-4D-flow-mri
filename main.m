clc
close all
clear

%% Orquestador del pipeline PWV
% Este es el unico archivo que se edita para correr el analisis: selecciona el
% dataset y fija los parametros. Todo el procesamiento vive en la carpeta
% 'funciones', agrupado por tipo de utilidad en subcarpetas (carga_datos,
% geometria, flujo, pwv, visualizacion, reportes). genpath agrega la carpeta y
% todas sus subcarpetas al path.

project_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(project_root, 'funciones')))
addpath(fullfile(project_root, '4D-Flow-Matlab-Toolbox-main', 'iso2mesh'))

%% Configuracion
config = struct();

% Dataset a procesar. Cambiar de phantom = editar esta ruta.
config.data_folder = ['C:\Users\ignac\OneDrive\Desktop\USM\Memoria\Data\' ...
    'Datos_Procesados\9 mm rest\MATLAB FILES'];

% Numero de cortes transversales equidistantes sobre la centerline.
config.number_of_aortic_sections = 30;

% Posiciones anatomicas I y IV (Urbina/Sotelo) ancladas por longitud de arco
% FIJA desde la raiz de la centerline, no por indice de plano. Como es el
% mismo modelo de aorta en los 4 phantoms, la distancia fisica I-IV es fija;
% anclar por arco la mantiene consistente entre datasets. Antes I=3/IV=30
% (indices) hacia que la distancia variara ~25% con la extension de la
% segmentacion (centerline 196-246 mm), sesgando la PWV. Ambos ajustables:
%   - I ~ AAo, ejemplo 35 mm sobre la raiz.
%   - IV ~ DAo a nivel coronario; se deja <= la menor longitud de centerline
%     de los datasets (~196 mm) para que el plano exista en todos.
config.position_I_arc_from_root_mm = 35;
config.position_IV_arc_from_root_mm = 185;

%% Pipeline: carga, geometria, cortes y flujo
data = load_pwv_dataset(config);
geometry = build_section_geometry(data, config);
sections = extract_aortic_sections(data, geometry);
[sections, flow] = project_section_velocities(data, geometry, sections);

%% Visualizacion y resumen de flujo
visualization_phase_index = ...
    visualize_sections_3d(data, geometry, sections, flow);
summarize_section_flow(geometry, flow, sections);
plot_flow_curves(data, geometry, flow);

%% Calculo y reporte de PWV
pwv = estimate_pwv(geometry, data, flow);
report_pwv_results(data, geometry, pwv, visualization_phase_index);
plot_pwv_results(geometry, pwv);
