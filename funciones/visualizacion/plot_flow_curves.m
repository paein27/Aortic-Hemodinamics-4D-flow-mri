function plot_flow_curves(data, geometry, flow)
%PLOT_FLOW_CURVES Curvas de flujo por corte en el tiempo.
%   Grafica el flujo (velocidad normal media por area transversal) de cada
%   corte a lo largo del ciclo.

    time_vector_seconds = data.time_vector_seconds;
    dataset_label = data.dataset_label;
    number_of_aortic_sections = geometry.number_of_aortic_sections;
    flow_rate_mL_per_s = flow.flow_rate_mL_per_s;

    figure('Name', char(dataset_label + ": curvas de flujo por corte"))
    plot(time_vector_seconds, flow_rate_mL_per_s', 'LineWidth', 1.2)
    xlabel('Tiempo [s]')
    ylabel('Flujo [mL/s]')
    title('Flujo = velocidad normal media por area transversal')
    legend(compose('Corte %d', 1:number_of_aortic_sections), ...
        'Location', 'eastoutside')
    grid on
end
