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
