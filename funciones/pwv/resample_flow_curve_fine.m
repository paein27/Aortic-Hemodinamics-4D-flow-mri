function [time_fine_s, flow_fine_mL_per_s, valid_curve] = ...
    resample_flow_curve_fine( ...
        time_vector_seconds, flow_curve_mL_per_s, number_of_fine_samples)

    % Remuestrea una curva de flujo a una grilla temporal fina usando pchip
    % (interpolacion monotona a trozos, sin overshoot). A 33 ms la subida
    % sistolica tiene 1-2 muestras; sobre la grilla fina hay decenas de
    % puntos, lo que permite localizar peak y pie con resolucion sub-fase y
    % estabilizar el ajuste de la pendiente de subida.

    time_vector_seconds = time_vector_seconds(:);
    flow_curve_mL_per_s = flow_curve_mL_per_s(:);
    valid_curve = false;
    time_fine_s = linspace(time_vector_seconds(1), ...
        time_vector_seconds(end), number_of_fine_samples)';
    flow_fine_mL_per_s = nan(number_of_fine_samples, 1);

    if nnz(isfinite(flow_curve_mL_per_s)) < 4
        return
    end

    flow_curve_mL_per_s = fillmissing( ...
        flow_curve_mL_per_s, 'linear', 'EndValues', 'nearest');
    flow_fine_mL_per_s = interp1(time_vector_seconds, ...
        flow_curve_mL_per_s, time_fine_s, 'pchip');
    valid_curve = any(isfinite(flow_fine_mL_per_s)) && ...
        (max(flow_fine_mL_per_s) - min(flow_fine_mL_per_s)) > eps;
end
