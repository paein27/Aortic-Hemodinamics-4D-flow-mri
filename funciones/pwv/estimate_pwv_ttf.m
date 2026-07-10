function [PWV_mps, TTF_s, fit_coefficients, R2] = estimate_pwv_ttf( ...
    section_distance_m, time_vector_seconds, flow_rate_mL_per_s)

    % Time To Foot (foot-to-foot tipo Markl). Para cada corte:
    %   1) remuestrea la curva a grilla fina (pchip);
    %   2) estima un baseline diastolico (minimo pre-sistolico);
    %   3) ajusta una recta al tramo 20-80% de la AMPLITUD sobre baseline;
    %   4) el pie es la interseccion de esa recta con el baseline, no con
    %      flujo cero: extrapolar a cero es una extrapolacion larga que
    %      amplifica el error de pendiente y sesga el pie.
    % Sobre grilla fina el tramo 20-80% tiene decenas de puntos y la
    % pendiente es estable; sobre los 25 puntos crudos tiene 1-2 muestras y
    % el pie sale aleatorio (pendiente negativa, R2 ~ 0).

    number_of_sections = size(flow_rate_mL_per_s, 1);
    number_of_fine_samples = 1000;
    TTF_s = nan(number_of_sections, 1);

    for section_index = 1:number_of_sections
        [time_fine_s, flow_fine_mL_per_s, valid_curve] = ...
            resample_flow_curve_fine(time_vector_seconds, ...
            flow_rate_mL_per_s(section_index,:)', number_of_fine_samples);

        if ~valid_curve
            continue
        end

        [peak_flow_mL_per_s, peak_fine_index] = max(flow_fine_mL_per_s);

        if peak_fine_index < 3
            continue
        end

        baseline_flow_mL_per_s = min(flow_fine_mL_per_s(1:peak_fine_index));
        systolic_amplitude_mL_per_s = ...
            peak_flow_mL_per_s - baseline_flow_mL_per_s;

        if ~isfinite(systolic_amplitude_mL_per_s) || ...
                systolic_amplitude_mL_per_s <= 0
            continue
        end

        threshold_20 = ...
            baseline_flow_mL_per_s + 0.20*systolic_amplitude_mL_per_s;
        threshold_80 = ...
            baseline_flow_mL_per_s + 0.80*systolic_amplitude_mL_per_s;
        upslope_curve_mL_per_s = flow_fine_mL_per_s(1:peak_fine_index);
        index_20 = find(upslope_curve_mL_per_s >= threshold_20, 1, 'first');
        index_80 = find(upslope_curve_mL_per_s >= threshold_80, 1, 'first');

        if isempty(index_20) || isempty(index_80) || index_80 <= index_20
            continue
        end

        % Un ascenso limpio cruza el umbral del 20% una sola vez. Multiples
        % cruces indican una subida oscilante (jet/reflexiones tras una
        % coartacion severa): el pie no es fiable y se descarta el plano.
        number_of_20_percent_crossings = nnz( ...
            diff(upslope_curve_mL_per_s >= threshold_20) == 1);

        if number_of_20_percent_crossings > 1
            continue
        end

        upslope_sample_indices = (index_20:index_80)';

        if numel(upslope_sample_indices) < 3
            continue
        end

        upslope_fit_coefficients = polyfit( ...
            time_fine_s(upslope_sample_indices), ...
            flow_fine_mL_per_s(upslope_sample_indices), ...
            1);
        upslope_slope = upslope_fit_coefficients(1);
        upslope_intercept = upslope_fit_coefficients(2);

        if ~isfinite(upslope_slope) || upslope_slope <= eps
            continue
        end

        % El tramo 20-80% de un ascenso sistolico limpio es casi lineal. Un R2
        % local bajo indica una onda distorsionada (planos en/tras la
        % coartacion): el pie no es confiable y se descarta el plano.
        fitted_upslope_mL_per_s = polyval( ...
            upslope_fit_coefficients, time_fine_s(upslope_sample_indices));
        upslope_residual_ss = sum(( ...
            flow_fine_mL_per_s(upslope_sample_indices) - ...
            fitted_upslope_mL_per_s).^2);
        upslope_total_ss = sum(( ...
            flow_fine_mL_per_s(upslope_sample_indices) - mean( ...
            flow_fine_mL_per_s(upslope_sample_indices))).^2);
        upslope_local_r2 = 1 - upslope_residual_ss / max(upslope_total_ss, eps);

        if ~isfinite(upslope_local_r2) || upslope_local_r2 < 0.95
            continue
        end

        % Pie = interseccion de la recta de subida con el baseline diastolico.
        foot_time_s = ...
            (baseline_flow_mL_per_s - upslope_intercept) / upslope_slope;

        if foot_time_s < time_fine_s(1) || ...
                foot_time_s > time_fine_s(peak_fine_index)
            continue
        end

        TTF_s(section_index) = foot_time_s;
    end

    [PWV_mps, fit_coefficients, R2] = fit_pwv_linear( ...
        section_distance_m, TTF_s);
end
