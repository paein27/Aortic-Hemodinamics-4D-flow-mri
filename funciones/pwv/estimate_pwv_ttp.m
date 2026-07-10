function [PWV_mps, TTP_s, fit_coefficients, R2] = estimate_pwv_ttp( ...
    section_distance_m, time_vector_seconds, flow_rate_mL_per_s)

    % Time To Peak con localizacion sub-fase del instante del maximo. La curva
    % se remuestrea a grilla fina (pchip) y el peak se refina ajustando una
    % parabola a una ventana alrededor del maximo; se toma el vertice.
    %
    % Por que NO basta el argmax de la curva fina: pchip es monotona a trozos y
    % no hace overshoot, asi que su maximo cae siempre sobre una muestra
    % original -> TTP quedaria cuantizada a multiplos de dt (33 ms) y el grafico
    % muestra escalones. La parabola recupera resolucion continua entre fases.

    number_of_sections = size(flow_rate_mL_per_s, 1);
    number_of_fine_samples = 1000;
    mean_time_step_seconds = mean(diff(time_vector_seconds));
    TTP_s = nan(number_of_sections, 1);

    for section_index = 1:number_of_sections
        [time_fine_s, flow_fine_mL_per_s, valid_curve] = ...
            resample_flow_curve_fine(time_vector_seconds, ...
            flow_rate_mL_per_s(section_index,:)', number_of_fine_samples);

        if ~valid_curve
            continue
        end

        [~, peak_fine_index] = max(flow_fine_mL_per_s);

        % Un peak pegado a un borde no es un maximo sistolico bien definido.
        if peak_fine_index <= 1 || peak_fine_index >= number_of_fine_samples
            continue
        end

        TTP_s(section_index) = refine_peak_time_parabolic( ...
            time_fine_s, flow_fine_mL_per_s, peak_fine_index, ...
            mean_time_step_seconds);
    end

    [PWV_mps, fit_coefficients, R2] = fit_pwv_linear( ...
        section_distance_m, TTP_s);
end
