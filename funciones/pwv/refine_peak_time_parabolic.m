function peak_time_s = refine_peak_time_parabolic( ...
    time_fine_s, flow_fine_mL_per_s, peak_fine_index, frame_step_seconds)

    % Ajusta una parabola por minimos cuadrados a una ventana de +/- ~1 fase
    % alrededor del maximo de la curva fina y devuelve el instante del vertice.
    % Si la parabola no es concava hacia abajo o el vertice cae fuera de la
    % ventana, retorna el tiempo del maximo muestral como respaldo.

    peak_time_s = time_fine_s(peak_fine_index);
    fine_step_seconds = time_fine_s(2) - time_fine_s(1);
    half_window_samples = max(2, round(frame_step_seconds / fine_step_seconds));
    window_start_index = max(1, peak_fine_index - half_window_samples);
    window_end_index = ...
        min(numel(time_fine_s), peak_fine_index + half_window_samples);
    window_indices = (window_start_index:window_end_index)';

    if numel(window_indices) < 3
        return
    end

    % Se centra el tiempo en el peak muestral para condicionar bien el ajuste.
    centered_time_s = time_fine_s(window_indices) - peak_time_s;
    parabola_coefficients = polyfit( ...
        centered_time_s, flow_fine_mL_per_s(window_indices), 2);
    quadratic_term = parabola_coefficients(1);
    linear_term = parabola_coefficients(2);

    if ~isfinite(quadratic_term) || quadratic_term >= 0
        return
    end

    vertex_offset_s = -linear_term / (2*quadratic_term);

    if abs(vertex_offset_s) > (centered_time_s(end) - centered_time_s(1))
        return
    end

    peak_time_s = peak_time_s + vertex_offset_s;
end
