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
