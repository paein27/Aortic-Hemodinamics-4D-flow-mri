function [PWV_mps, Delay_s, fit_coefficients, R2, ...
    pair_delay_s, pair_distance_m, maximum_correlation] = ...
    estimate_pwv_xcor_urbina_pair( ...
        section_distance_m, time_vector_seconds, flow_rate_mL_per_s, ...
        section_I, section_IV)

    number_of_sections = size(flow_rate_mL_per_s, 1);
    Delay_s = nan(number_of_sections, 1);
    fit_coefficients = [NaN NaN];
    R2 = NaN;
    PWV_mps = NaN;
    pair_delay_s = NaN;
    pair_distance_m = NaN;
    maximum_correlation = NaN;

    if section_I < 1 || section_I > number_of_sections || ...
            section_IV < 1 || section_IV > number_of_sections
        warning('Las secciones XCor I-IV estan fuera del rango de cortes.')
        return
    end

    flow_I = flow_rate_mL_per_s(section_I,:)';
    flow_IV = flow_rate_mL_per_s(section_IV,:)';

    if all(~isfinite(flow_I)) || all(~isfinite(flow_IV))
        warning('XCor I-IV no tiene curvas de flujo validas.')
        return
    end

    flow_I = fillmissing(flow_I, 'linear', 'EndValues', 'nearest');
    flow_IV = fillmissing(flow_IV, 'linear', 'EndValues', 'nearest');

    flow_I = flow_I - mean(flow_I, 'omitnan');
    flow_IV = flow_IV - mean(flow_IV, 'omitnan');

    mean_time_step_seconds = mean(diff(time_vector_seconds));
    maximum_lag_seconds = 0.30;
    maximum_lag_samples = min( ...
        floor(numel(time_vector_seconds)/2), ...
        max(1, round(maximum_lag_seconds/mean_time_step_seconds)));
    lag_sample_vector = (-maximum_lag_samples:maximum_lag_samples)';
    correlation_by_lag = nan(numel(lag_sample_vector), 1);

    for lag_index = 1:numel(lag_sample_vector)
        current_lag_samples = lag_sample_vector(lag_index);

        if current_lag_samples >= 0
            reference_segment = flow_I(1:end-current_lag_samples);
            section_segment = flow_IV(1+current_lag_samples:end);
        else
            reference_segment = flow_I(1-current_lag_samples:end);
            section_segment = flow_IV(1:end+current_lag_samples);
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
        warning('XCor I-IV no pudo calcular correlaciones validas.')
        return
    end

    [maximum_correlation, maximum_correlation_index] = ...
        max(correlation_by_lag);
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

    pair_delay_s = refined_lag_samples * mean_time_step_seconds;
    pair_distance_m = ...
        section_distance_m(section_IV) - section_distance_m(section_I);

    if ~isfinite(pair_delay_s) || pair_delay_s <= 0 || ...
            ~isfinite(pair_distance_m) || pair_distance_m <= 0
        warning('XCor I-IV entrego delay o distancia invalida.')
        return
    end

    PWV_mps = pair_distance_m / pair_delay_s;
    fit_coefficients = [pair_delay_s / pair_distance_m, ...
        -pair_delay_s / pair_distance_m * section_distance_m(section_I)];

    Delay_s(section_I) = 0;
    Delay_s(section_IV) = pair_delay_s;
end
