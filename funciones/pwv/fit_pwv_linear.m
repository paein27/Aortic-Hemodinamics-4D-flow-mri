function [PWV_mps, fit_coefficients, R2] = fit_pwv_linear( ...
    section_distance_m, characteristic_time_seconds)

    % Ajuste robusto de tiempo vs distancia. Si t = a*s + b, PWV = 1/a. Se
    % rechazan planos outlier (p.ej. pies mal detectados sobre ondas
    % distorsionadas por una coartacion severa) con robustfit; si no hay
    % Statistics Toolbox, se usa un descarte iterativo por MAD. Sin este
    % rechazo, unos pocos planos vuelcan la pendiente global (incluso a
    % negativa). Sobre datos limpios (sin outliers) equivale al OLS.

    PWV_mps = NaN;
    fit_coefficients = [NaN NaN];
    R2 = NaN;

    valid_sample_mask = ...
        isfinite(section_distance_m(:)) & ...
        isfinite(characteristic_time_seconds(:));

    if nnz(valid_sample_mask) < 3
        return
    end

    valid_section_distance_m = section_distance_m(valid_sample_mask);
    valid_section_distance_m = valid_section_distance_m(:);
    valid_characteristic_time_seconds = ...
        characteristic_time_seconds(valid_sample_mask);
    valid_characteristic_time_seconds = valid_characteristic_time_seconds(:);

    inlier_mask = true(size(valid_section_distance_m));

    if exist('robustfit', 'file') == 2
        [~, robust_stats] = robustfit( ...
            valid_section_distance_m, valid_characteristic_time_seconds);
        robust_residuals = robust_stats.resid;
        residual_scale = 1.4826*median(abs( ...
            robust_residuals - median(robust_residuals)));

        if ~isfinite(residual_scale) || residual_scale <= eps
            residual_scale = std(robust_residuals);
        end

        if isfinite(residual_scale) && residual_scale > eps
            inlier_mask = abs(robust_residuals) <= 3*residual_scale;
        end
    else
        for iteration = 1:5
            current_fit = polyfit( ...
                valid_section_distance_m(inlier_mask), ...
                valid_characteristic_time_seconds(inlier_mask), 1);
            residuals = valid_characteristic_time_seconds - ...
                polyval(current_fit, valid_section_distance_m);
            residual_scale = 1.4826*median(abs( ...
                residuals(inlier_mask) - median(residuals(inlier_mask))));

            if ~isfinite(residual_scale) || residual_scale <= eps
                break
            end

            updated_mask = abs(residuals) <= 3*residual_scale;

            if nnz(updated_mask) < 3 || isequal(updated_mask, inlier_mask)
                break
            end

            inlier_mask = updated_mask;
        end
    end

    if nnz(inlier_mask) < 3
        return
    end

    fit_coefficients = polyfit( ...
        valid_section_distance_m(inlier_mask), ...
        valid_characteristic_time_seconds(inlier_mask), 1);
    slope_seconds_per_meter = fit_coefficients(1);

    if ~isfinite(slope_seconds_per_meter) || ...
            abs(slope_seconds_per_meter) <= eps
        PWV_mps = NaN;
    elseif slope_seconds_per_meter <= 0
        PWV_mps = NaN;
        warning( ...
            ['Pendiente temporal no positiva (%.6f s/m) tras ajuste ', ...
            'robusto. PWV no fisiologica; se reporta NaN.'], ...
            slope_seconds_per_meter)
    else
        PWV_mps = 1/slope_seconds_per_meter;
    end

    fitted_time_seconds = polyval( ...
        fit_coefficients, valid_section_distance_m(inlier_mask));
    residual_sum_of_squares = sum(( ...
        valid_characteristic_time_seconds(inlier_mask) - ...
        fitted_time_seconds).^2);
    total_sum_of_squares = sum(( ...
        valid_characteristic_time_seconds(inlier_mask) - mean( ...
        valid_characteristic_time_seconds(inlier_mask))).^2);

    if total_sum_of_squares > 0
        R2 = 1 - residual_sum_of_squares/total_sum_of_squares;
    end
end
