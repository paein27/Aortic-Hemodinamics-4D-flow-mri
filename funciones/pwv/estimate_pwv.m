function pwv = estimate_pwv(geometry, data, flow)
%ESTIMATE_PWV Calcula la PWV por los cuatro metodos del pipeline.
%   Coordina TTP, TTF, cross-correlation global (Markl) y cross-correlation
%   del par anatomico I-IV (Urbina/Sotelo). La distancia usada en el ajuste
%   lineal es la posicion acumulada sobre la centerline, en metros; PWV es
%   1/pendiente del ajuste tiempo vs distancia. Devuelve la estructura 'pwv'
%   con los resultados de cada metodo.

    section_arc_positions_m = geometry.section_arc_positions_m;
    urbina_section_arc_positions_m = geometry.urbina_section_arc_positions_m;
    section_I_urbina = geometry.section_I_urbina;
    section_IV_urbina = geometry.section_IV_urbina;

    time_vector_seconds = data.time_vector_seconds;
    flow_rate_mL_per_s = flow.flow_rate_mL_per_s;

    [PWV_TTP_mps, TTP_s, fit_TTP, R2_TTP] = estimate_pwv_ttp( ...
        section_arc_positions_m(:), time_vector_seconds, flow_rate_mL_per_s);

    [PWV_TTF_mps, TTF_s, fit_TTF, R2_TTF] = estimate_pwv_ttf( ...
        section_arc_positions_m(:), time_vector_seconds, flow_rate_mL_per_s);

    xcor_reference_section_index = section_I_urbina;

    [PWV_XCor_mps, Delay_XCor_s, fit_XCor, R2_XCor] = ...
        estimate_pwv_xcor( ...
            section_arc_positions_m(:), ...
            time_vector_seconds, ...
            flow_rate_mL_per_s, ...
            xcor_reference_section_index);

    [PWV_Urbina_XCor_mps, Delay_Urbina_XCor_s, fit_Urbina_XCor, ...
        R2_Urbina_XCor, urbina_pair_delay_s, urbina_pair_distance_m, ...
        urbina_pair_maximum_correlation] = ...
        estimate_pwv_xcor_urbina_pair( ...
            urbina_section_arc_positions_m(:), ...
            time_vector_seconds, ...
            flow_rate_mL_per_s, ...
            section_I_urbina, ...
            section_IV_urbina);

    pwv = struct();
    pwv.xcor_reference_section_index = xcor_reference_section_index;

    pwv.PWV_TTP_mps = PWV_TTP_mps;
    pwv.TTP_s = TTP_s;
    pwv.fit_TTP = fit_TTP;
    pwv.R2_TTP = R2_TTP;

    pwv.PWV_TTF_mps = PWV_TTF_mps;
    pwv.TTF_s = TTF_s;
    pwv.fit_TTF = fit_TTF;
    pwv.R2_TTF = R2_TTF;

    pwv.PWV_XCor_mps = PWV_XCor_mps;
    pwv.Delay_XCor_s = Delay_XCor_s;
    pwv.fit_XCor = fit_XCor;
    pwv.R2_XCor = R2_XCor;

    pwv.PWV_Urbina_XCor_mps = PWV_Urbina_XCor_mps;
    pwv.Delay_Urbina_XCor_s = Delay_Urbina_XCor_s;
    pwv.fit_Urbina_XCor = fit_Urbina_XCor;
    pwv.R2_Urbina_XCor = R2_Urbina_XCor;
    pwv.urbina_pair_delay_s = urbina_pair_delay_s;
    pwv.urbina_pair_distance_m = urbina_pair_distance_m;
    pwv.urbina_pair_maximum_correlation = urbina_pair_maximum_correlation;
end
