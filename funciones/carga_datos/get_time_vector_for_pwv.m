function [time_vector_seconds, time_vector_source] = get_time_vector_for_pwv( ...
    time_file, number_of_time_phases)

    % Busca dentro de Time.mat la primera variable numerica compatible con
    % el numero de fases de VEL. Si el eje temporal parece estar en ms, lo
    % convierte a segundos y lo desplaza para iniciar en t = 0.

    time_vector_seconds = [];
    time_vector_source = '';

    if exist(time_file, 'file') == 2
        time_data = load(time_file);
        variable_names = fieldnames(time_data);

        for variable_index = 1:numel(variable_names)
            candidate_time_vector = time_data.(variable_names{variable_index});

            if isnumeric(candidate_time_vector) && ~isempty(candidate_time_vector)
                candidate_time_vector = squeeze(candidate_time_vector);

                if isvector(candidate_time_vector) && ...
                        numel(candidate_time_vector) >= number_of_time_phases

                    time_vector_seconds = double(candidate_time_vector(:));
                    time_vector_seconds = ...
                        time_vector_seconds(1:number_of_time_phases);
                    time_vector_source = ...
                        ['Time.mat -> variable: ', variable_names{variable_index}];
                    break
                end
            end
        end
    end

    if isempty(time_vector_seconds)
        error('No se encontro un eje temporal valido en %s.', time_file)
    end

    if max(abs(time_vector_seconds)) > 10
        time_vector_seconds = time_vector_seconds / 1000;
    end

    time_vector_seconds = time_vector_seconds - time_vector_seconds(1);

    if numel(time_vector_seconds) ~= number_of_time_phases
        error( ...
            'El eje temporal tiene %d puntos, pero VEL tiene %d fases.', ...
            numel(time_vector_seconds), number_of_time_phases)
    end

    if any(~isfinite(time_vector_seconds)) || ...
            any(diff(time_vector_seconds) <= 0)
        error('El eje temporal debe ser finito y estrictamente creciente.')
    end
end
