function face_mean_value = calculate_face_mean_vertex_values( ...
    vertex_values, face_node_indices)

    % Promedia los valores nodales de cada cara. Las caras triangulares se
    % tratan con solo sus tres vertices reales para evitar contar dos veces
    % el vertice repetido por qmeshcut.

    values_per_face = vertex_values(face_node_indices);
    face_mean_value = mean(values_per_face, 2);

    triangular_face_mask = face_node_indices(:,3) == face_node_indices(:,4);
    face_mean_value(triangular_face_mask) = ...
        mean(values_per_face(triangular_face_mask,1:3), 2);
end
