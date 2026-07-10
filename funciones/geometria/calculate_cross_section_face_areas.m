function face_area_mm2 = calculate_cross_section_face_areas( ...
    vertices_mm, face_node_indices)

    % qmeshcut puede devolver caras triangulares o cuadrilateras. Para
    % homogeneizar, cada cara se divide en dos triangulos: (1,2,3) y
    % (1,3,4). Si la cara es triangular, qmeshcut repite el indice 3 en la
    % columna 4, por lo que el segundo triangulo tiene area cero.

    first_vertex_mm = vertices_mm(face_node_indices(:,1),:);
    second_vertex_mm = vertices_mm(face_node_indices(:,2),:);
    third_vertex_mm = vertices_mm(face_node_indices(:,3),:);
    fourth_vertex_mm = vertices_mm(face_node_indices(:,4),:);

    first_triangle_area_mm2 = 0.5 * vecnorm( ...
        cross(second_vertex_mm - first_vertex_mm, ...
        third_vertex_mm - first_vertex_mm, 2), 2, 2);
    second_triangle_area_mm2 = 0.5 * vecnorm( ...
        cross(third_vertex_mm - first_vertex_mm, ...
        fourth_vertex_mm - first_vertex_mm, 2), 2, 2);

    face_area_mm2 = first_triangle_area_mm2 + second_triangle_area_mm2;
end
