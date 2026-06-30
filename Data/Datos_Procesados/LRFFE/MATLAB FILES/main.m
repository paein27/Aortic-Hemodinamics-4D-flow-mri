clc
clear
close all 

path(path,'../../../../Proyecto/4D-Flow-Matlab-Toolbox-main/iso2mesh/iso2mesh')

load 'FE Mesh'\elem.mat
load 'FE Mesh'\faces.mat
load 'FE Mesh'\nodes.mat
load 'FE Laplace'\Laplace.mat
load 'FE Velocity'\VEL.mat
load 'FE Area'\Area.mat
load 'FE Centerline'\Centerline.mat

figure, 
% ----------------- NUEVO -------------------------------
patch('faces',faces,'Vertices',nodes,...
      'edgecolor','k',...
      'facecolor','interp',...
      'FaceVertexCData',Laplace,...
      'FaceAlpha',0.3)

hold on
plot3(Centerline(:,1), Centerline(:,2), Centerline(:,3), ...
    'r-', 'LineWidth', 3)
hold off
% ----------------- NUEVO -------------------------------
daspect([1,1,1])
axis off
colorbar()
view(0,-90)

figure,
% ----------------- NUEVO -------------------------------
patch('faces',faces,'Vertices',nodes,...
      'edgecolor','k',...
      'facecolor','interp',...
      'FaceVertexCData',Laplace,...
      'FaceAlpha',0.3)

hold on
plot3(Centerline(:,1), Centerline(:,2), Centerline(:,3), ...
    'r-', 'LineWidth', 3)
hold off
% ----------------- NUEVO -------------------------------
daspect([1,1,1])
axis off
colorbar()
view(0,-90)


slices = 100;
x = linspace(min(Laplace(Laplace>0)), max(Laplace(Laplace<1)),slices);
area_vector = zeros(1,length(x));
velocity = zeros(length(x), 20);
mag_vel = squeeze(sqrt(sum(VEL.^2,2)));

figure, 
patch('faces',faces,'Vertices',nodes,'edgecolor','k','facecolor','r','Facealpha',0.1)
hold on 

for n = 1:slices
    [cutpos,~,facedata,elemid] = qmeshcut(elem,nodes,Laplace,x(n));
    selected_elements = elem(elemid,:);
    id_nodes = unique(selected_elements(:));
    velocity(n,:) = mean(mag_vel(id_nodes,:),1)*100;

    [cutpos,~,~,elemid] = qmeshcut(faces,nodes,Laplace,x(n));
    selected_faces = faces(elemid,:);
    id_nodes = unique(selected_faces(:));

    area_vector(n) = mean(Area(id_nodes));
    

    % selected_nodes = unique(selected_elements(:));
    plot3(nodes(id_nodes,1),nodes(id_nodes,2),nodes(id_nodes,3),'*y')
    patch('faces',selected_elements,'Vertices',nodes,'edgecolor','k','facecolor','b')
end
hold off
daspect([1,1,1])
axis off
colorbar()
view(0,-90)


flujo = velocity.*repmat(area_vector',1,20);

figure,
plot(flujo')
xlabel('cardiac_phase')
ylabel('flujo cm3/s')
grid on
