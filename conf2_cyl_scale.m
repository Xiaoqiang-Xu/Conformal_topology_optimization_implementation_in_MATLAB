%% Conformally map input mesh to  rectangle
function [F,tuv,d1,d2]=conf2_cyl_scale(F, xyz, extra)
bd = compute_bd(F);
plot_mesh(F,xyz);
corner_id = [7700;3281;4020;9692]; %%random conformal mapping
u_1 =  euclidean_ricci_flow_rec(F,xyz,corner_id); %% using ricci flow calculate the uniformization metric

uv= euclidean_embed(F,u_1,2);
XY=[real(uv),imag(uv)];

%% transform the triangle meshed figure
theta=(atan((XY(corner_id(2),2)-XY(corner_id(1),2))/(XY(corner_id(2),1)-XY(corner_id(1),1))));
T=[cos(theta), sin(theta); -sin(theta),cos(theta)];

tuv=(T*XY')'; % transformed uv matrix

d1=tuv(corner_id(1),1)-0;
d2=tuv(corner_id(1),2)-0;

tuv=[abs(tuv(:,1)-d1),abs(tuv(:,2)-d2)]; % put to the (0,0) point


d1=abs(tuv(corner_id(2),1)-tuv(corner_id(1),1));
d2=abs(tuv(corner_id(3),2)-tuv(corner_id(1),2));

tuv=[tuv(:,1)./d2,tuv(:,2)./d2];
d1=d1/d2;
d2=1;

