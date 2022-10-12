%% Topology optimization using dimension reduction level set method (DR-LSM), Xiaoqiang Xu, 2022 
clc;
close all;
clear all;
addpath('C:\Users\xxqhi\OneDrive\initial folder\Documents\CMADO_codes\DR_LSM_Heat_Conduction\Conformal_topology_optimization_implementation_in_MATLAB');
addpath(genpath('C:\Users\xxqhi\OneDrive\initial folder\Documents\CMADO_codes\DR_LSM_Heat_Conduction\Conformal_topology_optimization_implementation_in_MATLAB\level_set'));
addpath(genpath('C:\Users\xxqhi\OneDrive\initial folder\Documents\CMADO_codes\DR_LSM_Heat_Conduction\Conformal_topology_optimization_implementation_in_MATLAB\geometry_processing_package_master'));
parent_dir_name = 'D:\XLSM\figs';
%% Conformal mapping
input_mesh = 'sphere_random_surface4_29.obj'; 
[F, xyz,extra] = read_obj(input_mesh); % in the geometry processing tool box
[F,tuv,d1,d2]=conf2_cyl_scale(F, xyz, extra);
figure 
plot_mesh(F,tuv);

%% Create the grid.
g2.dim = 2;
g2.min = [0; 0];
g2.dx =1/200; 
g2.max = [d1; d2];
g2.bdry = @addGhostExtrapolate;
g2 = processGrid(g2);

%% Integration parameters.
ItNum = 1;
ReinitInterval =5;
TotalItNum =300;
Vol_target=0.4;

%% Initial level set Phi
nohx=4;
nohy=5;
r=0.04;
hdx=g2.max(1)/(nohx); % distance between neighboring circles in x direction
hdy=g2.max(2)/(nohy); % distance between neighboring circles in y direction
 LSgrid.Phi=ones(113,201); % assign the initial values for Phi, 1 everywhere
for ix=1:1:nohx
    for iy=1:1:nohy
        xc=hdx/2+hdx*(ix-1);
        yc=hdy/2+hdy*(iy-1);
        LSgrid.Phi1=sqrt((g2.xs{1}-xc).^2 + (g2.xs{2}-yc).^2) - r;
         LSgrid.Phi=min(LSgrid.Phi,LSgrid.Phi1); 
    end
end
x=g2.xs{1,1}(:,1)'; % x coordinates of all the fixed grid points
y=g2.xs{2,1}(1,:);  % y coordinates of all the fixed grid points
Xmax=max(x);
Ymax=max(y);
LSgrid.Phi = min(LSgrid.Phi, g2.xs{2});
LSgrid.Phi = min(LSgrid.Phi, Ymax-g2.xs{2});
LSgrid.Phi = min(LSgrid.Phi, g2.xs{1});
LSgrid.Phi = min(LSgrid.Phi, Xmax-g2.xs{1});
figure
contourf(g2.xs{1},g2.xs{2},LSgrid.Phi,[0,0])
axis equal;
colormap winter
Lambda=1;mu_L=5;alpha=0.99;

%% Calculate riemannian metric and intoplate onto grid 
G=riemannian_metric(F,tuv,xyz);
G_grid=griddata(tuv(:,1),tuv(:,2),G,g2.xs{1},g2.xs{2},'nearest');

%% Main loop
TC1=0.01;
while ItNum < TotalItNum
ItNum
%Reinitialize phi
if ItNum == 1 || mod(ItNum, ReinitInterval) == 0
  [LSgrid.Phi, g2] = LevelSetReinit2D(LSgrid.Phi, g2);
end

% Calculate curvature
[curvature_u, gradMag ] = curvatureSecond(g2, LSgrid.Phi);

% Set the upper and lower bound of the curvature
curvature_u = min(curvature_u, max(1./ g2.dx)); 
curvature_u = max(curvature_u, -max(1./ g2.dx));
LSgrid.k = curvature_u;

%% Output the Phi function on 2D grid to the 2D FEA

tmpx = (g2.vs{1}).'; tmpy = (g2.vs{2}).'; % tmpx tmpy refers to a row vector containing x and y coordinate information respectively..

tmpdata = LSgrid.Phi;

fid = fopen('Phi.txt', 'wt');

fprintf(fid,'%% grid \n');

fprintf(fid, '%e \t', tmpx);

fprintf(fid,'\n');

fprintf(fid, '%e \t', tmpy);

fprintf(fid,'\n%% data \n');

fprintf(fid, '%e \t', tmpdata);

fclose(fid);

clear tmpx tmpy tmpdata;

%% Output the Conformal factor file from the 2D mapped mesh rather than the 2D grid
% Write conformal factor as .txt
% only needs to be run once since conformal factor does not change during the optimization

tmpx = tuv(:,1)'; tmpy = tuv(:,2)';
tmpdata = G;
tmpdata=[tmpx;tmpy;tmpdata];
fid = fopen('Conformal_factor2.txt', 'wt');
fprintf(fid,'%%x %%y %%conformal_factor\r\n','x','y','phi');
fprintf(fid,'%f %f %f\r\t',tmpdata);
fclose(fid);
clear tmpx tmpy tmpdata;


%% Update the Augmented Lagrangian parameters
if ItNum>1
Lambda=max(0,Lambda+1/mu_L*(FEAResult.Vol-FEAResult.Vol_ref*Vol_target))
% Lambda=Lambda+1/mu_L*(FEAResult.Vol-2*Vol_target)
mu_L=alpha*mu_L
end

%% Call COMSOL FEA solver 
% Calculate the displacement and strain energy desnsity
FEAResult = FEM_2D_sphere_random_surface5_3_5_6(g2);
VolRatio=FEAResult.Vol/FEAResult.Vol_ref;
disp(('Calculating thermal compliance'));
TC=FEAResult.TC;
if abs((TC-TC1))/TC1<1e-6
    break;
end
disp(('Calculating shape gradient and velocity field...'));

disp(' ');
b=FEAResult.b;
T=FEAResult.T;
K=FEAResult.K;
Tx=FEAResult.Tx;
Ty=FEAResult.Ty;

LSgrid.Vn=-K.*1./G_grid.*(Tx.*Tx+Ty.*Ty);

% avoid singularity
[m1, n1] = size(LSgrid.Vn);

for i = 1: m1
   for j = 1: n1
       if g2.xs{1}(i,j)<=0.020 && g2.xs{2}(i,j)<=0.58 && g2.xs{2}(i,j)>=0.41
          LSgrid.Vn(i,j) = 0;  
       end
   end
end

for i = 1: m1
   for j = 1: n1
       if g2.xs{1}(i,j)>=0.54 && g2.xs{2}(i,j)<=0.58 && g2.xs{2}(i,j)>=0.41
          LSgrid.Vn(i,j) = 0;  
       end
   end
end



for i = 1: m1
   for j = 1: n1
       if g2.xs{2}(i,j)<=0.030 && g2.xs{1}(i,j)<=0.33 && g2.xs{1}(i,j)>=0.22
          LSgrid.Vn(i,j) = 0;  
       end
   end
end

for i = 1: m1
   for j = 1: n1
       if g2.xs{2}(i,j)>=0.97 && g2.xs{1}(i,j)<=0.33 && g2.xs{1}(i,j)>=0.22
          LSgrid.Vn(i,j) = 0;  
       end
   end
end

% LSgrid.Vn =50*(VolRatio-Vol_target)-TE; % using lambda for volume

LSgrid.Vn =Lambda+1/mu_L*(FEAResult.Vol-FEAResult.Vol_ref*Vol_target)+LSgrid.Vn; 
Vn_vol=Lambda+1/mu_L*(FEAResult.Vol-FEAResult.Vol_ref*Vol_target)  
% using Augmented lagrangian multiplier

LSgrid.k = LSgrid.k ./ max(1./ g2.dx); % normalize the curvature as well as the Vn
LSgrid.Vn = LSgrid.Vn - 3* LSgrid.k;
LSgrid.Vn = LSgrid.Vn ./ max(abs(LSgrid.Vn(:))); % normalize, divided by the maximum absolute value 
LSgrid.Vn = reshape(LSgrid.Vn, size(g2.xs{1}));
[m1, n1] = size(LSgrid.Vn);


%non-design area setting
[m1, n1] = size(LSgrid.Vn);

for i = 1: m1
   for j = 1: n1
       if g2.xs{1}(i,j)<=0.020 && g2.xs{2}(i,j)<=0.58 && g2.xs{2}(i,j)>=0.41
          LSgrid.Vn(i,j) = 0;  
       end
   end
end

for i = 1: m1
   for j = 1: n1
       if g2.xs{1}(i,j)>=0.54 && g2.xs{2}(i,j)<=0.58 && g2.xs{2}(i,j)>=0.41
          LSgrid.Vn(i,j) = 0;  
       end
   end
end



for i = 1: m1
   for j = 1: n1
       if g2.xs{2}(i,j)<=0.030 && g2.xs{1}(i,j)<=0.33 && g2.xs{1}(i,j)>=0.22
          LSgrid.Vn(i,j) = 0;  
       end
   end
end

for i = 1: m1
   for j = 1: n1
       if g2.xs{2}(i,j)>=0.97 && g2.xs{1}(i,j)<=0.33 && g2.xs{1}(i,j)>=0.22
          LSgrid.Vn(i,j) = 0;  
       end
   end
end

%% Modify velocity field by conformal factor, multiply e^-\lambda
LSgrid.Vn=1./sqrt(G_grid).*LSgrid.Vn; 

%% Level-set solution
% if(nargin < 1)
% 
% accuracy = 'medium';
% 
% end

% Set up spatial approximation scheme.

schemeFunc = @termNormal;

schemeData.speed = LSgrid.Vn;

schemeData.grid = g2;

% Set up time approximation scheme.

integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.

schemeData.derivFunc = @upwindFirstENO2;

integratorFunc = @odeCFL2;

 

singleStep = 0;

if(singleStep)

integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');

end

%--------------------------------------------------------------------------

y0 = LSgrid.Phi(:);

% % How far to step?

% tSpan = [ tNow, min(tMax, tNow + tPlot) ];

% tSpan = [0 1];
tSpan = [0 0.01]; %%% 3?1  %%0.2
% Take a timestep.

% [ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...

% integratorOptions, schemeData);

[ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
          integratorOptions, schemeData);

% tNow = t(end);

% Get back the correctly shaped data array

data = reshape(y, g2.shape);

LSgrid.Phi = data;
  
if ItNum == 1 || mod(ItNum, 1) == 0

h = figure; 
set(h, 'visible','off');
T=contourf(g2.xs{1}, g2.xs{2}, LSgrid.Phi, [0 0]);
grid on
grid minor
axis equal;
filename=[num2str(ItNum)];
savefig();
FileName=[parent_dir_name,'\Fig_', num2str(ItNum),'.jpg'];
saveas(h,FileName);
delete(h);      
end

h3 = figure;
set(h3, 'visible','off');
tmp = 1:ItNum;
VR(ItNum) = VolRatio; % volume ratio
Obj(ItNum) = TC;
[AX,H1,H2] = plotyy(tmp, Obj, tmp, VR ,'plot');
set(get(AX(1),'Ylabel'),'String','Obj');
set(get(AX(2),'Ylabel'),'String','Volume Ratio (VR)');
title(['Obj(--) = ',num2str(Obj(end)),', VR(-) = ',num2str(VR(end))]);
set(H1,'LineStyle','--');
set(H2,'LineStyle','-');
grid on
FileName=[parent_dir_name,'\Obj.jpg'];
saveas(h3,FileName);
delete(h3);
ItNum = ItNum +1;
%% Convergence criterion
TC1=TC;
end



