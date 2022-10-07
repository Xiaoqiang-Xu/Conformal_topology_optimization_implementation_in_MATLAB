function FEAResult = FEM_2D_sphere_random_surface5_3_5_6(g2)
%
% FEM_2D_sphere_random_surface5_3_5_6.m
%
% Model exported on May 7 2021, 02:37 by COMSOL 5.3.0.223.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\xxqhi\Dropbox\XLSM\Code_Xiaoqiang\R3_DR_LSM');

model.label('FEM_2D_sphere_random_surface5_3_5_6.mph');

model.comments(['Untitled\n\n']);

model.param.set('K0', '0.001');
model.param.set('K1', '10');
model.param.set('rho0', '0.001');
model.param.set('rho1', '1');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.file.create('res4');
model.file.create('res5');

model.result.table.create('tbl1', 'Table');

model.component('comp1').func.create('int1', 'Interpolation');
model.component('comp1').func.create('int2', 'Interpolation');
model.component('comp1').func('int1').set('sourcetype', 'model');
model.component('comp1').func('int1').set('modelres', 'res4');
model.component('comp1').func('int1').set('importedname', 'Phi.txt');
model.component('comp1').func('int1').set('importedstruct', 'Grid');
model.component('comp1').func('int1').set('importeddim', '2D');

model.file('res4').resource('C:\Users\xxqhi\OneDrive - Stony Brook University\XLSM\Code_Xiaoqiang\R3_DR_LSM\Phi.txt');

model.component('comp1').func('int1').set('source', 'file');
model.component('comp1').func('int1').set('nargs', '1');
model.component('comp1').func('int1').set('struct', 'grid');
model.component('comp1').func('int2').set('sourcetype', 'model');
model.component('comp1').func('int2').set('modelres', 'res5');
model.component('comp1').func('int2').set('importedname', 'Conformal_factor2.txt');
model.component('comp1').func('int2').set('importedstruct', 'Spreadsheet');
model.component('comp1').func('int2').set('importeddim', '2D');

model.file('res5').resource('C:\Users\xxqhi\OneDrive - Stony Brook University\XLSM\Code_Xiaoqiang\R3_DR_LSM\Conformal_factor2.txt');

model.component('comp1').func('int2').set('source', 'file');
model.component('comp1').func('int2').set('nargs', '2');
model.component('comp1').func('int2').set('struct', 'spreadsheet');

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('size', [0.5625 1]);
model.component('comp1').geom('geom1').create('pt1', 'Point');
model.component('comp1').geom('geom1').feature('pt1').setIndex('p', '0.256', 0, 0);
model.component('comp1').geom('geom1').feature('pt1').setIndex('p', '1', 1, 0);
model.component('comp1').geom('geom1').create('pt2', 'Point');
model.component('comp1').geom('geom1').feature('pt2').setIndex('p', '0.3048', 0, 0);
model.component('comp1').geom('geom1').feature('pt2').setIndex('p', '1', 1, 0);
model.component('comp1').geom('geom1').create('pt3', 'Point');
model.component('comp1').geom('geom1').feature('pt3').setIndex('p', '0.2579', 0, 0);
model.component('comp1').geom('geom1').feature('pt3').setIndex('p', '0', 1, 0);
model.component('comp1').geom('geom1').create('pt4', 'Point');
model.component('comp1').geom('geom1').feature('pt4').setIndex('p', '0.3066', 0, 0);
model.component('comp1').geom('geom1').feature('pt4').setIndex('p', '0', 1, 0);
model.component('comp1').geom('geom1').create('pt5', 'Point');
model.component('comp1').geom('geom1').feature('pt5').setIndex('p', '0', 0, 0);
model.component('comp1').geom('geom1').feature('pt5').setIndex('p', '0.4467', 1, 0);
model.component('comp1').geom('geom1').create('pt6', 'Point');
model.component('comp1').geom('geom1').feature('pt6').setIndex('p', '0', 0, 0);
model.component('comp1').geom('geom1').feature('pt6').setIndex('p', '0.5539', 1, 0);
model.component('comp1').geom('geom1').create('pt7', 'Point');
model.component('comp1').geom('geom1').feature('pt7').setIndex('p', '0.5625', 0, 0);
model.component('comp1').geom('geom1').feature('pt7').setIndex('p', '0.4457', 1, 0);
model.component('comp1').geom('geom1').create('pt8', 'Point');
model.component('comp1').geom('geom1').feature('pt8').setIndex('p', '0.5625', 0, 0);
model.component('comp1').geom('geom1').feature('pt8').setIndex('p', '0.5529', 1, 0);
model.component('comp1').geom('geom1').run;

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('hs', 'flc1hs(phi,0.01)', 'heaviside function');
model.component('comp1').variable('var1').set('phi', 'int1(x,y)', 'level set function');
model.component('comp1').variable('var1').set('K2', 'hs*K1+(1-hs)*K0', 'interpolated K');
model.component('comp1').variable('var1').set('rho2', 'hs*rho1+(1-hs)*rho0', 'interpolated rho');
model.component('comp1').variable('var1').set('Thermal_compliance', '(K2*(Tx*Tx+Ty*Ty))');
model.component('comp1').variable('var1').set('b', '10');
model.component('comp1').variable('var1').set('conformal_factor', 'int2(x,y)');

model.component('comp1').material.create('mat1', 'Common');

model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl('intop1').selection.set([1]);

model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');
model.component('comp1').physics('ht').create('temp1', 'TemperatureBoundary', 1);
model.component('comp1').physics('ht').feature('temp1').selection.set([6 7]);
model.component('comp1').physics('ht').create('lihs1', 'LineHeatSource', 0);
model.component('comp1').physics('ht').feature('lihs1').selection.set([5]);
model.component('comp1').physics('ht').create('hf1', 'HeatFluxBoundary', 1);
model.component('comp1').physics('ht').feature('hf1').selection.set([3 11]);

model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');

model.component('comp1').probe.create('var1', 'GlobalVariable');

model.result.table('tbl1').label('Probe Table 1');

model.component('comp1').view('view1').axis.set('xmin', -0.41373562812805176);
model.component('comp1').view('view1').axis.set('xmax', 0.9762356281280518);
model.component('comp1').view('view1').axis.set('ymin', -0.048424094915390015);
model.component('comp1').view('view1').axis.set('ymax', 1.0484241247177124);
model.component('comp1').view('view1').axis.set('abstractviewlratio', -0.7355300188064575);
model.component('comp1').view('view1').axis.set('abstractviewrratio', 0.7355300188064575);
model.component('comp1').view('view1').axis.set('abstractviewbratio', -0.04999998211860657);
model.component('comp1').view('view1').axis.set('abstractviewtratio', 0.04999995231628418);
model.component('comp1').view('view1').axis.set('abstractviewxscale', 0.0015759310917928815);
model.component('comp1').view('view1').axis.set('abstractviewyscale', 0.0015759310917928815);

model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'K2' '0' '0' '0' 'K2' '0' '0' '0' 'K2'});
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', '1');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho2');

model.component('comp1').physics('ht').prop('PhysicalModelProperty').set('dz', '0.01[m]');
model.component('comp1').physics('ht').feature('solid1').set('k', {'K2'; '0'; '0'; '0'; 'K2'; '0'; '0'; '0'; 'K2'});
model.component('comp1').physics('ht').feature('solid1').set('rho', 'rho2');
model.component('comp1').physics('ht').feature('temp1').set('T0', 0);
model.component('comp1').physics('ht').feature('lihs1').set('Ql', 10);
model.component('comp1').physics('ht').feature('lihs1').set('Pl', 10);
model.component('comp1').physics('ht').feature('lihs1').active(false);
model.component('comp1').physics('ht').feature('hf1').set('q0', 20.9328);

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('hmax', 0.0077);
model.component('comp1').mesh('mesh1').run;

model.component('comp1').probe('var1').set('expr', 'intop1((K2*(Tx*Tx+Ty*Ty)))');
model.component('comp1').probe('var1').set('unit', '');
model.component('comp1').probe('var1').set('descractive', true);
model.component('comp1').probe('var1').set('descr', 'intop1((K2*(Tx*Tx+Ty*Ty)))');
model.component('comp1').probe('var1').set('table', 'tbl1');
model.component('comp1').probe('var1').set('window', 'window1');

model.component('comp1').physics('ht').feature('solid1').set('k_mat', 'userdef');
model.component('comp1').physics('ht').feature('solid1').set('rho_mat', 'userdef');

model.study.create('std1');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('d1', 'Direct');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').create('i2', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i2').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.result.dataset.create('dset2', 'Solution');
model.result.dataset('dset2').set('probetag', 'var1');
model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical('gev1').set('data', 'dset2');
model.result.numerical('gev1').set('probetag', 'var1');
model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup2D');
model.result.create('pg3', 'PlotGroup1D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg2').create('con1', 'Contour');
model.result('pg3').set('probetag', 'window1_default');
model.result('pg3').create('tblp1', 'Table');
model.result('pg3').feature('tblp1').set('probetag', 'var1');

model.component('comp1').probe('var1').genResult([]);

model.sol('sol1').attach('std1');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'd1');
model.sol('sol1').feature('s1').feature('fc1').set('initstep', 0.01);
model.sol('sol1').feature('s1').feature('fc1').set('minstep', 1.0E-6);
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 50);
model.sol('sol1').feature('s1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('i1').label('Algebraic Multigrid Solver (ht)');
model.sol('sol1').feature('s1').feature('i2').label('Geometric Multigrid Solver (ht)');
model.sol('sol1').runAll;

model.result.dataset('dset2').label('Probe Solution 2');
model.result.numerical('gev1').set('unit', {'K^2'});
model.result.numerical('gev1').setResult;
model.result('pg1').label('Temperature (ht)');
model.result('pg1').feature('surf1').label('Surface');
model.result('pg1').feature('surf1').set('colortable', 'ThermalLight');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg2').label('Isothermal Contours (ht)');
model.result('pg2').feature('con1').label('Contour');
model.result('pg2').feature('con1').set('colortable', 'ThermalLight');
model.result('pg2').feature('con1').set('resolution', 'normal');
model.result('pg3').set('xlabel', 'intop1((K2*(Tx*Tx+Ty*Ty)))');
model.result('pg3').set('ylabel', 'intop1((K2*(Tx*Tx+Ty*Ty)))');
model.result('pg3').set('xlabelactive', false);
model.result('pg3').set('ylabelactive', false);

xx = [g2.xs{1}(:), g2.xs{2}(:)];
xx = xx.';

FEAResult.T= mphinterp(model,'T','coord',xx);
FEAResult.K= mphinterp(model,'K2','coord',xx);
FEAResult.b= mphinterp(model,'b','coord',xx);
FEAResult.Tx= mphinterp(model,'Tx','coord',xx);
FEAResult.Ty= mphinterp(model,'Ty','coord',xx);


FEAResult.T = reshape(FEAResult.T,size(g2.xs{1}));
FEAResult.K = reshape(FEAResult.K,size(g2.xs{1}));
FEAResult.b = reshape(FEAResult.b,size(g2.xs{1}));
FEAResult.Tx = reshape(FEAResult.Tx,size(g2.xs{1}));
FEAResult.Ty = reshape(FEAResult.Ty,size(g2.xs{1}));

%% other information
FEAResult.Vol=mphint2(model,'hs*int2(x,y)','surface');
FEAResult.Vol_ref=mphint2(model,'1*int2(x,y)','surface');
FEAResult.TC=mphint2(model,'K2*(Tx*Tx+Ty*Ty)','surface');

end
