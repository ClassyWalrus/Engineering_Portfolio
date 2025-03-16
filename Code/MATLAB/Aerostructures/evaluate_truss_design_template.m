function evaluate_truss_design_template(user_x500)
%
% example usage: evaluate_truss_design_template('smith0000')

warning('error','MATLAB:singularMatrix');
warning('error','MATLAB:nearlySingularMatrix');
warning('error','MATLAB:load:variableNotFound');

report_name = ['truss_', user_x500, '.txt'];
report = fopen(report_name, 'w');

fprintf(report, 'Truss design analysis report.\n');
fprintf(report, '\n');
fprintf(report, 'Report for: %s\n', user_x500);
fprintf(report, 'Generated at %s\n', datetime);
fprintf(report, '\n');

mat_file_name = ['truss_', user_x500, '.mat'];
try
  PDorig=load(mat_file_name,'PD');
catch
  fprintf(report, 'Unable to open file: %s.\n', mat_file_name);
  fprintf(report, '\n');
  fprintf(report, 'Truss design does *NOT* satisfy all requirements.\n');
  fclose(report);
  return
end

try
  PDorig=PDorig.PD;
catch
  fprintf(report, 'Unable to find variable "PD" in file: %s.\n', ...
                  mat_file_name);
  fprintf(report, '\n');
  fprintf(report, 'Truss design does *NOT* satisfy all requirements.\n');
  fclose(report);
  return
end

fail = false;  % assume success


fprintf(report, 'Loaded design file: %s\n', mat_file_name);
fprintf(report, '\n');

fprintf(report, 'Constraint 2.1.1:\n');
fprintf(report, '\tAssuming conformance\n');
fprintf(report, '\t** OK **\n');
fprintf(report, '\n');

try
  PD.N=PDorig.N;
  if(~isscalar(PD.N))
    fprintf(report, 'Wrong shape for "PD.N".\n');  fprintf(report, '\n');
    fprintf(report, 'Truss design does *NOT* satisfy all requirements.\n');
    fclose(report);
    return
  end
catch
  fprintf(report, 'Unable to find "PD.N".\n');  fprintf(report, '\n');
  fprintf(report, 'Truss design does *NOT* satisfy all requirements.\n');
  fclose(report);
  return
end

try
  PD.NE=PDorig.NE;
  if(~isscalar(PD.NE))
    fprintf(report, 'Wrong shape for "PD.NE".\n');  fprintf(report, '\n');
    fprintf(report, 'Truss design does *NOT* satisfy all requirements.\n');
    fclose(report);
    return
  end
catch
  fprintf(report, 'Unable to find "PD.NE".\n');
  fprintf(report, '\n');
  fprintf(report, 'Truss design does *NOT* satisfy all requirements.\n');
  fclose(report);
  return
end

fprintf(report, 'Constraint 2.1.2:\n');
fprintf(report, '\tNumber of Nodes: %i\n', PD.N);
fprintf(report, '\tNumber of Bars:  %i\n', PD.NE);
fprintf(report, '\t** OK **\n');
fprintf(report, '\n');

try
  PD.NodePos = PDorig.NodePos();
catch
  fprintf(report, 'Unable to find "PD.NodePos".\n');
  fprintf(report, '\n');
  fprintf(report, 'Truss design does *NOT* satisfy all requirements.\n');
  fclose(report);
  return
end
if ((size(PD.NodePos,1) ~= PD.N) | (size(PD.NodePos,2) ~= 3))
  fprintf(report, 'Wrong size for "PD.NodePos".\n');
  fprintf(report, '\n');
  fprintf(report, 'Truss design does *NOT* satisfy all requirements.\n');
  fclose(report);
  return
else
  PD.NodePos = PD.NodePos(1:PD.N,1:3);
end

%%%%%%%
%
% ... consider criteria that can be evaluated before computing solution
%
%
%%%%%%%
fprintf(report, '\n');

fprintf(report, 'Performing analysis of loaded/deformed structure:\n');

try
  PDans=PD_truss_static(PD);
catch
  fprintf(report, 'Unable to successfully execute "PD_truss_static".\n');
  fprintf(report, '\n');
  fail = true;
  fprintf(report, '\t** VIOLATION **\n');
  fprintf(report, 'Truss design does *NOT* satisfy all requirements.\n');
  fclose(report);
  return
end

if (sum(isnan(PDans.ElmStress)) > 0)
  fprintf(report, 'PD_truss_static computes "NaN" for some stresses.\n');
  fprintf(report, 'There is something wrong with the design.\n');
  fprintf(report, '\n');
  fprintf(report, '\t** VIOLATION **\n');
  fail = true;
end
fprintf(report, '\n');

%%%%%%%
%
% ... consider criteria that can only be evaluated after computing solution
%
%
%%%%%%%

truss_mass = 0.0;
for i=1:PD.NE
  bconn = PD.ElmConnect(i,:);
  blen = norm(PD.NodePos(bconn(2),:) - PD.NodePos(bconn(1),:));
  bmat = PD.MatsSets(PD.ElmMats(i));

  bar_mass = blen * bmat.A * bmat.rho;

  truss_mass = truss_mass + bar_mass;
end

fprintf(report, '---------------------------------------------------------\n');
fprintf(report, '\n');
fprintf(report, 'Truss mass: %e\n', truss_mass);
fprintf(report, '\n');

if (fail)
    fprintf(report, 'Truss design does *NOT* satisfy all requirements.\n');
else
    fprintf(report, 'Truss design *DOES* satisfy all requirements.\n');
end

fclose(report);

return;
end
