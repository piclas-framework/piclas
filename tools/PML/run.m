function time=run(name,copypath,command)
disp('-------------------------------------------------------------------');
disp(['Running ' name]);
disp(['Running Command: ' copypath command])
tic
[status, result] = unix([copypath command]);
% if status == 127 % Problem mit MPI, MATLAB lädt nicht automatisch die Bibliotheken von MPI wie das die Shell tut
%   disp('Running Command: /home/stephen/Diplomarbeit/PMLalgorithm/./flexi_single parameter_flexi.ini 1>std.out 2>err.out')
%   [status, result] = unix(['/home/stephen/Diplomarbeit/PMLalgorithm/./flexi_single parameter_flexi.ini 1>std.out 2>err.out']);
% end
time=toc;
if status ~=0
  disp(['status = ' num2str(status)]);
  disp(['result = ' result]);
  error(['Error in ' name ': Calculation failed']);
end
disp(['End ' name ': ' num2str(time)]);
disp('-------------------------------------------------------------------');
end