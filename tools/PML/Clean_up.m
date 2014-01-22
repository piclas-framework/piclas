function Clean_up(matpath,copy,input)
%--------------------------------------------------------------------------
% Clean up
%--------------------------------------------------------------------------
disp('Cleaning up.....');
if copy==1
  unix('rm flexi');
  unix('rm postrec');
  unix('rm preproctool');
  unix('rm recordpoints');
  unix('rm visu3D');
end


if 1==2
unix(['rm -rf ' input.ProjectName '_State']); % State Files
unix(['rm ' input.ProjectName '_RPSet.h5']);  % Recordpoint setup
unix(['rm ' input.ProjectName '_mesh.h5']);   % Gitter
unix('rm *.plt');                             % TecPlot Files
unix(['rm -rf ' input.ProjectName '_RP']);    % Record Points


  unix('rm -rf visit_files');                   % Visu3D Files
  unix('rm -rf visit_output');                  % Bilder
  unix('rm *.out');                             % Output Files
  unix('rm -rf Probes');                        % RP Diagrams
  unix('rm Database.csv');                      % Energie
end

cd([matpath]);
disp('Done............');
end
