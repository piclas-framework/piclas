function [  ] = write_datafile( filename, X, Y )
%-------------------------------------------------------------------------%
% short function to write a 2d table file
%
% use with care
%
% Script written by P. Ortwein.
% Last modified:      2012-21-5 
%-------------------------------------------------------------------------%

format = {'%5.0f',' %0.5e'};

outid = fopen(filename, 'w+');
dlmwrite (filename,[X(:),Y(:)],'delimiter', '\t','roffset',0,'-append',...
            'precision', char(format(2)), 'newline','unix');
%-------------------------------------------------------------------------%
end

