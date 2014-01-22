function difference = calcDifferenceInterpolation(f1, f2)


% f1 = [x1 y1]    f2 = [x1 y1]
%       :  :            :  :
%      [xN yN]         [xK yK]

interpol_temp=zeros(length(f1),1); % interpolierte kurve auf f2

for J=1:length(f1)
  interpol_temp(J)=interpolate(f2,f1(J,1));
end

difference = [f1(:,1) interpol_temp-f1(:,2)];

end