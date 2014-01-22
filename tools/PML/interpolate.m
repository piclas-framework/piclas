function interpolate_1d=interpolate(E,t)  %E(:,1)=t und E(:,2)=E_abs
index=find(E(:,1)==t);
if ~isempty(index)
  interpolate_1d=E(index,2);
else
  index1=max(find(E(:,1)<=t));
  t1=E(index1,1);
  E1=E(index1,2);
  index2=min(find(E(:,1)>t));
  t2=E(index2,1);
  E2=E(index2,2);
  interpolate_1d=((E2-E1)/(t2-t1))*(t-t1)+E1;
end
end