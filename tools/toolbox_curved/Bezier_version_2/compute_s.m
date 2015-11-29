function compute_s()
%sarray=zeros(n,2);
% sarray(1,1)=-0.75;
% sarray(1,2)=0.25;
% sarray(2,1)=-0.5;
% sarray(2,2)=0.5;
% sarray(3,1)=-0.5;
% sarray(3,2)=0.0;
% sarray(4,1)=-1.0;
% sarray(4,2)=1.0;
sarray(1,1)=-0.6667;
sarray(2,1)=-1.0000;
sarray(3,1)=-0.4202;
sarray(4,1)=-0.0283;
sarray(5,1)=0.1790;
sarray(6,1)=-0.9292;
sarray(7,1)=-0.2904;

sarray(1,2)=1.0000;
sarray(2,2)=-0.1218;
sarray(3,2)=-0.0680;
sarray(4,2)=1.0000;
sarray(5,2)=0.3801;
sarray(6,2)=-0.9277;
sarray(7,2)=-0.2902;


dimension=size(sarray,1);

smin=sarray(dimension,1);
smax=sarray(dimension,2);
for i=dimension-1:-1:1;
  smin=LinIntPol1D(sarray(i,1),sarray(i,2),smin);
  smax=LinIntPol1D(sarray(i,1),sarray(i,2),smax);
end
disp(['smin: ',num2str(smin),' smax: ',num2str(smax),' mean: ',num2str(0.5*(smin+smax))])

end

function result=LinIntPol1D(y1,y2,x)
m=(y2-y1)/2;
result=m*(x+1)+y1;
end 