clear
close
clc


data=h5read('twt_RP_000.0000000002000000_reference.h5','/RP_Data');

dimensions=size(data);

Nvar=dimensions(1);
N_RP=dimensions(2);
N=dimensions(3);

time=zeros(N,1);
RP_data=zeros(N,N_RP,Nvar-1);

for I=1:N
    time(I)=data(1,1,I);
    for J=1:N_RP
        for K=1:Nvar-1
            RP_data(I,J,K)=data(K+1,J,I);
        end
    end
end

figure; hold on
style='rgbcmrgbcmrgbcm';
for J=1:N_RP
    plot(time,RP_data(:,J,1));
end

