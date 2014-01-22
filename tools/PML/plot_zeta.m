clear ; clc; close all;
%load('zeta.dat')
disp(['Loading zeta.dat ...']);
fid = fopen('zeta.dat','r');
zeta=textscan(fid,'%f %f %f %f %f %f');
zeta=cell2mat(zeta);
disp(['Loading zeta.dat ... done']);
% f=figure
% x=sortrows(zeta,1);
% plot(zeta(:,1),zeta(:,1),'ks');hold on;
% y=sortrows(zeta,2);
% plot(zeta(:,2),zeta(:,2),'rx')
% z=sortrows(zeta,3);
% plot(zeta(:,3),zeta(:,3),'bo')
%%
f=figure; set(f, 'color', 'white'); set(f, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.3, 0.7]); MS=5; LW=0.75; LabSiz=14; FonSiz=13;

xmin=min(abs(zeta(:,1))); xmax=max(abs(zeta(:,1))); 
ymin=min(abs(zeta(:,2))); ymax=max(abs(zeta(:,2))); 
zmin=min(abs(zeta(:,3))); zmax=max(abs(zeta(:,3)));

round(zmax)


subplot(3,1,1); % x-axis nearest
Data1=zeta( abs(zeta(:,2))<=ymin*1.01 & abs(zeta(:,3))<=zmin*1.01 ,:);
zeta_max=max(Data1(:,4))
Data1=sortrows(Data1,1);
set(gca, 'ColorOrder', [0 0 1]);      zeta_x = plot(Data1(:,1)/round(zmax),Data1(:,4)/zeta_max,'.-');hold on;
set(gca, 'ColorOrder', [1 0 0]);      zeta_y = plot(Data1(:,1)/round(zmax),Data1(:,5)/zeta_max,'.-');
set(gca, 'ColorOrder', [0 0.4706 0]); zeta_z = plot(Data1(:,1)/round(zmax),Data1(:,6)/zeta_max,'.-');
legend([zeta_x zeta_y zeta_z], 'zeta_x','zeta_y','zeta_z','Location','EastOutside')
title([' y = ' num2str(ymin) ', z = ' num2str(zmin)],'FontSize',FonSiz*1.6);
ylabel(['\zeta / \zeta_{max} [-]']); xlabel('x / x_{max} [m]');
set(gca,'FontSize',FonSiz);

subplot(3,1,2); % y-axis nearest
Data2=zeta( abs(zeta(:,1))<=xmin*1.01 & abs(zeta(:,3))<=zmin*1.01 ,:);
Data2=sortrows(Data2,2);
set(gca, 'ColorOrder', [0 0 1]);      zeta_x = plot(Data2(:,2)/round(zmax),Data2(:,4)/zeta_max,'.-');hold on;
set(gca, 'ColorOrder', [1 0 0]);      zeta_y = plot(Data2(:,2)/round(zmax),Data2(:,5)/zeta_max,'.-');
set(gca, 'ColorOrder', [0 0.4706 0]); zeta_z = plot(Data2(:,2)/round(zmax),Data2(:,6)/zeta_max,'.-');
legend([zeta_x zeta_y zeta_z], 'zeta_x','zeta_y','zeta_z','Location','EastOutside')
title([' x = ' num2str(xmin) ', z = ' num2str(zmin)],'FontSize',FonSiz*1.6);
ylabel(['\zeta / \zeta_{max} [-]']); xlabel('y / y_{max} [m]');
set(gca,'FontSize',FonSiz);

subplot(3,1,3); % z-axis nearest
Data3=zeta( abs(zeta(:,2))<=ymin*1.01 & abs(zeta(:,1))<=xmin*1.01 ,:);
Data3=sortrows(Data3,3);
set(gca, 'ColorOrder', [0 0 1]);      zeta_x = plot(Data3(:,3)/round(zmax),Data3(:,4)/zeta_max,'.-');hold on;
set(gca, 'ColorOrder', [1 0 0]);      zeta_y = plot(Data3(:,3)/round(zmax),Data3(:,5)/zeta_max,'.-');
set(gca, 'ColorOrder', [0 0.4706 0]); zeta_z = plot(Data3(:,3)/round(zmax),Data3(:,6)/zeta_max,'.-');
legend([zeta_x zeta_y zeta_z], 'zeta_x','zeta_y','zeta_z','Location','EastOutside')
title([' x = ' num2str(xmin) ', y = ' num2str(ymin)],'FontSize',FonSiz*1.6);
ylabel(['\zeta / \zeta_{max} [-]']); xlabel('z / z_{max} [m]');
set(gca,'FontSize',FonSiz);

addpath /home/stephen/Diplomarbeit/PMLalgorithm/export_fig;
export_fig([pwd '/Zeta.pdf']);
close(gcf);

pml = 5
shape = 'sinus'
switch pml
  case 1
    xline = Data1(Data1(:,2)==ymin & Data1(:,3) == zmin,1)
    gauss = xline(89:end)/1.2-5 % intervall ist 1
    xline = xline(89:end)
    xline  = (gauss(end)-gauss(1))*(xline-xline(1))/(xline(end)-xline(1))+gauss(1)  % trafo auf [0.0199, 0.9801]
    
    yline     = Data1(Data1(:,2)==ymin & Data1(:,3) == zmin,4);
    yline     = yline(89:end)/zeta_max;
    
    plot(xline,yline,'.'); set(f, 'color', 'white'); xlabel('\xi / \delta [m]');ylabel(['\zeta / \zeta_{max} [-]']);
    export_fig([pwd '/Zeta_single.pdf']);
    close(gcf);
    
    p7_1PML_sinus = [xline yline]
    p7_1PML_linear = [xline yline]
    p7_1PML_constant = [xline yline]
  case 5
    
    xline     = Data1(Data1(:,2)==ymin & Data1(:,3) == zmin,1)
    gauss = xline(81:88)/1.2*0.2 % intervall ist 0.2
    xline     = xline(121:end)
    n=1
    for I =1:5
      for J=1:8
        index=(I-1)*8+J;
        xline(index) = gauss(J)+(I-1)*0.2
      end
    end
    %xline     = (xline-xline(1))/(xline(end)-xline(1))           % trafo auf [0, 1]
    %xline     = (gauss(end)-gauss(1))*(xline-xline(1))/(xline(end)-xline(1))+gauss(1)  % trafo auf [0.0199, 0.9801]
    yline     = Data1(Data1(:,2)==ymin & Data1(:,3) == zmin,4);
    yline     = yline(121:end)/zeta_max;
    
    plot(xline,yline,'.'); set(f, 'color', 'white'); xlabel('\xi / \delta [m]');ylabel(['\zeta / \zeta_{max} [-]']);
    export_fig([pwd '/Zeta_single.pdf']);
    close(gcf);
    
    p7_5PML_sinus = [xline yline]
    p7_5PML_linear = [xline yline]
    p7_5PML_constant = [xline yline]
end
save([pwd '/p7_' num2str(pml) 'PML_' shape '.mat'],['p7_' num2str(pml) 'PML_' shape]);








return
f=figure;hold on;
subplot(2,2,1);
set(f, 'color', 'white');
h=scatter3(zeta(find(zeta(:,4)>0),1),zeta(find(zeta(:,4)>0),2),zeta(find(zeta(:,4)>0),3),2,'blue','filled');
legend(['Zeta_X = ' num2str(mean(zeta(find(zeta(:,4)>0),4)))],'Location','NorthOutside')
view(40,30)

subplot(2,2,2);
%l=figure;hold on;
%set(l, 'color', 'white')
i=scatter3(zeta(find(zeta(:,5)>0),1),zeta(find(zeta(:,5)>0),2),zeta(find(zeta(:,5)>0),3),2,'red','filled');
legend(['Zeta_Y = ' num2str(mean(zeta(find(zeta(:,5)>0),5)))],'Location','NorthOutside')
view(40,30)

subplot(2,2,3);
%m=figure;hold on;
%set(m, 'color', 'white')
j=scatter3(zeta(find(zeta(:,6)>0),1),zeta(find(zeta(:,6)>0),2),zeta(find(zeta(:,6)>0),3),2,'black','filled');
legend(['Zeta_Z = ' num2str(mean(zeta(find(zeta(:,6)>0),6)))],'Location','NorthOutside')
view(40,30)

subplot(2,2,4);hold on
%m=figure;;
%set(m, 'color', 'white')
h=scatter3(zeta(find(zeta(:,4)>0),1),zeta(find(zeta(:,4)>0),2),zeta(find(zeta(:,4)>0),3),2,'blue','filled');
i=scatter3(zeta(find(zeta(:,5)>0),1),zeta(find(zeta(:,5)>0),2),zeta(find(zeta(:,5)>0),3),2,'red','filled');
j=scatter3(zeta(find(zeta(:,6)>0),1),zeta(find(zeta(:,6)>0),2),zeta(find(zeta(:,6)>0),3),2,'black','filled');
legend(['Zeta = ' num2str(mean(zeta(find(zeta(:,6)>0),6)))],'Location','NorthOutside')
view(40,30)


addpath /home/stephen/Diplomarbeit/Flexi/PMLalgorithm/export_fig;
export_fig('Zeta.pdf');