function getPlotFramework(x1,x2,x3,x4,x5,x6,x7,x8)
global ha
global input;global x_p1;global y_p1;global z_p1;global x_p2;global y_p2;global z_p2
%ha=subplot(10,1,1:6);  hold on;set(gcf, 'color', 'white');view(79.5,22);grid on;
axis tight
if true(input.removeAx), axis off; set(gca,'XTick',[]); set(gca,'YTick',[]); set(gca,'ZTick',[]); end

%xlim([-3.5 4.5]);ylim([-3.5 4.5]);zlim([-4. 5]);
text(x1(1)-0.5,x1(2),x1(3)-0.1,'x_1');
text(x2(1),x2(2)-0.1,x2(3)-0.1,'x_2');
text(x3(1)+0.1,x3(2),x3(3)-0.1,'x_3');
text(x4(1),x4(2)+0.05,x4(3)-0.1,'x_4');
text(x5(1),x5(2),x5(3)+0.1,'x_5');
text(x6(1)+0.3,x6(2),x6(3)+0.1,'x_6');
text(x7(1),x7(2),x7(3)+0.1,'x_7');
text(x8(1)+0.3,x8(2),x8(3)+0.1,'x_8');
xlabel('x');ylabel('y');zlabel('z');
quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[0.1;0;0],[0;0.05;0],[0;0;0.1],0,'LineWidth',1)
plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'ko-')
plot3([x3(1) x2(1)],[x3(2) x2(2)],[x3(3) x2(3)],'ko-')
plot3([x3(1) x4(1)],[x3(2) x4(2)],[x3(3) x4(3)],'ko-')
plot3([x1(1) x4(1)],[x1(2) x4(2)],[x1(3) x4(3)],'ko-')

plot3([x1(1) x5(1)],[x1(2) x5(2)],[x1(3) x5(3)],'ko-')
plot3([x6(1) x2(1)],[x6(2) x2(2)],[x6(3) x2(3)],'ko-')
plot3([x3(1) x7(1)],[x3(2) x7(2)],[x3(3) x7(3)],'ko-')
plot3([x8(1) x4(1)],[x8(2) x4(2)],[x8(3) x4(3)],'ko-')

plot3([x6(1) x5(1)],[x6(2) x5(2)],[x6(3) x5(3)],'ko-')
plot3([x6(1) x7(1)],[x6(2) x7(2)],[x6(3) x7(3)],'ko-')
plot3([x8(1) x7(1)],[x8(2) x7(2)],[x8(3) x7(3)],'ko-')
plot3([x8(1) x5(1)],[x8(2) x5(2)],[x8(3) x5(3)],'ko-')
end
