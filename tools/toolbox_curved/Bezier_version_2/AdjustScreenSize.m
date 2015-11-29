function AdjustScreenSize()

monitor='small';
set(0,'Units','pixels')
scnsize = get(0,'ScreenSize');
position = get(gcf,'Position');
outerpos = get(gcf,'OuterPosition');
borders = outerpos - position;
edge = -borders(1)/2;
if strcmp(monitor,'large'),pos1=[outerpos(1)*5.5,scnsize(4),scnsize(3)/3-edge,scnsize(4)/1.25];
  pos2=[outerpos(1)*6.7,scnsize(4),scnsize(3)/3-edge,scnsize(4)/1.25];
  pos3=[outerpos(1)*5.5,-outerpos(1)/1.5,scnsize(3)/3-edge,scnsize(4)/2];
  pos4=[outerpos(1)*6.7,-outerpos(1)/1.5,scnsize(3)/3-edge,scnsize(4)/2];
end;

if strcmp(monitor,'small'),pos1=[1400,550,scnsize(3)/6,scnsize(4)/2.2];
  pos2=[1950,550,scnsize(3)/6,scnsize(4)/2.2];
  pos3=[1400,50,scnsize(3)/7,scnsize(4)/2.2];
  pos4=[1800,50,scnsize(3)/7,scnsize(4)/2.2];
end;
set(gcf,'OuterPosition',pos1)
end