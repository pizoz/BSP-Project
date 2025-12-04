function figResize(px, py, pwidth, pheigth)
% ---------------------------------------------------------------------------------------------
% Figure resizing and positioning
% px      = orizzontal position in fraction of screen width
% py      = vertical position in fraction of screen heigth
% pwidth  = width in fraction of screen width
% pheigth = heigth in fraction of screen heigth
% example: figResize(0, 1, 1, .3);
%
% ---------------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% ---------------------------------------------------------------------------------------------


oldpos=get(gcf,'Position');
fgxpos=oldpos(1); fgypos=oldpos(2); fgwidth=oldpos(3); fgheigth=oldpos(4);

set(gcf,'Color','white');
hTitle=14+ 10;
hXlabel=12+ 4;
hYlabel=10;
hTlabel=hTitle+hXlabel;
bdwidth=5;
topbdwidth=70;
set(0,'Units','pixels')
scnsize = get(0,'ScreenSize');

if(nargin>3 && ~isempty(pheigth)), fgheigth = scnsize(4)*pheigth; end
if(nargin>2 && ~isempty(pwidth)), fgwidth = scnsize(3)*pwidth - 2*bdwidth; end
if(nargin>1 && ~isempty(py)), fgypos = scnsize(4)*py - fgheigth - (topbdwidth+bdwidth); end
if(nargin>0 && ~isempty(px)), fgxpos = scnsize(3)*px + bdwidth; end
pos = [fgxpos, fgypos, fgwidth, fgheigth ];
% change position and size of the figure.
pause(0.1);   % it need on my PC
set(gcf,'Position', pos);

% % get and change the position and length of the main ordinate axis.
% pos = get(gca,'Position');
% hrTitle= hTitle/fgheigth;
% hrXlabel= hXlabel/fgheigth;  hrXlabel=0;
% hrYlabel= hYlabel/fgwidth;   hrYlabel=0;
% % change the position and length of the main ordinate axis.
% set(gca,'Position',[pos(1)+hrYlabel, pos(2)+hrXlabel pos(3), pos(4)-hrTitle]);

end %== function ================================================================
%

