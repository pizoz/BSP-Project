function DrawVertMarker(xmk,color,style,marker)
% -------------------------------------------------------------------------------------
%   DrawVertMarker(xmk,color,style,marker)
%
%   Draw a vertical line on the current axis of current figure.
%   xmk    : position on x axis
%   color  : color
%   style  : line style
%   marker : at start and end of the line
%
%   Version 1.00, Date: 01/04/2002
% -------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% -------------------------------------------------------------------------------------

ylimV=get(gca,'Ylim');
hold on
if(~isempty(xmk))
    xmg=[xmk,xmk];
    ymg = ones(size(xmg)); ymg(:,1)=ylimV(1); ymg(:,2)=ylimV(2);
    line(xmg',ymg','Color',color,'LineStyle',style,'Marker',marker);    % draw vertical lines
end
end %== function ================================================================
%

