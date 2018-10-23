function [ax,h]=subtitle(text)
% from https://www.mathworks.com/matlabcentral/answers/...
% 100459-how-can-i-insert-a-title-over-a-group-of-subplots
%Centers a title over a group of subplots.
%Returns a handle to the title and the handle to an axis.
% [ax,h]=subtitle(text)
%           returns handles to both the axis and the title.
% ax=subtitle(text)
%           returns a handle to the axis only.
ax=axes('Units','Normal','Position',[.075 .075 .9 .00],'Visible','off');
set(get(ax,'Title'),'Visible','on', 'FontSize', 14)
title(text);
if (nargout < 2)
    return
end
h=get(ax,'Title');