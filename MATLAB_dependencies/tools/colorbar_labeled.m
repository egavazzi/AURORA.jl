function handle=colorbar_labeled(label,lg,varargin)
% COLORBAR_LABELED - colorbar with label for range linear or log
% scale. 
% 
% Calling:
%   handle=colorbar_labeled(label,lg)
% Input:
%   label - string with text
%   lg    - string [ {'linear'} | 'log' ] for scale type

if nargin < 2 || isempty(lg),
  lg = 'linear';
end

% Determine color limits by context.  If any axes child is an image
% use scale based on size of colormap, otherwise use current CAXIS.

t = caxis;

h = gca;

% Search for existing colorbar
ch = get(gcf,'children');
ax = [];
try
  for i = 1:length(ch),
    d = get(ch(i),'userdata');
    if numel(d) == 1 && d == h
      ax = ch(i);
      break
    end
  end
catch
  % Dont know what to do around here...
end

if strcmp(get(gcf,'NextPlot'),'replace'),
  set(gcf,'NextPlot','add')
  wasNextplotReplace = 1;
else
  wasNextplotReplace = 0;
end


if isempty(ax)
  pos = get(h,'Position');
% stripe=0.03; edge=0.08; 
% [az,el]=view;
% if all([az,el]==[0 90]), space=0.02; else, space=.1; end
% set(h,'Position',[pos(1) pos(2) pos(3)*(1-stripe-edge-space) pos(4)])
% rect=[pos(1)+(1-stripe-edge)*pos(3) pos(2) stripe*pos(3) pos(4)];
  rect = [pos(1)+1.02*pos(3) pos(2) .035*pos(3) pos(4)];

  % Create axes for stripe
  ax = axes('Position', rect);
else
  axes(ax)
end

if ~isempty(varargin)
  fszi = find_in_cellstr('fontsize',varargin);
  if ~isempty(fszi)
    fsz = varargin{fszi+1}-2;
  else
    fsz = get(gca,'fontsize');
  end
else
  fsz = get(gca,'fontsize');
end
% Create color stripe
cmap = colormap;
n = size(cmap,1);
if strcmp(lg,'log')
  tt=((0:n)-.5)'*diff(t)/(n-1)+t(1);
  surface([0 1],10.^tt,[tt tt]),shading flat
  set(ax,'CLim',t,'ylim',10.^t,'yscale','log','layer','top','fontsize',fsz)
  %surface([0 1],tt,[tt tt]),shading flat
  %set(ax,'CLim',t,'ylim',t,'yscale','log','layer','top','fontsize',fsz)
else
  %image([0 1],t,[1:n]'),
  tt=((0:n)-.5)'*diff(t)/(n-1)+t(1);
  surface([0 1],tt,[tt tt]),shading flat
  set(ax,'CLim',t,'ylim',t,'yscale','linear','layer','top','fontsize',fsz)
% set(ax,'TickDir','in')
  set(ax,'ylim',t,'Ydir','normal')
% Justify color axis
% set(ax,'yticklabel',strjust(get(ax,'yticklabel')))
end

set(ax,'userdata',h,'YaxisLoc','right','xtick',[])
  ylabel(label,'Rotation',-90,'VerticalAlignment','baseline',varargin{:})
set(gcf,'CurrentAxes',h)

if wasNextplotReplace == 1
  set(gcf,'Nextplot','Replace')
end

if nargout > 0, 
  handle = ax;
end
