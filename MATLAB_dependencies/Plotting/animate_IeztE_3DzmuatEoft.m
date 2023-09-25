function output = animate_IeztE_3DzmuatEoft(t,h_atm,E,Ie_ztE,dE,BeamSA,iE,cx_lims,theta_strs,movieOut)
% animate_IeztE_3DzmuatEoft - animation of time-varying electron flux
% stepping in time. 
% animate_IeztE_3DzmuatEoft animates the altitude-pitch-angle-variation of
% electron number-fluxes at a fixed time-step by time-step from the
% start to the finish. The electron spectra are plotted in log-scale
% with pcolor.
% 
% Calling:
%  clims = animate_IeztE_3DzmuatEoft(t,h_atm,E,Ie_ztE,dE,BeamSA,iE,cx_lims, theta_strs[,movieOut])
% Input:
%  t          - time-scale (s), double array [1 x n_t]
%  h_atm      - altitudes (km), double array [n_z x 1]
%  E          - energies (eV), double array [1 x n_E]
%  Ie_ztE     - electron number-flux, double array [n_z,n_t,n_E]
%  dE         - energy bin sizes (eV), double array [1 x n_E]
%  BeamSA     - solid angle sizes of the beams, [n_beams x 1]
%  iE         - index for the energy to plot, scalar integer in [1 n_E]
%  cx_lims    - limit for colour-scale double array [1 x 2],
%               optional argument, if left empty the fluxes of each
%               beam will be scaled thusly: caxis([-4 0]+max(caxis))
%  theta_strs - titles for the subplots, cell-array with
%               pitch-angle-boundaries, has to be at least as long
%               as the number of sub-plots plotted.
%  movieOut   - Optional argument, if MovieName in matlab-supported
%               video-format then animate_IeEztE_3DtEofz will attempt
%               to write the displayed sequence to a video-file
%               with that filename using the VideoWriter
%               functionality. See VideoWriter for details. If
%               numerical value that is converted to TRUE, a matlab
%               movie will be returned.
% Output:
%  clims - 3-D array with color-scales for all sub-plots plotted.
%  M     - Matlab-movie with the displayed sequence.
% 
% Example: See Etrp_example_7stream.m

%   Copyright © 2020 Björn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later

nZ = numel(h_atm);
clims = [];
doMovie = 0;

if nargin > 9 && ~isempty(movieOut)
  if isstr(movieOut)
    try
      vidObj = VideoWriter(movieOut);
      open(vidObj);
      doMovie = 1;
    catch
      doMovie = 0;
      disp(['Failed to initialize VideoWriter Object for movie: ',movieOut])
    end
  elseif movieOut
    doMovie = 2;
  end
end

subplot(1,1,1)
cblh = colorbar_labeled('log10(#/m2/s/eV/ster)');
cbhp = get(cblh,'position');
clf
for i_t = 1:numel(t),
  Ie_curr = squeeze(Ie_ztE(:,i_t,iE));
  Ie_curr = reshape(Ie_curr,numel(h_atm),[]);
  pcolor(0.5:1:(0.5+numel(BeamSA)),...
         h_atm/1e3,...
         asinh((Ie_curr(:,[1:end,end])./repmat(BeamSA([1:end,end]),size(h_atm))))),
  shading flat
  if isempty(cx_lims)
    clims(i_t,:) = caxis;
    %caxis([-4 0]+max(caxis))
  else
    caxis(cx_lims)
  end
  set(gca,...
      'tickdir','out',...
      'xscale','linear',...
      'xtick',[0.5:1:(0.5+numel(BeamSA))],...
      'xticklabel',theta_strs,...
      'box','off')
  cblh = colorbar_labeled('log10(#/m2/s/eV/ster)');
  set(cblh,'position',cbhp+[-0.02,0,0,0])
  ylabel('height (km)')
  xlabel('pitch-angles')
  title(sprintf('%4.3f (s)',t(i_t)))
  drawnow
  if doMovie == 1 || doMovie == 2
    if doMovie == 1
      writeVideo(vidObj,getframe(gcf));
    elseif doMovie == 2
      M(j) = getframe(gcf);
    end
  end
end

if doMovie == 2
  output = M;
elseif doMovie == 1
  close(vidObj);
  output = [];
else
  output = clims;
end
