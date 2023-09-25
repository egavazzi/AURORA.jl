function output = animate_IeztE_pitchangleflux(t,h_atm,E,dE,IeZTE,BeamW,mu_lims,iZ,movieOut,cxin)
% animate_IeztE_pitchangleflux - animation of time-varying electron energy flux
%  animate_IeztE_pitchangleflux animates the time-varying electron
%  fluxes varying with energy - pitch-angle at four altitudes
%  stepping in time. The electron-fluxes are plotted in log-scale
%  with polarPcolor in four subplots.
%
% Calling:
%  output = animate_IeztE_pitchangleflux(t,h_atm,E,dE,IeZTE,BeamW,mu_lims,iZ,movieOut)
% Input:
%  t          - time-scale (s), double array [1 x n_t]
%  h_atm      - altitudes (km), double array [n_z x 1]
%  E          - energies (eV), double array [1 x n_E]
%  Ie_ztE     - electron number-flux, double array [n_z*n_beams,n_t,n_E]
%  dE         - energy bin sizes (eV), double array [1 x n_E]
%  BeamSA     - solid angle sizes of the beams, [n_beams x 1]
%  mu_lims    - beam-boundaries, (cos(theta_b)), [n_beams x 1] or
%               [1 x n_beams]
%  movieOut   - Optional argument, if MovieName in matlab-supported
%               video-format then animate_IeEztE_3DtEofz will attempt
%               to write the displayed sequence to a video-file
%               with that filename using the VideoWriter
%               functionality. See VideoWriter for details. If
%               numerical value that is converted to TRUE, a matlab
%               movie will be returned.
% Output:
%  M     - Matlab-movie with the displayed sequence.


%   Copyright � 20190506 Bj�rn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later

doMovie = 0;

if nargin > 7 && ~isempty(movieOut)
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

if nargin < 10 || isempty(cxin)
  cx2use = ([-5 0]+12.25);
else
  cx2use = cxin;
end
nZ = numel(h_atm);
%% Extracting electron flux at 4 altitudes

for iz = numel(iZ):-1:1,
  Ie_alts{iz} = IeZTE(iZ(iz):nZ:end,:,:);
end

%% Graphics settings:
ftsz = 14;
colormap(parula)
xtxtpos = [0,1.3]; % Altitude-label above change 1.3 to -1.3 for below
suptpos = [1.1,1.4];
cblpos = [-0.7,-0.3];

%% Animationing
for it = 1:numel(t),
  caxis(cx2use)
  clf
  for ip = 1:4
    IePAD{ip}(:,:,it) = squeeze(Ie_alts{5-ip}([1:end,end],it,:))./...
                          (dE(1:size(IeZTE,3))'*BeamW([1:end,end]))';
    subplot(2,2,ip)
    polarPcolor(log10(E(1:size(IeZTE,3))),...
                acos(mu_lims)*180/pi,...
                log10(max(cx2use(1),squeeze(Ie_alts{5-ip}([1:end,end],it,:))./...
                          (dE(1:size(IeZTE,3))'*BeamW([1:end,end]))')),...
                'Ncircles',3,...
                'Nspokes',7,...
                'Rvals',[1 2 3],...
                'Rticklabel',{'10^1','10^2','10^3'});
    if ip == 1
      text(suptpos(1),suptpos(2),sprintf('time: %4.3f (s)',t(it)),'fontsize',ftsz)
    end
    alt_str = sprintf('z: %2.1f (km)',h_atm(iZ(5-ip)));
    % text(xtxtpos(1),xtxtpos(2),'height: 610 km','fontsize',ftsz)
    text(xtxtpos(1),xtxtpos(2),alt_str,'fontsize',ftsz)
    text(cblpos(1),cblpos(2),'(e^{-1}/m^2/s/ster/eV)','rotation',90)
    caxis(cx2use)
  end
  drawnow
  if doMovie == 1 || doMovie == 2
    if doMovie == 1
      writeVideo(vidObj,getframe(gcf));
    elseif doMovie == 2
      M(it) = getframe(gcf);
    end
  end
end

if doMovie == 2
  output = M;
elseif doMovie == 1
  output = [];
  close(vidObj);
else
  output = IePAD;
end
