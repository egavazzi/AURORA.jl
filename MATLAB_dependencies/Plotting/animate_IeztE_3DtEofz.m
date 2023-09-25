function output = animate_IeztE_3DtEofz(t,h_atm,E,Ie_ztE,dE,BeamSA,cx_lims,spp, theta_strs,movieOut)
% animate_IeEztE_3DtzofE - animation of time-varying electron energy flux
% stepping in time. 
% animate_IeztE_3DEzoft animates the Energy and altitude-variation of
% electron energy-fluxes time-step by time-step from the start to the
% finish. The electron spectra are plotted in log-scale
% with pcolor in a number of subplots.
% 
% Calling:
%  clims = animate_IeztE_3DtEofz(t,h_atm,E,Ie_ztE,dE,BeamSA,cx_lims,spp, theta_strs)
% Input:
%  t,         - time-scale (s), double array [1 x n_t]
%  h_atm      - altitudes (km), double array [n_z x 1]
%  E          - energies (eV), double array [1 x n_E]
%  Ie_ztE     - electron number-flux, double array [n_z,n_t,n_E]
%  dE         - energy bin sizes (eV), double array [1 x n_E]
%  BeamSA     - solid angle sizes of the beams, [n_beams x 1]
%  cx_lims    - limit for colour-scale double array [1 x 2],
%               optional argument, if left empty the fluxes of each
%               beam will be scaled thusly: caxis([-4 0]+max(caxis))
%  spp        - sub-plot-position, integer array [n_beams x 3] with
%               sub-plot layout and subplot position row-by-row,
%               the function is designed for two rows of sub-plots
%               with the idea that downward fluxes should be
%               plotted in the top row and the upward in the bottom
%               row.
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

%   Copyright © 2018 Björn Gustavsson, <bjorn.gustavsson@uit.no>
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

z_prev = 1e23;
for i_z = numel(h_atm):-1:1,
  clf
  if h_atm(i_z)/1e3 - z_prev < -1
    z_prev = h_atm(i_z)/1e3;
    if ~isempty(cx_lims)
      subplot(1,1,1)
      caxis(cx_lims)
      cblh = colorbar_labeled('log10(#e/m2/s/eV/ster)');
    end
    for i1 = 1:size(spp,1),     
      % disp([i1, spp(i1,:)])
      % log10(max(0,real(squeeze(Ie_ztE(i_z+(i1-1)*numel(h_atm),:,:))./(ones(size(h_atm))*dE)/BeamSA(i1))))),
      subplot(spp(i1,1),spp(i1,2),spp(i1,3))
      %dx = -0.02*mod(spp(i1,3)-1,spp(i1,2));
      %set(gca,'position',get(gca,'position')+[dx,0,0,0])
      if i1*numel(h_atm) <= size(Ie_ztE,1)
        pcolor(t,...
               E,...
               log10(max(0,real(squeeze(Ie_ztE(i_z+(i1-1)*numel(h_atm),:,:))'./repmat(dE(:),size(t))/BeamSA(i1))))),
        shading flat
        if isempty(cx_lims)
          clims(i_z,i1,:) = caxis;
          caxis([-4 0]+max(caxis))
        else
          caxis(cx_lims)
        end
        set(gca,'tickdir','out','yscale','log','ytick',[1,10,1000,10000,1e5])
        if mod(i1,spp(1,2))==0 || i1 == size(spp,1)
          % colorbar_labeled('log10(#e/m2/s/eV/ster)')
        elseif isempty(cx_lims)
          colorbar_labeled('')      
        end
        if mod(spp(i1,3),spp(1,2))==1 %i1 == 1 || i1 == 4
          ylabel('energy (eV)')
        else
          set(gca,'yticklabel','')        
        end
        if spp(i1,3) > (spp(1,2)*(spp(1,1)-1)) % spp(1,2)
          xlabel('time (s)')
        else
          set(gca,'xticklabel','')
        end
        if spp(i1,3) == ceil(spp(1,2)/2) % was i1 == 2
                                  % if i1 == 2
          sgtitle(sprintf('%4.3f (km)',h_atm(i_z)/1e3))
%         else
        end
          % title(sprintf('pitch-angle ~ %s',theta_strs{i1}))
        tstr = sprintf('\\theta_B ~ %s',theta_strs{i1});
        title(tstr)
%           title(['\theta-B ~ ',theta_strs{i1}])
%         end
      else
        
      end
    end
%     [~,iS] = sort(spp(:,3)); 
%     for i1 = 1:size(spp,1),                             
%       subplot(spp(iS(i1),1),spp(iS(i1),2),spp(iS(i1),3))
%       dx = -0.02*mod(spp(iS(i1),3)-1,spp(iS(i1),2));
%       set(gca,'position',get(gca,'position')+[dx,0,0,0])
%     end
%     set(cblh,'position',get(cblh,'position')+[dx,0,0,0])
    drawnow
    if doMovie == 1 || doMovie == 2
      if doMovie == 1
        writeVideo(vidObj,getframe(gcf));
      elseif doMovie == 2
        M(j) = getframe(gcf);
      end
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
