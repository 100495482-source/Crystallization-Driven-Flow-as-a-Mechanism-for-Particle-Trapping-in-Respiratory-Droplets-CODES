function drawTrackSet(ax, xx, yy, tidx, vv, tlist, climV, col, lw, alph, cmap)
%DRAWTRACKSET  Draw a set of trajectories as a single patch object.
%
%   SYNTAX
%     drawTrackSet(ax, xx, yy, tidx, vv, tlist, climV, col, lw, alph, cmap)
%
%   INPUTS
%     ax     Target axes
%     xx,yy  Spot coordinates, sorted by track and then by frame
%     tidx   Track index of each spot
%     vv     Speed at each spot (um/s); may contain NaN
%     tlist  Track indices to draw
%     climV  Colour limits, when colouring by speed
%     col    Flat colour, when not colouring by speed
%     lw     Line width
%     alph   Edge alpha, for the flat-colour case
%     cmap   Colormap for speed colouring, or [] for a flat colour
%
%   WHY A SINGLE PATCH
%     A droplet has thousands of tracks. Calling plot once per track produces
%     thousands of graphics objects and makes the figure unusably slow to
%     render and export. Building one patch whose faces are the individual
%     segments draws the whole set in a single object.
%
%     Segments are only created between consecutive spots BELONGING TO THE
%     SAME TRACK, which is what the tidx comparison below enforces; otherwise
%     the patch would connect the end of one trajectory to the start of the
%     next.

if isempty(tlist), return; end
m = ismember(tidx, tlist);
if ~any(m), return; end

X = xx(m); Y = yy(m); Tq = tidx(m); V = vv(m);
n = numel(X);  i1 = (1:n-1)';
ok = Tq(i1) == Tq(i1+1);   % no segment across a track boundary
i1 = i1(ok);
if isempty(i1), return; end

if isempty(cmap)
    patch(ax, 'Faces',[i1 i1+1], 'Vertices',[X Y], ...
          'EdgeColor',col, 'EdgeAlpha',alph, 'FaceColor','none', 'LineWidth',lw);
else
    cd = V;  cd(~isfinite(cd)) = climV(1);   % NaN speeds go to the low end
    patch(ax, 'Faces',[i1 i1+1], 'Vertices',[X Y], ...
          'FaceVertexCData',cd, 'EdgeColor','interp', 'FaceColor','none', ...
          'LineWidth',lw);
    colormap(ax, cmap); clim(ax, climV);
end
end
