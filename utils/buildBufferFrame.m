function O = buildBufferFrame(cen, Rbuf, k, Ny, Nx)
%BUILDBUFFERFRAME  Region of influence of each crystal in one frame.
%
%   Returns a uint8 label image: 0 outside every buffer, c inside the buffer
%   belonging to crystal c. No pixel is ever shared and none is left orphaned.
%
%   SYNTAX
%     O = buildBufferFrame(cen, Rbuf, k, Ny, Nx)
%
%   INPUTS
%     cen   Centroids, [frames x NCrys x 2], in pixels
%     Rbuf  Fixed buffer radius per crystal, [1 x NCrys], in pixels
%     k     Frame index (1-based)
%     Ny,Nx Image size
%
%   HOW OVERLAPS ARE RESOLVED
%     Where two or more discs overlap, each pixel is assigned to the circle
%     that minimises d^2 - r^2, its power with respect to that circle. For a
%     pair of circles that criterion traces exactly the line through their
%     two intersection points, that is, their common chord. With equal radii
%     it reduces to the perpendicular bisector; with unequal radii the chord
%     shifts towards the smaller circle, which is the sensible behaviour: a
%     small crystal should not claim territory from a large one.
%
%     This is the power diagram of the discs, and it partitions the contested
%     region cleanly, which matters because the accumulation statistics count
%     each particle once and cannot tolerate double counting.
%
%   WHY THE BUFFER IS NOT STORED PER FRAME
%     Saving the mask for every frame of every droplet would run to gigabytes.
%     Only the parameters are stored (centroid per frame, radius per crystal)
%     and this function rebuilds any frame on demand.
%
%   See also ACCUM3_BUILDBUFFERZONES.

O = zeros(Ny, Nx, 'uint8');
cxk = reshape(cen(k,:,1), 1, []);
cyk = reshape(cen(k,:,2), 1, []);
pres = find(isfinite(Rbuf(:)') & isfinite(cxk) & isfinite(cyk));
if isempty(pres), return; end

% Bounding box of all the discs, so the whole image is not scanned.
xmin = inf; xmax = -inf; ymin = inf; ymax = -inf;
for c = pres
    xmin = min(xmin, cxk(c)-Rbuf(c));  xmax = max(xmax, cxk(c)+Rbuf(c));
    ymin = min(ymin, cyk(c)-Rbuf(c));  ymax = max(ymax, cyk(c)+Rbuf(c));
end
c0 = max(1, floor(xmin));  c1 = min(Nx, ceil(xmax));
r0 = max(1, floor(ymin));  r1 = min(Ny, ceil(ymax));
if c1 < c0 || r1 < r0, return; end

[XX, YY] = meshgrid(c0:c1, r0:r1);
best = inf(size(XX));  own = zeros(size(XX), 'uint8');
for c = pres
    d2 = (XX - cxk(c)).^2 + (YY - cyk(c)).^2;
    pw = d2 - Rbuf(c)^2;                 % power with respect to this circle
    m  = d2 <= Rbuf(c)^2 & pw < best;    % inside the disc and closer
    best(m) = pw(m);  own(m) = c;
end
O(r0:r1, c0:c1) = own;
end
