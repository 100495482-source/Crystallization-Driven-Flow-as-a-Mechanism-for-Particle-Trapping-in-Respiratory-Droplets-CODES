function cm = getCmap(name, trim)
%GETCMAP  A named colormap with its top end trimmed off.
%
%   SYNTAX
%     cm = getCmap(name, trim)
%
%   INPUTS
%     name  Colormap function name, e.g. 'turbo'. Falls back to jet if the
%           name is not available in this MATLAB release.
%     trim  Fraction of the colormap to keep, from the bottom. 0.78 keeps the
%           lower 78 % and discards the rest.
%
%   WHY TRIM
%     The bright end of turbo is nearly white, which disappears against the
%     white figure background used throughout the thesis. Generating a longer
%     colormap and keeping only the lower part preserves 256 usable levels
%     while removing the ones that would be invisible.

try
    n = round(256/max(min(trim,1),0.3));  base = feval(name, n);
catch
    n = round(256/max(min(trim,1),0.3));  base = jet(n);
end
cm = base(1:256,:);
end
