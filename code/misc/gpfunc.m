function L = gpfunc(varargin)

HYP = varargin{1};
im  = varargin{2};   
%keyboard
par = varargin(3:end-1);

hypstruc = varargin{end};
%T = varargin{end};
%try
HYP = rewrap(hypstruc,HYP);
%catch
%    keyboard
%end
%keyboard
[L dnlZ] = gp(HYP, im, par{:});
L = L;
