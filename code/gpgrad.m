function dnlZ = gpgrad(varargin)

HYP = varargin{1}';

im = varargin{2};   
par = varargin(3:end-1);

hypstruc = varargin{end};
%T = varargin{end};
%par        = {'meanConst','covBranchingProcess_2A','likGauss',X2,Y1,hyp};

HYP = rewrap(hypstruc,HYP);
[L dnlZ] = gp(HYP, im, par{:});
dnlZ = unwrap(dnlZ)';