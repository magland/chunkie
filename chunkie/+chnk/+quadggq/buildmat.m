function [sysmat] = buildmat(chnkr,kern,opdims,type,auxquads,ilist)
%CHNK.QUADGGQ.BUILDMAT build matrix for given kernel and chnkr 
% description of boundary, using special quadrature for self
% and neighbor panels.
%
% Input:
%   chnkr - chunker object describing boundary
%   kern  - kernel function. By default, this should be a function handle
%           accepting input of the form kern(srcinfo,targinfo), where srcinfo
%           and targinfo are in the ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.n - unit normals (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%   opdims - (2) dimension of the kernel, for scalar kernels opdims(1:2) = 1;
%
% Optional input: quantities in brackets indicate default settings
%  type - string ('log'), type of singularity of kernel. Type
%          can take on the following arguments:
%             log => logarithmically singular kernels
%             pv => principal value singular kernels
%             hs => hypersingular kernels
%  auxquads - struct (chnk.quadggq.setuplogquads), structure containing
%             auxilliary quadrature nodes, weights and related
%             interpolation matrices.
%  ilist - cell array of integer arrays ([]), list of panel interactions that 
%          should be ignored when constructing matrix entries or quadrature
%          corrections. 
%
% Ouput:
%   sysmat - the system matrix for discretizing integral operator whose kernel 
%            is defined by kern with a density on the domain defined by chnkr
% 
% NB: if type and auxquads are both present then auxquads will be used 
%       independent of type.

disp('buildmat test 1: entered chnk.quadggq.buildmat');
pause(1);

if (nargin < 3)
    error('not enough arguments in chnk.quadggq.buildmat');
end

if (nargin <6)
    ilist = [];
end

k = chnkr.k;
nch = chnkr.nch;
r = chnkr.r;
adj = chnkr.adj;
d = chnkr.d;
d2 = chnkr.d2;
n = chnkr.n;

fprintf('buildmat test 2: k=%d, nch=%d, opdims=[%d,%d]\n', k, nch, opdims(1), opdims(2));
pause(1);

data = [];
if (chnkr.hasdata)
    data = chnkr.data;
end

disp('buildmat test 3: data resolved');
pause(1);

if nargin < 4 || isempty(type)
    type = 'log';
end
if nargin < 5 || isempty(auxquads)
    disp('buildmat test 4a: auxquads empty, calling chnk.quadggq.setup');
    pause(1);
    auxquads = chnk.quadggq.setup(k,type);
end

fprintf('buildmat test 4: type=%s, auxquads ready\n', type);
pause(1);

temp = eye(opdims(2));

xs1 = auxquads.xs1;
wts1 = auxquads.wts1;
xs0 = auxquads.xs0;
wts0 = auxquads.wts0;

ainterp1 = auxquads.ainterp1;
ainterp1kron = kron(ainterp1,temp);

ainterps0 = auxquads.ainterps0;
ainterps0kron = cell(k,1);
for j = 1:k
    ainterps0kron{j} = kron(ainterps0{j},temp);
end

disp('buildmat test 5: kron interpolation matrices built');
pause(1);

% do smooth weight for all
wts = chnkr.wstor;

disp('buildmat test 6: calling chnk.quadnative.buildmat for smooth fill');
pause(1);
sysmat = chnk.quadnative.buildmat(chnkr,kern,opdims,1:nch,1:nch,wts);
disp('buildmat test 7: smooth fill returned');
pause(1);

% overwrite nbor and self
for j = 1:nch

    fprintf('buildmat test 8: j=%d / %d (overwrite loop)\n', j, nch);
    pause(1);

    jmat = 1 + (j-1)*k*opdims(2);
    jmatend = j*k*opdims(2);

    ibefore = adj(1,j);
    iafter = adj(2,j);

    fprintf('buildmat test 9: j=%d, ibefore=%d, iafter=%d\n', j, ibefore, iafter);
    pause(1);

    % neighbors

    if ibefore > 0
        if ~isempty(ilist) && ismember(ibefore,ilist) && ismember(j,ilist)
        % skip construction if both chunks are in the "bad" chunk list
            fprintf('buildmat test 10: j=%d, skipping ibefore=%d (in ilist)\n', j, ibefore);
            pause(1);
        else
            fprintf('buildmat test 10: j=%d, calling nearbuildmat for ibefore=%d\n', j, ibefore);
            pause(1);
            submat = chnk.quadggq.nearbuildmat(r,d,n,d2,data,ibefore,j, ...
                kern,opdims,xs1,wts1,ainterp1kron,ainterp1);
            fprintf('buildmat test 11: j=%d, ibefore nearbuildmat returned\n', j);
            pause(1);

            imat = 1 + (ibefore-1)*k*opdims(1);
            imatend = ibefore*k*opdims(1);

            sysmat(imat:imatend,jmat:jmatend) = submat;
        end
    end

    if iafter > 0
      if ~isempty(ilist) && ismember(iafter,ilist) && ismember(j,ilist)
        % skip construction if both chunks are in the "bad" chunk list
        fprintf('buildmat test 12: j=%d, skipping iafter=%d (in ilist)\n', j, iafter);
        pause(1);
      else
        fprintf('buildmat test 12: j=%d, calling nearbuildmat for iafter=%d\n', j, iafter);
        pause(1);
        submat = chnk.quadggq.nearbuildmat(r,d,n,d2,data,iafter,j, ...
            kern,opdims,xs1,wts1,ainterp1kron,ainterp1);
        fprintf('buildmat test 13: j=%d, iafter nearbuildmat returned\n', j);
        pause(1);

        imat = 1 + (iafter-1)*k*opdims(1);
        imatend = iafter*k*opdims(1);

        sysmat(imat:imatend,jmat:jmatend) = submat;
      end
    end

    % self
    if ~isempty(ilist) && ismember(j,ilist)
      % skip construction if the chunk is in the "bad" chunk list
      fprintf('buildmat test 14: j=%d, skipping self (in ilist)\n', j);
      pause(1);
    else
      fprintf('buildmat test 14: j=%d, calling diagbuildmat (self)\n', j);
      pause(1);
      submat = chnk.quadggq.diagbuildmat(r,d,n,d2,data,j,kern,opdims,...
          xs0,wts0,ainterps0kron,ainterps0);
      fprintf('buildmat test 15: j=%d, diagbuildmat returned\n', j);
      pause(1);

      imat = 1 + (j-1)*k*opdims(1);
      imatend = j*k*opdims(1);

      sysmat(imat:imatend,jmat:jmatend) = submat;
    end

    fprintf('buildmat test 16: j=%d done\n', j);
    pause(1);

end

disp('buildmat test 17: overwrite loop complete, returning');
pause(1);
	 

end
