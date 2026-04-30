function submat = buildmat(chnkr,kern,opdims,i,j,wts)
%CHNK.QUADSMOOTH.BUILDMAT build matrix for far interactions with this kernel
% assuming that the smooth rule is sufficient
% 
%


disp('quadnative buildmat test 1: entered chnk.quadnative.buildmat');
pause(1);

if nargin < 4
    i = 1:chnkr.nch;
end
if nargin < 5
    j = 1:chnkr.nch;
end

if nargin < 6
    [~,wts] = lege.exps(chnkr.k);
end

fprintf('quadnative buildmat test 2: nargs resolved, length(i)=%d, length(j)=%d\n', length(i), length(j));
pause(1);

% grab specific boundary data

r = chnkr.rstor;
d = chnkr.dstor;
d2 = chnkr.d2stor;
n = chnkr.nstor;
data = chnkr.datastor;

disp('quadnative buildmat test 3: stor arrays grabbed');
pause(1);

[dim,k,~] = size(r);
fprintf('quadnative buildmat test 4: dim=%d, k=%d\n', dim, k);
pause(1);

rs = r(:,:,j); rt = r(:,:,i); ds = d(:,:,j); dt = d(:,:,i); nt = n(:,:,i);
ns = n(:,:,j);
d2s = d2(:,:,j); d2t = d2(:,:,i);
rs = reshape(rs,dim,k*length(j)); rt = reshape(rt,dim,k*length(i));
ds = reshape(ds,dim,k*length(j)); dt = reshape(dt,dim,k*length(i));
d2s = reshape(d2s,dim,k*length(j)); d2t = reshape(d2t,dim,k*length(i));

disp('quadnative buildmat test 5: src/targ arrays sliced and reshaped');
pause(1);

srcinfo = []; srcinfo.r = rs; srcinfo.d = ds; srcinfo.d2 = d2s; srcinfo.n = ns;
targinfo = []; targinfo.r = rt; targinfo.d = dt; targinfo.d2 = d2t; targinfo.n = nt;

if (not (isempty(data)))
    srcinfo.data = data(:,:,j);
    targinfo.data = data(:,:,i);
end

disp('quadnative buildmat test 6: srcinfo/targinfo built');
pause(1);

dsnrms = sqrt(sum(ds.^2,1));
ws = repmat(wts(:),length(j),1);

dsdt = dsnrms(:).*ws;

dsdtndim2 = repmat(dsdt(:).',opdims(2),1);
dsdtndim2 = dsdtndim2(:);

disp('quadnative buildmat test 7: weight vector built, calling kern');
pause(1);

submat = bsxfun(@times,kern(srcinfo,targinfo),(dsdtndim2(:)).');

disp('quadnative buildmat test 8: kern returned, submat built, returning');
pause(1);

end
