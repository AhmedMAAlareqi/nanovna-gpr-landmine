function bscan_analyze(label, bg)
%BSCAN_ANALYZE  Radargrams before/after processing for one sandbox B-scan
%   bscan_analyze('water3')   (reads bscan_water3.mat and bscan_kernel.mat)
%
%   Panels:  (1) raw radargram
%            (2) B1  mean-trace subtraction
%            (3) B2  SVD, first component removed
%            (4) SR-SPD residual: median trace removed + surface-change subspace
%                [K, uK, u^2K] projected out per trace  (u = f/fc - 1)
%            (5) SR-SPD sparse reflectivity map (red = like metal/water,
%                blue = lower permittivity than the sand: plastic-mine case)
%            (6) detection statistic vs position for B1, B2 and SR-SPD
%   Appends one line per scan to bscan_results.csv.
%
%   Optional background scans (fixed container/geometry pattern removed first):
%     bscan_analyze('water3', {'empty_a','empty_b'})   % target vs mean of empties
%     bscan_analyze('empty_a', {'empty_b'})            % fair false-alarm check (leave one out)

if nargin < 1, label = input('Scan label (e.g. water3): ', 's'); end
label = regexprep(strtrim(label), '^(bscan_)+', '');
P.nd = 2;              % surface-change terms (1 = small height changes, 2 = larger)
P.win = [0.15 2.0];    % target window below the surface, ns (bottom echo comes later)
P.M = 2;               % reflectors per trace in the sparse fit
P.show = [-2 5];       % time range shown, ns relative to the surface

B  = load(['bscan_' label '.mat']);  f = B.f(:).';  X = B.S;  x = B.x_cm(:).';  [T, n] = size(X);
tag = '';  rawTitle = 'Raw radargram';
if nargin >= 2 && ~isempty(bg)
    if ischar(bg), bg = {bg}; end
    Xb = zeros(size(X));
    for q = 1:numel(bg)
        Bq = load(['bscan_' regexprep(bg{q}, '^(bscan_)+', '') '.mat']);
        if ~isequal(size(Bq.S), size(X)), error('Background %s has a different size.', bg{q}); end
        Xb = Xb + Bq.S / numel(bg);
    end
    X = X - Xb;  tag = [' [bg: ' strjoin(bg, '+') ']'];  rawTitle = 'Input after background subtraction';
end
Kf = load('bscan_kernel.mat');  K = Kf.K(:).';
N  = 4096;  w = hann1(n);  t = (0:N-1)/(N*(f(2)-f(1)));
td = @(Y) ifft(Y .* w, N, 2);
[~, ik] = max(abs(td(K)));  t0 = t(ik);                 % surface time (from the kernel)
tr = (t - t0)*1e9;  ti = tr >= P.show(1) & tr <= P.show(2);
inw = tr >= P.win(1) & tr <= P.win(2);

% ---- baselines ----
R1 = X - mean(X, 1);                                   % B1
[U, Sg, V] = svd(X, 'econ');  R2 = X - U(:,1)*Sg(1,1)*V(:,1)';   % B2

% ---- SR-SPD ----
Mref = median(real(X), 1) + 1j*median(imag(X), 1);     % robust common part (coupling + mean ground)
u  = f/mean(f) - 1;
A  = K.';  for m = 1:P.nd, A = [A, (u.^m .* K).']; end %#ok<AGROW>
Pr = eye(n) - A*pinv(A);                               % removes surface moving / changing strength
R3 = (Pr * (X - Mref).').';
tau = (P.win(1):0.02:P.win(2))*1e-9;
At  = Pr * (K.' .* exp(-1j*2*pi*f.'*tau));             % shifted kernels, same projection
An  = At ./ vecnorm(At);
Cmap = zeros(T, numel(tau));  stat = zeros(1, T);  pol = zeros(1, T);  tpk = zeros(1, T);
for k = 1:T
    y = R3(k,:).';  r = y;  sel = [];
    for it = 1:P.M
        c = abs(An' * r);  c(sel) = 0;  [~, j] = max(c);  sel(end+1) = j; %#ok<AGROW>
        coef = At(:, sel) \ y;  r = y - At(:, sel)*coef;
    end
    Cmap(k, sel) = coef.';
    [pk, j] = max(abs(coef));  stat(k) = pk;  pol(k) = cos(angle(coef(j)));  tpk(k) = tau(sel(j))*1e9;
end
statS = movmean(stat, 3);                               % lateral smoothing (3 traces)

% ---- baseline statistics: peak |time signal| inside the window ----
s1 = max(abs(td(R1)) .* inw, [], 2).';  s2 = max(abs(td(R2)) .* inw, [], 2).';
outer = abs(x) >= 10;  if ~any(outer), outer = true(1, T); end
con = @(s) 20*log10(max(s) / median(s(outer)));        % peak vs outer traces, dB
[~, kx] = max(statS);

% ---- figure ----
fig = figure('Color', 'w', 'Position', [60 60 1400 780]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact');
title(tl, sprintf('B-scan "%s"%s: %s', label, tag, B.info.desc), 'Interpreter', 'none');
showdB(nexttile, x, tr(ti), X(:, :),  td, ti, rawTitle);
showdB(nexttile, x, tr(ti), R1,       td, ti, sprintf('B1 mean subtraction (contrast %.1f dB)', con(s1)));
showdB(nexttile, x, tr(ti), R2,       td, ti, sprintf('B2 SVD (contrast %.1f dB)', con(s2)));
showdB(nexttile, x, tr(ti), R3,       td, ti, 'SR-SPD residual (ground removed)');

ax = nexttile;  tg = linspace(P.show(1), P.show(2), 400);  img = zeros(numel(tg), T);
for k = 1:T
    for j = find(Cmap(k,:))
        a = abs(Cmap(k,j)) * cos(angle(Cmap(k,j)));
        img(:, k) = img(:, k) + a * exp(-((tg - tau(j)*1e9)/0.08).^2).';
    end
end
imagesc(ax, x, tg, img);  m = max(abs(img(:)));  if m == 0, m = 1; end
caxis(ax, [-m m]);  colormap(ax, rb());  colorbar(ax);
xlabel(ax, 'Position (cm)');  ylabel(ax, 'Time below surface (ns)');
title(ax, 'SR-SPD reflectors (red +, blue -)');
yline(ax, P.win(1), 'g--');  yline(ax, P.win(2), 'g--');

ax = nexttile;  hold(ax, 'on');  grid(ax, 'on');
plot(ax, x, 20*log10(s1/median(s1(outer))), 'o-', 'DisplayName', 'B1');
plot(ax, x, 20*log10(s2/median(s2(outer))), 's-', 'DisplayName', 'B2');
plot(ax, x, 20*log10(statS/median(statS(outer))), 'k^-', 'LineWidth', 1.5, 'DisplayName', 'SR-SPD');
xlabel(ax, 'Position (cm)');  ylabel(ax, 'Statistic vs outer traces (dB)');
title(ax, 'Detection statistic');  legend(ax, 'Location', 'best');
outname = label;  if ~isempty(tag), outname = [label '_bg']; end
saveas(fig, ['bscan_' outname '_analysis.png']);

% ---- summary ----
cS = con(statS);
fprintf('\nB-scan "%s"%s  (%s)\n', label, tag, B.info.desc);
fprintf('  contrast (peak vs outer traces):  B1 %.1f dB | B2 %.1f dB | SR-SPD %.1f dB\n', con(s1), con(s2), cS);
fprintf('  SR-SPD peak at x = %+.1f cm, %.2f ns below surface, polarity %+.2f\n', x(kx), tpk(kx), pol(kx));
fprintf('  drift check at first mark: %.1f dB\n', B.info.drift_dB);
newfile = ~isfile('bscan_results.csv');
fid = fopen('bscan_results.csv', 'a');
if newfile, fprintf(fid, 'label,desc,B1_dB,B2_dB,SRSPD_dB,peak_x_cm,peak_ns,polarity,drift_dB\n'); end
fprintf(fid, '%s,"%s",%.2f,%.2f,%.2f,%.1f,%.2f,%.2f,%.1f\n', [label tag], B.info.desc, con(s1), con(s2), cS, ...
        x(kx), tpk(kx), pol(kx), B.info.drift_dB);
fclose(fid);
fprintf('  figure: bscan_%s_analysis.png   |   results appended to bscan_results.csv\n', outname);
end

% ================================================================
function showdB(ax, x, trs, R, td, ti, ttl)
D = 20*log10(abs(td(R)) + 1e-12);  D = D(:, ti).';
imagesc(ax, x, trs, D);  mx = max(D(:));  caxis(ax, [mx-50 mx]);  colormap(ax, parula);  colorbar(ax);
xlabel(ax, 'Position (cm)');  ylabel(ax, 'Time rel. surface (ns)');  title(ax, ttl);
end

function w = hann1(n)
w = 0.5 - 0.5*cos(2*pi*(0:n-1)/(n-1));
end

function c = rb()
c = [linspace(0.15,1,32)' linspace(0.35,1,32)' ones(32,1); ones(32,1) linspace(1,0.35,32)' linspace(1,0.15,32)'];
end
