% ANALYZE_SAND2  Remove surface (ground-bounce) changes, then look below the surface
%   Idea (SR-SPD, contribution C1/C2): the echo of the sand surface has the shape
%   of the tray-on-sand echo K. Digging/smoothing only moves the surface a few mm
%   and changes its strength, so that change lies in the span of
%       K,  (f/fc-1) K,  (f/fc-1)^2 K ...   (complex coefficients: the first
%   order term covers BOTH a small shift (imaginary part) and a tilt (real part))
%   We project that part out of every difference D = Stgt - Sref and then search
%   the remainder for an echo 0.3 - 2.5 ns BELOW the surface.
%
%   nd = number of derivative terms (1 = small height changes, 2 = larger ones).
%   More terms remove more ground but also start to remove shallow targets.

nd          = 1;
kernel_file = 'radar_tray_sand.mat';
noise_file  = 'radar_sand_empty1.mat';
G = { 'surface',      {'water_surf','oil_surf'}
      'dig control',  {'dig_1','dig_2','dig_3','dig_4'}
      'water buried', {'water_3cm_1','water_3cm_2'}
      'oil buried',   {'oil_3cm_1','oil_3cm_2'}
      'drift (end)',  {'sand_empty2'} };
tw = [0.2 1.2];                        % target window: top of an object 0-10 cm deep (ns)

K0 = load(kernel_file);  f = K0.f(:).';  K = K0.Stgt - K0.Sref;
fc = mean(f);  u = f/fc - 1;           % centred, well conditioned
A = K.';
for m = 1:nd, A = [A, (u.^m .* K).']; end %#ok<AGROW>
proj = @(D) D - (A*(A\D.')).';         % remove the surface-change subspace

n = numel(f);  w = 0.5 - 0.5*cos(2*pi*(0:n-1)/(n-1));
dt = (-5:0.01:15)*1e-9;  P = exp(1j*2*pi*f(:)*dt);
nk = sqrt(sum(w.*abs(K).^2));
mfc = @(R) ((conj(K).*R.*w) * P) / nk;
inWin = dt >= tw(1)*1e-9 & dt <= tw(2)*1e-9;

E  = load(noise_file);  Rn = proj(E.Stgt - E.Sref);
sigma = sqrt(mean(abs(mfc(Rn)).^2));

fprintf('\nSurface-change model: K + %d order term(s); target window %.1f-%.1f ns\n', nd, tw);
fprintf('%-14s %-13s %10s %11s %9s %9s %10s\n', 'group', 'file', 'explained', 'resid SNR', 'delay ns', 'polarity', 'global ns');
res = struct('g', {}, 'name', {}, 'snr', {}, 'pol', {}, 'expl', {});
figure('Color', 'w');  tiledlayout(2, 1);
ax1 = nexttile;  hold(ax1, 'on');  grid(ax1, 'on');
for gi = 1:size(G, 1)
    for fi = 1:numel(G{gi,2})
        nm = G{gi,2}{fi};  fn = ['radar_' nm '.mat'];
        if ~isfile(fn), fprintf('%-14s %-13s   (not found)\n', G{gi,1}, nm); continue; end
        X = load(fn);  D = X.Stgt - X.Sref;  R = proj(D);
        expl = 100*(1 - sum(abs(R).^2)/sum(abs(D).^2));
        z = mfc(R);  [pk, k] = max(abs(z) .* inWin);
        snr = 20*log10(pk/sigma);  pol = cos(angle(z(k)));
        [~, kg] = max(abs(z));
        edge = '';  if abs(dt(k)*1e9 - tw(1)) < 0.02 || abs(dt(k)*1e9 - tw(2)) < 0.02, edge = '  (at window edge)'; end
        fprintf('%-14s %-13s %9.1f%% %10.1f %9.2f %+9.2f %10.2f%s\n', G{gi,1}, nm, expl, snr, dt(k)*1e9, pol, dt(kg)*1e9, edge);
        res(end+1) = struct('g', G{gi,1}, 'name', nm, 'snr', snr, 'pol', pol, 'expl', expl); %#ok<SAGROW>
        plot(ax1, dt*1e9, 20*log10(abs(z)/sigma), 'LineWidth', 1.2, 'DisplayName', strrep(nm, '_', '\_'));
    end
end
xline(ax1, tw(1), 'g--', 'HandleVisibility', 'off');  xline(ax1, tw(2), 'g--', 'HandleVisibility', 'off');
xlabel(ax1, 'Delay relative to the surface echo (ns)');  ylabel(ax1, 'Residual MF / noise (dB)');
legend(ax1, 'Location', 'eastoutside');  title(ax1, sprintf('After removing surface changes (nd = %d)', nd));

dig = [res(strcmp({res.g}, 'dig control')).snr];
fprintf('\nGroup means (residual SNR dB):\n');
for gi = 1:size(G, 1)
    s = [res(strcmp({res.g}, G{gi,1})).snr];
    if ~isempty(s), fprintf('  %-14s %6.1f  (n = %d)\n', G{gi,1}, mean(s), numel(s)); end
end
if ~isempty(dig)
    thr = max(dig);  fprintf('\nThreshold = highest dig control = %.1f dB\n', thr);
    for g = {'water buried', 'oil buried'}
        s = [res(strcmp({res.g}, g{1})).snr];
        if ~isempty(s)
            mg = min(s) - thr;
            if mg > 3, verdict = 'SEPARATED'; else, verdict = 'NOT separated'; end
            fprintf('  %-13s margin above digs: %+.1f dB -> %s\n', g{1}, mg, verdict);
        end
    end
end

ax2 = nexttile;  cols = lines(size(G, 1));
for r = 1:numel(res)
    gi = find(strcmp(G(:,1), res(r).g));
    bar(ax2, r, res(r).snr, 'FaceColor', cols(gi,:));  hold(ax2, 'on');
end
if ~isempty(dig), yline(ax2, max(dig), 'r--', 'max dig control'); end
set(ax2, 'XTick', 1:numel(res), 'XTickLabel', strrep({res.name}, '_', '\_'));  xtickangle(ax2, 45);
ylabel(ax2, 'Residual SNR (dB)');  grid(ax2, 'on');  title(ax2, 'Targets vs dig controls after ground removal');
saveas(gcf, sprintf('analyze_sand2_nd%d.png', nd));
fprintf('global ns = where the strongest residual is (if far outside the window: container bottom / multiples)\n');
fprintf('\nexplained = share of the change that is just the surface moving / changing strength\n');
fprintf('polarity  : +1 like metal/water (higher permittivity), -1 lower permittivity than sand\n');
