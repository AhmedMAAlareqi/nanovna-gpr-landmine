% ANALYZE_SAND  Compare buried targets with the dig controls (no hardware needed)
%   Kernel  : tray lying flat on the sand (radar_tray_sand.mat)
%   Noise   : sand_empty1 (nothing changed after the baseline)
%   For every file: raw difference peak, matched-filter SNR, delay and
%   POLARITY relative to the metal tray (+1 = like metal/water, -1 = opposite:
%   target permittivity LOWER than the sand).

kernel_file = 'radar_tray_sand.mat';
noise_file  = 'radar_sand_empty1.mat';
G = { 'surface',      {'water_surf','oil_surf'}
      'dig control',  {'dig_1','dig_2','dig_3','dig_4'}
      'water buried', {'water_3cm_1','water_3cm_2'}
      'oil buried',   {'oil_3cm_1','oil_3cm_2'}
      'drift (end)',  {'sand_empty2'} };
win_ns = [-1.0 2.0];                 % search window relative to the tray echo (surface .. ~8 cm)

K0 = load(kernel_file);  f = K0.f(:).';  K = K0.Stgt - K0.Sref;
E  = load(noise_file);   De = E.Stgt - E.Sref;
n  = numel(f);  w = 0.5 - 0.5*cos(2*pi*(0:n-1)/(n-1));
dt = (-10:0.01:20)*1e-9;  P = exp(1j*2*pi*f(:)*dt);
nk = sqrt(sum(w.*abs(K).^2));
mfc = @(D) ((conj(K).*D.*w) * P) / nk;               % complex matched-filter output
sigma = sqrt(mean(abs(mfc(De)).^2));
inWin = dt >= win_ns(1)*1e-9 & dt <= win_ns(2)*1e-9;
N = 4096;  t = (0:N-1)/(N*(f(2)-f(1)));  kr = t > 7.5e-9 & t < 14e-9;

fprintf('\n%-14s %-13s %9s %9s %9s %9s\n', 'group', 'file', 'raw dB', 'SNR dB', 'delay ns', 'polarity');
res = struct('g', {}, 'name', {}, 'snr', {}, 'pol', {});
figure('Color', 'w');  tiledlayout(2, 1);
ax1 = nexttile; hold(ax1, 'on'); grid(ax1, 'on');
plot(ax1, dt*1e9, 20*log10(abs(mfc(De))/sigma), 'k:', 'DisplayName', 'sand\_empty1 (noise)');
for gi = 1:size(G, 1)
    for fi = 1:numel(G{gi,2})
        nm = G{gi,2}{fi};  fn = ['radar_' nm '.mat'];
        if ~isfile(fn), fprintf('%-14s %-13s   (not found)\n', G{gi,1}, nm); continue; end
        X = load(fn);  D = X.Stgt - X.Sref;
        raw = max(20*log10(abs(ifft(D.*w, N)) + 1e-12) .* kr + (~kr)*(-300));
        z = mfc(D);  [pk, k] = max(abs(z) .* inWin);
        snr = 20*log10(pk/sigma);  pol = cos(angle(z(k)));
        fprintf('%-14s %-13s %9.1f %9.1f %9.2f %+9.2f\n', G{gi,1}, nm, raw, snr, dt(k)*1e9, pol);
        res(end+1) = struct('g', G{gi,1}, 'name', nm, 'snr', snr, 'pol', pol); %#ok<SAGROW>
        plot(ax1, dt*1e9, 20*log10(abs(z)/sigma), 'LineWidth', 1.2, 'DisplayName', strrep(nm, '_', '\_'));
    end
end
xline(ax1, win_ns(1), 'g--', 'HandleVisibility', 'off');  xline(ax1, win_ns(2), 'g--', 'HandleVisibility', 'off');
xlabel(ax1, 'Delay relative to the tray echo (ns)');  ylabel(ax1, 'MF output / noise (dB)');
legend(ax1, 'Location', 'eastoutside');  ylim(ax1, [-20 60]);  title(ax1, 'Matched filter, sandbox');

% ---- summary: targets vs dig controls ----
dig = [res(strcmp({res.g}, 'dig control')).snr];
fprintf('\nGroup means (SNR dB):\n');
for gi = 1:size(G, 1)
    s = [res(strcmp({res.g}, G{gi,1})).snr];
    if ~isempty(s), fprintf('  %-14s %6.1f  (n = %d)\n', G{gi,1}, mean(s), numel(s)); end
end
if ~isempty(dig)
    thr = max(dig);
    fprintf('\nThreshold = highest dig control = %.1f dB\n', thr);
    for g = {'water buried', 'oil buried'}
        s = [res(strcmp({res.g}, g{1})).snr];
        if ~isempty(s)
            fprintf('  %-13s margin above digs: %+.1f dB (min target - max dig) -> %s\n', g{1}, min(s) - thr, ...
                    ternary(min(s) - thr > 3, 'SEPARATED', 'NOT separated'));
        end
    end
end

ax2 = nexttile;  cols = lines(size(G, 1));
for r = 1:numel(res)
    gi = find(strcmp(G(:,1), res(r).g));
    bar(ax2, r, res(r).snr, 'FaceColor', cols(gi,:)); hold(ax2, 'on');
end
if ~isempty(dig), yline(ax2, max(dig), 'r--', 'max dig control'); end
set(ax2, 'XTick', 1:numel(res), 'XTickLabel', strrep({res.name}, '_', '\_'));  xtickangle(ax2, 45);
ylabel(ax2, 'SNR (dB)');  grid(ax2, 'on');  title(ax2, 'Targets vs dig controls');
saveas(gcf, 'analyze_sand.png');
fprintf('\nPolarity: +1 = same sign as the metal tray (target permittivity HIGHER than sand, e.g. water)\n');
fprintf('          -1 = opposite sign (target permittivity LOWER than sand: plastic-mine case)\n');

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
