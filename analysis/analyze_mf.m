% ANALYZE_MF  Matched-filter analysis of saved radar tests (no hardware needed)
%   Uses the measured metal-tray echo as the system kernel (antenna dispersion
%   + NanoVNA response), searches each test for that echo shape, and compares
%   with an empty run to give an SNR. Same idea as SR-SPD contribution C1.
%
%   HOW TO USE
%   1) Put this file in the same folder as the radar_*.mat files.
%   2) Check the names in SETTINGS below (type  dir *.mat  to see your files).
%   3) Type  analyze_mf  in the Command Window.

%% ---------------- SETTINGS ----------------
tray_file  = 'radar_test1.mat';      % metal tray echo (kernel), tray at 40 cm
tray_cm    = 40;                     % distance of the tray in tray_file
target_cm  = 30;                     % distance where the test targets were placed
empty_file = 'radar_EMPTY.mat';      % target away in both sweeps (noise level)
test_files = {'radar_water_h.mat', 'radar_oil.mat', 'radar_empty2.mat', ...
              'radar_oil_h.mat', 'radar_oil_h2.mat', 'radar_oil_h3.mat'};
%% ------------------------------------------

fprintf('\n.mat files in this folder:\n'); disp(ls('*.mat'));
must = {tray_file, empty_file};
for i = 1:2
    if ~isfile(must{i}), error('File "%s" not found. Edit SETTINGS.', must{i}); end
end

T = load(tray_file);   f = T.f(:).';   K = T.Stgt - T.Sref;
E = load(empty_file);  De = E.Stgt - E.Sref;

n  = numel(f);
w  = 0.5 - 0.5*cos(2*pi*(0:n-1)/(n-1));
dt = (-10:0.01:20)*1e-9;                  % delay relative to the tray echo
P  = exp(1j*2*pi*f(:)*dt);
nk = sqrt(sum(w.*abs(K).^2));
mf = @(D) abs((conj(K).*D.*w) * P) / nk;  % matched-filter output vs delay

zE    = mf(De);
sigma = sqrt(mean(zE.^2));                % noise/drift level

% expected delay of a target at target_cm relative to the tray at tray_cm
c = 299792458; half = 0.12; delta = 0.06;
rt = @(d) 2*sqrt((d/100 + delta)^2 + half^2)/c;
d0 = rt(target_cm) - rt(tray_cm);
inWin  = dt > d0 - 1.0e-9 & dt < d0 + 1.5e-9;   % search window around the target
outWin = dt > d0 + 3.5e-9;                       % far away: false-alarm check

figure('Color', 'w'); hold on; grid on
plot(dt*1e9, 20*log10(zE/sigma), 'k:', 'LineWidth', 1, 'DisplayName', 'empty (noise)');
fprintf('\nExpected target delay: %.2f ns relative to the tray echo\n', d0*1e9);
fprintf('%-22s %10s %12s %14s\n', 'file', 'SNR (dB)', 'delay (ns)', 'max outside');
for i = 1:numel(test_files)
    fn = test_files{i};
    if ~isfile(fn), fprintf('%-22s  (not found, skipped)\n', fn); continue; end
    X = load(fn);  z = mf(X.Stgt - X.Sref);
    [pk, k] = max(z .* inWin);
    snr  = 20*log10(pk/sigma);
    snrO = 20*log10(max(z(outWin))/sigma);
    fprintf('%-22s %10.1f %12.2f %14.1f\n', fn, snr, dt(k)*1e9, snrO);
    plot(dt*1e9, 20*log10(z/sigma), 'LineWidth', 1.3, 'DisplayName', strrep(fn, '_', '\_'));
end
xline((d0-1.0e-9)*1e9, 'g--', 'HandleVisibility', 'off');
xline((d0+1.5e-9)*1e9, 'g--', 'HandleVisibility', 'off');
xlabel('Delay relative to the tray echo (ns)');
ylabel('Matched-filter output / noise (dB)');
title('Matched filter with measured tray kernel'); legend('Location', 'best');
ylim([-20 60]);
fprintf('\nSNR > 10 dB: confident | 6-10 dB: marginal | < 6 dB: not detected\n');
fprintf('"max outside" should stay below ~6 dB (false-alarm check)\n');
saveas(gcf, 'analyze_mf.png');
