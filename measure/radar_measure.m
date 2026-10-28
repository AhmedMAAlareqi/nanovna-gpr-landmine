% RADAR_MEASURE  Baseline-once protocol for the NanoVNA-H4 + 2 LPDA radar
%   STEP 1 (once per session): label  baseline  -> NOTHING in front, 1 sweep block.
%   STEP 2 (each target)     : label  tray30 / water / oil / emptyB
%                              -> place the target FIRST, then run; only ONE measurement.
%   The difference is always "target scene minus empty baseline", so the
%   order of placing/removing cannot go wrong any more.
%   Files are compatible with analyze_mf (f, Sref = baseline, Stgt = target).

port  = "COM7";   nav = 16;
f1 = 600e6; f2 = 1500e6; npts = 201;
c = 299792458;  t_sys = 7.61e-9;  half = 0.12;  delta = 0.06;

label = strtrim(input('Label (baseline, tray30, water, oil, emptyB ...): ', 's'));
label = regexprep(label, '^(radar_)+', '');
if isempty(label), disp('No label, stopped.'); return; end
isBase = strcmpi(label, 'baseline');
if ~isBase && ~isfile('radar_baseline.mat')
    error('Record the baseline first: run radar_measure with label  baseline');
end
fname = ['radar_' label '.mat'];
if isfile(fname)
    a = input(sprintf('%s already exists. Overwrite? (y/n): ', fname), 's');
    if ~strcmpi(a, 'y'), disp('Stopped.'); return; end
end

if isBase
    input('BASELINE: NOTHING in front of the antennas, you BEHIND them. Press Enter: ', 's');
else
    input(sprintf('TARGET "%s" must be IN PLACE now. Go BEHIND the antennas. Press Enter: ', label), 's');
end
for k = 5:-1:1, fprintf('   starting in %d s\n', k); pause(1); end
beep; disp('   MEASURING - do not move ...');
clear v
v = serialport(port, 115200, "Timeout", 30);  configureTerminator(v, "CR");
[f, S] = sweepS21(v, nav, f1, f2, npts);
clear v
beep; disp('   Done.');

n = numel(f);  w = 0.5 - 0.5*cos(2*pi*(0:n-1)/(n-1));  N = 4096;
t = (0:N-1) / (N*(f(2)-f(1)));
if isBase
    Sref = S; Stgt = S;                          %#ok<NASGU>
    hb = 20*log10(abs(ifft(S.*w, N)) + 1e-12);
    fprintf('\nBASELINE saved. Peak (7.5-13 ns) %.1f dB (empty setup is about -68..-70)\n', ...
            max(hb(t > 7.5e-9 & t < 13e-9)));
    save(fname, 'f', 'Sref', 'Stgt');
    return
end

B = load('radar_baseline.mat');  Sref = B.Sref;  Stgt = S;
H = [ifft(Sref.*w, N); ifft(Stgt.*w, N); ifft((Stgt - Sref).*w, N)];
HdB = 20*log10(abs(H) + 1e-12);

figure('Color', 'w');
plot(t*1e9, HdB, 'LineWidth', 1.3); xlim([0 30]); ylim([-130 -50]); grid on
xlabel('Time (ns)'); ylabel('Amplitude (dB)');
legend('Empty baseline', 'Target scene', 'Difference');
title(['NanoVNA radar vs baseline: ' label], 'Interpreter', 'none');

k = t > 7.5e-9 & t < 13e-9;  tk = t(k);
[difPk, j] = max(HdB(3, k));  te = tk(j);
R = c*(te - t_sys)/2;  d_est = sqrt(max(R^2 - half^2, 0)) - delta;
fprintf('\nCHECK  difference peak %.1f dB at %.2f ns | est. distance %.1f cm\n', difPk, te*1e9, d_est*100);
fprintf('Guide: tray -60..-65 | water about -95 | empty below about -103 (dB)\n');
early = max(HdB(3, t < 6e-9));      % nothing in front can arrive this early
fprintf('Early-time check (0-6 ns): %.1f dB (must be below about -100)\n', early);
fprintf('Early level is %.1f dB below the peak (strong echoes leak ~25 dB early: antenna dispersion)\n', difPk - early);
if early > -100 && early > difPk - 20
    fprintf(['WARNING: the SYSTEM changed since the baseline (cables/antennas touched, NanoVNA\n' ...
             '         restarted, or too much time passed). Record a NEW baseline, then repeat.\n']);
end
save(fname, 'f', 'Sref', 'Stgt');
saveas(gcf, ['radar_' label '.png']);
fprintf('Saved %s\n', fname);

% ---------------------------------------------------------------
function [f, S] = sweepS21(v, nav, f1, f2, npts)
S = 0;
for a = 1:nav
    flush(v);
    writeline(v, sprintf("scan %d %d %d 7", round(f1), round(f2), npts));
    txt = ""; t1 = tic;
    while ~contains(txt, "ch>") && toc(t1) < 120
        pause(0.05);
        if v.NumBytesAvailable > 0
            txt = txt + char(read(v, v.NumBytesAvailable, "uint8"));
        end
    end
    L = splitlines(txt); M = zeros(0, 5);
    for k = 1:numel(L)
        x = sscanf(char(L(k)), '%f').';
        if numel(x) == 5, M(end+1, :) = x; end %#ok<AGROW>
    end
    if size(M, 1) ~= npts
        error('Got %d of %d points. Check the NanoVNA sweep settings.', size(M,1), npts);
    end
    S = S + (M(:,4) + 1i*M(:,5)).' / nav;
end
f = M(:,1).';
end
