% RADAR_TEST3  NanoVNA-H4 + 2 LPDA bistatic radar test, checks against an empty baseline
%   A = Tx on PORT1, B = Rx on PORT2, feed strips 24 cm apart.
%   Type the label WITHOUT "radar_" (e.g. tray30, water, oil, emptyA).
%   After the run it prints the reference and difference peaks and warns
%   if the reference looks contaminated or the file name already exists.

port  = "COM7";   nav = 16;
f1 = 600e6; f2 = 1500e6; npts = 201;
c = 299792458;  t_sys = 6.85e-9;  half = 0.12;  delta = 0.06;
ref_limit = -60.5;      % reference peak (7.5-13 ns) above this = something in the beam

label = strtrim(input('Test label (e.g. tray30, water, oil, emptyA): ', 's'));
label = regexprep(label, '^(radar_)+', '');          % remove "radar_" if typed
if isempty(label), label = 'test'; end
fname = ['radar_' label '.mat'];
if isfile(fname)
    a = input(sprintf('%s already exists. Overwrite? (y/n): ', fname), 's');
    if ~strcmpi(a, 'y'), disp('Stopped. Choose another label.'); return; end
end

clear v
v = serialport(port, 115200, "Timeout", 30);  configureTerminator(v, "CR");
input('1) REFERENCE: target AWAY, YOU behind the antennas. Press Enter: ', 's');
countdown(5); disp('   MEASURING REFERENCE - do not move ...'); beep;
[f, Sref] = sweepS21(v, nav, f1, f2, npts);
beep; disp('   Reference done.');
input('2) NOW place the target, go back BEHIND the antennas, then press Enter: ', 's');
countdown(5); disp('   MEASURING TARGET - do not move ...'); beep;
[~, Stgt] = sweepS21(v, nav, f1, f2, npts);
beep; disp('   Measurement done.');
clear v

n = numel(f);  w = 0.5 - 0.5*cos(2*pi*(0:n-1)/(n-1));  N = 4096;
t = (0:N-1) / (N*(f(2)-f(1)));
H = [ifft(Sref.*w, N); ifft(Stgt.*w, N); ifft((Stgt - Sref).*w, N)];
HdB = 20*log10(abs(H) + 1e-12);

figure('Color', 'w');
plot(t*1e9, HdB, 'LineWidth', 1.3); xlim([0 30]); ylim([-130 -50]); grid on
xlabel('Time (ns)'); ylabel('Amplitude (dB)');
legend('Reference', 'Measurement', 'Difference');
title(['NanoVNA bistatic radar: ' label], 'Interpreter', 'none');

k = t > 7.5e-9 & t < 13e-9;
refPk  = max(HdB(1, k));
[difPk, j] = max(HdB(3, k));  tk = t(k);  te = tk(j);
R = c*(te - t_sys)/2;  d_est = sqrt(max(R^2 - half^2, 0)) - delta;
fprintf('\nCHECK  reference peak %.1f dB | difference peak %.1f dB at %.2f ns | est. distance %.1f cm\n', ...
        refPk, difPk, te*1e9, d_est*100);
if refPk > ref_limit
    fprintf('WARNING: reference peak is high -> something was in the beam during the REFERENCE. Repeat.\n');
end
if isfile('radar_baseline.mat') && ~strcmp(label, 'baseline')
    B = load('radar_baseline.mat');
    dB = 20*log10(abs(ifft((Sref - B.Sref).*w, N)) + 1e-12);
    chg = max(dB(t > 5e-9 & t < 20e-9));
    fprintf('Reference vs baseline change: %.1f dB (should be below about -100)\n', chg);
    if chg > -100
        fprintf('WARNING: the REFERENCE differs from the empty baseline -> something (or someone) was in front. Repeat.\n');
    end
end
fprintf('Guide: tray -60..-65 | water about -95 | empty below about -103 (dB)\n');

save(fname, 'f', 'Sref', 'Stgt');
saveas(gcf, ['radar_' label '.png']);
fprintf('Saved %s in %s\n', fname, pwd);

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

function countdown(n)
for k = n:-1:1
    fprintf('   starting in %d s - stay behind the antennas\n', k); pause(1);
end
end
