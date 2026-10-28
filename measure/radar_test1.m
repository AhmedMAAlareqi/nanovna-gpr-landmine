% RADAR_TEST1  NanoVNA-H4 + 2 LPDA bistatic radar test (tray echo)
%   A = Tx on PORT1 (S11), B = Rx on PORT2 (S21), 24 cm apart.
%   Measures S21 without and with a metal tray, converts to time domain,
%   and prints the delay of the direct coupling and of the tray echo.
%   NanoVNA must show: START 600 MHz, STOP 1.5 GHz, 201p, C1, BW 100 Hz.

port = "COM7";      % NanoVNA serial port
nav  = 4;           % sweeps averaged per measurement
f1 = 600e6; f2 = 1500e6; npts = 201;

clear v
v = serialport(port, 115200, "Timeout", 30);
configureTerminator(v, "CR");

input('1) Tray AWAY (>1.5 m, face down), hands away. Press Enter: ', 's');
tic; [f, Sref] = sweepS21(v, nav, f1, f2, npts);
fprintf('   reference done (%.1f s per sweep)\n', toc/nav);

input('2) Tray UPRIGHT at 30 cm, centred. Step away. Press Enter: ', 's');
[~, Stgt] = sweepS21(v, nav, f1, f2, npts);
fprintf('   target done\n');
clear v

% ---- time domain (Hann window, zero-padded IFFT) ----
n = numel(f);
w = 0.5 - 0.5*cos(2*pi*(0:n-1)/(n-1));
N = 4096;
t = (0:N-1) / (N*(f(2)-f(1)));
H = [ifft(Sref.*w, N); ifft(Stgt.*w, N); ifft((Stgt - Sref).*w, N)];
HdB = 20*log10(abs(H) + 1e-12);

figure('Color', 'w');
plot(t*1e9, HdB, 'LineWidth', 1.3); xlim([0 30]); grid on
xlabel('Time (ns)'); ylabel('Amplitude (dB)');
legend('No target', 'With tray', 'Difference (tray echo)');
title('NanoVNA bistatic radar test 1');

win = t < 30e-9;
[~, k0] = max(abs(H(1, win)));
[~, kd] = max(abs(H(3, win)));
fprintf('Direct coupling peak %.2f ns, tray echo %.2f ns, gap %.2f ns\n', ...
        t(k0)*1e9, t(kd)*1e9, (t(kd) - t(k0))*1e9);
fprintf('Expected gap for 30 cm: about 1.3-1.5 ns\n');

save radar_test1.mat f Sref Stgt
saveas(gcf, 'radar_test1.png');

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
