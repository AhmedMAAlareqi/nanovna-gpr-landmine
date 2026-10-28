% BSCAN_MEASURE  Step-by-step B-scan of the sandbox (NanoVNA-H4 + 2 LPDA, COM7)
%   Container 48 x 38 x 26 cm, slid FRONT-TO-BACK under the fixed antennas.
%   Positions: -15 ... +15 cm in 2.5 cm steps (13 traces). 0 = target under the antennas.
%   First run: records a kernel (tray flat on the sand) for THIS geometry.
%   Saves bscan_<label>.mat  (f, x_cm, S [traces x freqs], Sd drift trace, info)

port = "COM7";  nav = 8;                    % 8 averages ~ 6.5 s per position
f1 = 600e6; f2 = 1500e6; npts = 201;
x_cm = -15:2.5:15;

if exist('v', 'var'), clear v; end
v = serialport(port, 115200, "Timeout", 30);  configureTerminator(v, "CR");

% ---------- kernel for this geometry (once, or when the antenna height changes) ----------
if ~isfile('bscan_kernel.mat') || strcmpi(strtrim(input('Record a NEW kernel (tray on sand)? (y/n) [n]: ', 's')), 'y')
    input('KERNEL 1/2: pointer on 0, smooth EMPTY sand, you behind the antennas. Enter: ', 's');
    countdown(3);  [f, S0] = sweep_avg(v, 16, f1, f2, npts);
    input('KERNEL 2/2: pizza tray lying FLAT on the sand, centred. Enter: ', 's');
    countdown(3);  [~, S1] = sweep_avg(v, 16, f1, f2, npts);
    K = S1 - S0;  save('bscan_kernel.mat', 'f', 'K', 'S0', 'S1');
    fprintf('Kernel saved. Now REMOVE the tray and prepare the scan.\n\n');
end

% ---------- scan ----------
label = regexprep(strtrim(input('Scan label (e.g. empty, water3, water10, oil3): ', 's')), '^(bscan_)+', '');
if isempty(label), disp('No label, stopped.'); clear v; return; end
desc  = input('Describe target, depth, sand state (free text): ', 's');
fname = ['bscan_' label '.mat'];
if isfile(fname) && ~strcmpi(strtrim(input(sprintf('%s exists. Overwrite? (y/n): ', fname), 's')), 'y')
    disp('Stopped.'); clear v; return;
end

T = numel(x_cm);  S = zeros(T, npts);  k = 1;  t0 = tic;
fprintf('\nSlide the container so the pointer is on each mark. Step back before Enter.\n');
while k <= T
    a = lower(strtrim(input(sprintf('[%2d/%d] pointer on %+5.1f cm | Enter = measure, r = redo previous, q = quit: ', ...
        k, T, x_cm(k)), 's')));
    if strcmp(a, 'q'), break; end
    if strcmp(a, 'r') && k > 1, k = k - 1; continue; end
    countdown(3);
    [f, S(k,:)] = sweep_avg(v, nav, f1, f2, npts);
    h = 20*log10(abs(ifft(S(k,:).*hann1(npts), 4096)) + 1e-12);
    fprintf('        done - strongest echo %.1f dB\n', max(h));
    k = k + 1;
end
if k <= T
    disp('Scan incomplete - saving the measured positions only.');
    S = S(1:k-1, :);  x_cm = x_cm(1:k-1);
end

% ---------- drift / repositioning check ----------
Sd = [];  dr = NaN;
if ~isempty(S) && isempty(strtrim(input('Drift check: pointer back on the FIRST mark, Enter (s = skip): ', 's')))
    countdown(3);  [~, Sd] = sweep_avg(v, nav, f1, f2, npts);
    dr = 20*log10(max(abs(ifft((Sd - S(1,:)).*hann1(npts), 4096))));
    fprintf('Change at the first mark after the scan: %.1f dB (below about -95 dB is good)\n', dr);
end
clear v

info = struct('label', label, 'desc', desc, 'nav', nav, 'time', datetime('now'), ...
              'duration_s', toc(t0), 'drift_dB', dr);
save(fname, 'f', 'x_cm', 'S', 'Sd', 'info');
fprintf('Saved %s  (%d positions, %.1f min). Run  bscan_analyze(''%s'')\n', fname, size(S,1), toc(t0)/60, label);

% ================================================================
function [f, S] = sweep_avg(v, nav, f1, f2, npts)
beep;  S = 0;
for a = 1:nav
    flush(v);
    writeline(v, sprintf("scan %d %d %d 7", round(f1), round(f2), npts));
    txt = "";  t1 = tic;
    while ~contains(txt, "ch>") && toc(t1) < 120
        pause(0.05);
        if v.NumBytesAvailable > 0, txt = txt + char(read(v, v.NumBytesAvailable, "uint8")); end
    end
    L = splitlines(txt);  M = zeros(0, 5);
    for q = 1:numel(L)
        x = sscanf(char(L(q)), '%f').';
        if numel(x) == 5, M(end+1, :) = x; end %#ok<AGROW>
    end
    if size(M, 1) ~= npts, error('Got %d of %d points. Check NanoVNA settings.', size(M,1), npts); end
    S = S + (M(:,4) + 1i*M(:,5)).' / nav;
end
f = M(:,1).';  beep;
end

function countdown(n)
for s = n:-1:1, fprintf('        %d ...\n', s); pause(1); end
end

function w = hann1(n)
w = 0.5 - 0.5*cos(2*pi*(0:n-1)/(n-1));
end
