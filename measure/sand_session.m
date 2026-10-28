% SAND_SESSION  Guided sandbox session: dry sand, water/oil bottles, dig controls
%   NanoVNA-H4 (COM7) + 2 LPDA pointing DOWN at the sand.
%   Walks you through every step, measures, and saves radar_<label>.mat
%   (f, Sref = empty-sand baseline, Stgt = scene) for analyze_sand.m.
%
%   At each step:  Enter = measure | s = skip | b = new baseline first | q = quit
%   You can restart at any step (it then reuses radar_baseline.mat).

port = "COM7";  nav = 16;
f1 = 600e6; f2 = 1500e6; npts = 201;
c = 299792458;  t_sys = 7.61e-9;  half = 0.12;  delta = 0.06;
rebase_min = 10;                     % suggest a new baseline after this many minutes

S = {
 'baseline',    'Smooth, level sand. NOTHING in or on it. You behind the antennas.'
 'sand_empty1', 'Nothing changed (noise check).'
 'tray_sand',   'Pizza tray lying FLAT on the sand, centred under the antennas.'
 'water_surf',  'Tray away. Water bottle LYING on the sand, left-right, centred.'
 'oil_surf',    'Water bottle away. Oil bottle LYING on the sand, same place.'
 'dig_1',       'Oil bottle away. Dig 12 cm trench at centre, refill with NOTHING, smooth.'
 'water_3cm_1', 'Dig the trench, water bottle inside (left-right), 3 cm sand on top, smooth.'
 'dig_2',       'Bottle out, refill with NOTHING, smooth.'
 'oil_3cm_1',   'Oil bottle in the trench, 3 cm sand on top, smooth.'
 'dig_3',       'Bottle out, refill with NOTHING, smooth.'
 'water_3cm_2', 'Water bottle buried again, 3 cm on top, smooth.'
 'oil_3cm_2',   'Water out. Oil bottle buried, 3 cm on top, smooth.'
 'dig_4',       'Bottle out, refill with NOTHING, smooth.'
 'sand_empty2', 'Nothing changed since dig_4 (drift check).'};
n = size(S, 1);

fprintf('\nSteps:\n'); for k = 1:n, fprintf('  %2d  %s\n', k, S{k,1}); end
k = input(sprintf('Start at step (1-%d) [1]: ', n));
if isempty(k), k = 1; end

if exist('v', 'var'), clear v; end            % release the port from a previous run
v = serialport(port, 115200, "Timeout", 30);  configureTerminator(v, "CR");

if k > 1
    if ~isfile('radar_baseline.mat'), error('No radar_baseline.mat: start at step 1.'); end
    B = load('radar_baseline.mat');  Sref = B.Sref;  f = B.f;
    tBase = tic;  fprintf('Using the existing radar_baseline.mat (age unknown: consider b).\n');
end

while k <= n
    lab = S{k,1};  txt = S{k,2};  isBase = strcmp(lab, 'baseline');
    fprintf('\n[%d/%d]  %s\n   %s\n', k, n, upper(lab), txt);
    if ~isBase && exist('tBase', 'var') && toc(tBase)/60 > rebase_min
        fprintf('   NOTE: baseline is %.0f min old -> type b to record a new one first.\n', toc(tBase)/60);
    end
    a = lower(strtrim(input('   Enter = measure | s = skip | b = new baseline | q = quit: ', 's')));
    if strcmp(a, 'q'), break; end
    if strcmp(a, 's'), k = k + 1; continue; end
    if strcmp(a, 'b')
        input('   Make the sand SMOOTH and EMPTY (targets out), go behind the antennas, Enter: ', 's');
        [f, Sref] = measure(v, nav, f1, f2, npts);  tBase = tic;
        Stgt = Sref;  save('radar_baseline.mat', 'f', 'Sref', 'Stgt');
        fprintf('   New baseline saved. Now set up step %d again.\n', k);
        continue
    end

    fname = ['radar_' lab '.mat'];
    if isfile(fname) && ~isBase
        r = lower(strtrim(input(sprintf('   %s exists. Overwrite? (y/n): ', fname), 's')));
        if ~strcmp(r, 'y'), k = k + 1; continue; end
    end

    if ~isBase && ~exist('Sref', 'var')
        B = load('radar_baseline.mat');  Sref = B.Sref;  tBase = tic;
    end
    [f, X] = measure(v, nav, f1, f2, npts);
    if isBase
        Sref = X;  Stgt = X;  tBase = tic;
        save('radar_baseline.mat', 'f', 'Sref', 'Stgt');  save('radar_sand_base.mat', 'f', 'Sref', 'Stgt');
        w = hann1(numel(f));  h = 20*log10(abs(ifft(X.*w, 4096)) + 1e-12);
        t = (0:4095)/(4096*(f(2)-f(1)));
        fprintf('   BASELINE saved. Strongest echo (ground) %.1f dB at %.2f ns\n', max(h), t(h == max(h))*1e9);
    else
        Stgt = X;  info.label = lab;  info.note = txt;  info.time = datetime('now');
        info.base_age_min = toc(tBase)/60;
        save(fname, 'f', 'Sref', 'Stgt', 'info');
        check(f, Sref, Stgt, t_sys, c, half, delta);
    end
    k = k + 1;
end
clear v
fprintf('\nSession finished. Run analyze_sand to compare targets with the dig controls.\n');

% ================================================================
function [f, S] = measure(v, nav, f1, f2, npts)
for s = 5:-1:1, fprintf('   starting in %d s - stay behind the antennas\n', s); pause(1); end
beep; disp('   MEASURING - do not move ...');
S = 0;
for a = 1:nav
    flush(v);
    writeline(v, sprintf("scan %d %d %d 7", round(f1), round(f2), npts));
    txt = ""; t1 = tic;
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
f = M(:,1).';
beep; disp('   Done.');
end

function check(f, Sref, Stgt, t_sys, c, half, delta)
N = 4096;  w = hann1(numel(f));  t = (0:N-1)/(N*(f(2)-f(1)));
d = 20*log10(abs(ifft((Stgt - Sref).*w, N)) + 1e-12);
k = t > 7.5e-9 & t < 14e-9;  tk = t(k);
[pk, j] = max(d(k));  te = tk(j);
R = c*(te - t_sys)/2;  dist = sqrt(max(R^2 - half^2, 0)) - delta;
early = max(d(t < 6e-9));
fprintf('   CHECK  difference peak %.1f dB at %.2f ns (~%.0f cm) | early %.1f dB\n', pk, te*1e9, dist*100, early);
if early > -100 && early > pk - 20
    fprintf('   WARNING: system changed since baseline -> type b at the next step and repeat this one.\n');
end
end

function w = hann1(n)
w = 0.5 - 0.5*cos(2*pi*(0:n-1)/(n-1));
end
