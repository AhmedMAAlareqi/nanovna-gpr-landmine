% ANALYZE_SAND3  Sparse deconvolution of the sandbox runs (core of SR-SPD)
%   The matched filter cannot separate a target 0.3-1 ns below the surface,
%   because the antennas smear every echo over ~3 ns. Here each difference
%   D = Stgt - Sref is written as
%       D = [surface change: K, (f/fc-1)K]  +  sum_i c_i * K * exp(-j*2*pi*f*tau_i)
%   i.e. a few shifted copies of the measured echo K (tray on sand).
%   The copies are found one by one (greedy / OMP). The result is a short list
%   of "reflectors": delay tau_i, strength |c_i| relative to the metal tray,
%   and polarity cos(angle(c_i)) (+1 like metal/water, -1 lower permittivity).

M  = 6;                                % max number of reflectors per run
tw = [0.2 1.5];                        % target window below the surface (ns)
kernel_file = 'radar_tray_sand.mat';
G = { 'surface',      {'water_surf','oil_surf'}
      'dig control',  {'dig_1','dig_2','dig_3','dig_4'}
      'water buried', {'water_3cm_1','water_3cm_2'}
      'oil buried',   {'oil_3cm_1','oil_3cm_2'}
      'drift (end)',  {'sand_empty2'} };

K0 = load(kernel_file);  f = K0.f(:).';  K = K0.Stgt - K0.Sref;
u  = f/mean(f) - 1;
Gm = [K.', (u.*K).'];                                   % surface-change subspace
tau = (-1.5:0.02:8)*1e-9;  tau(abs(tau) < 0.1e-9) = []; % 0 is covered by Gm
A  = K.' .* exp(-1j*2*pi*f.'*tau);                      % shifted copies of K
An = A ./ vecnorm(A);

fprintf('\n%-14s %-13s %11s %9s %9s %8s\n', 'group', 'file', 'target dB', 'tau ns', 'polarity', 'n refl');
res = struct('g', {}, 'name', {}, 'tgt', {}, 'pol', {});
cm = [linspace(0.15,1,32)' linspace(0.35,1,32)' ones(32,1); ones(32,1) linspace(1,0.35,32)' linspace(1,0.15,32)'];
figure('Color', 'w');  ax = axes;  hold(ax, 'on');  grid(ax, 'on');  names = {};
for gi = 1:size(G, 1)
    for fi = 1:numel(G{gi,2})
        nm = G{gi,2}{fi};  fn = ['radar_' nm '.mat'];
        if ~isfile(fn), fprintf('%-14s %-13s   (not found)\n', G{gi,1}, nm); continue; end
        X = load(fn);  D = (X.Stgt - X.Sref).';  E0 = norm(D)^2;
        r = D - Gm*(Gm\D);  sel = [];  coef = [];
        for it = 1:M
            c = abs(An' * r);  c(sel) = 0;  [~, j] = max(c);  sel(end+1) = j; %#ok<SAGROW>
            B = [Gm, A(:, sel)];  coef = B \ D;  r = D - B*coef;
            if norm(r)^2 < 0.01*E0, break; end
        end
        ca = coef(3:end);  ta = tau(sel);
        inw = (ta >= tw(1)*1e-9 & ta <= tw(2)*1e-9).';
        if any(inw)
            [pk, j] = max(abs(ca) .* inw);  tgt = 20*log10(pk);  tt = ta(j)*1e9;  pol = cos(angle(ca(j)));
        else
            tgt = -Inf;  tt = NaN;  pol = NaN;
        end
        fprintf('%-14s %-13s %11.1f %9.2f %+9.2f %8d\n', G{gi,1}, nm, tgt, tt, pol, numel(sel));
        res(end+1) = struct('g', G{gi,1}, 'name', nm, 'tgt', tgt, 'pol', pol); %#ok<SAGROW>
        names{end+1} = strrep(nm, '_', '\_'); %#ok<SAGROW>
        sz = 15 + 3*max(0, 50 + 20*log10(abs(ca.')));
        scatter(ax, ta*1e9, numel(names)*ones(size(ta)), sz, cos(angle(ca.')), 'filled', 'MarkerEdgeColor', [.3 .3 .3]);
    end
end
colormap(ax, cm);  caxis(ax, [-1 1]);  cb = colorbar(ax);  cb.Label.String = 'polarity (-1 lower \epsilon, +1 like metal/water)';
xline(ax, tw(1), 'g--');  xline(ax, tw(2), 'g--');
set(ax, 'YTick', 1:numel(names), 'YTickLabel', names, 'YDir', 'reverse');
xlabel(ax, 'Delay below the surface echo (ns)');  xlim(ax, [-1.5 8]);
title(ax, 'Deconvolved reflectors (marker size = strength)');
saveas(gcf, 'analyze_sand3.png');

dig = [res(strcmp({res.g}, 'dig control')).tgt];
fprintf('\nStrongest reflector inside %.1f-%.1f ns, relative to the metal tray (dB):\n', tw);
for gi = 1:size(G, 1)
    s = [res(strcmp({res.g}, G{gi,1})).tgt];
    if ~isempty(s), fprintf('  %-14s mean %6.1f   values: %s\n', G{gi,1}, mean(s(isfinite(s))), mat2str(round(s, 1))); end
end
if ~isempty(dig)
    thr = max(dig);  fprintf('\nThreshold = strongest dig-control reflector = %.1f dB\n', thr);
    for g = {'water buried', 'oil buried'}
        s = [res(strcmp({res.g}, g{1})).tgt];
        if ~isempty(s)
            mg = min(s) - thr;
            if mg > 3, verdict = 'SEPARATED'; else, verdict = 'NOT separated'; end
            fprintf('  %-13s margin: %+.1f dB -> %s\n', g{1}, mg, verdict);
        end
    end
end
