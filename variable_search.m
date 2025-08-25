%% === README ===
%{ 
    This MATLAB code searches the optimal component values (capacitance, resistance, frequency)
   that gives us the maximum phase/gain difference between two different positions: Metal
   wire fully inserted and not inserted. 

    4 different plots are drawn as output of this code. 2 of them show the
   change in gain and phase with respect to the inductance, with the optimal
   component values. (Optimal component values are printed on the terminal)
   Other two plots again show the gain and phase change with respect to the 
   inductance, however for multiple combinations of the component values.

     User should only change the variables in the "VARIABLE RANGE SETUP"
    section, according to the explanations.

%}

clc; clear; close all;
%% === VARIABLE RANGE SETUP ===

Rin_val = 50;   %Internal resistance of the signal generator (DO NOT CHANGE THIS VALUE) 
Vin_val = 1;    %Amplitude of the sinusoidal signal obtained from the signal generator
Rl_val  = 15;   %Series resistance of the solenoid (UNKNOWN VALUE, NEED TO MEASURE)

L_measured_min = 3000e-6;     %Minimum value of the measured solenoid inductance
L_measured_max = 3400e-6;     %Maximum value of the measured solenoid inductance

% L sweep for x-axis (Henry)
L_val_min = 2500e-6;  %Minimum value of the inductance range
L_val_max = 4000e-6;  %Maximum value of the inductance range
L_vals  = linspace(L_val_min, L_val_max, 10000);    %Inductance range setup

% Frequency range (rad/s)
w_val_min = 2*pi*100e3;     %Minimum value of the frequency range
w_val_max = 2*pi*1e6;   %Maximum value of the frequency range
w_vals  = linspace(w_val_min, w_val_max, 500);   %Frequency range setup

% Capacitance range (Farad)
C_val_min = 0.001e-9;     %Minimum value of the capacitance range
C_val_max = 10e-9;      %Maximum value of the capacitance range
C_vals   = linspace(C_val_min, C_val_max, 1000);     %Capacitance range setup

% Resistance range (Ohm)
R_val_min = 220;    %Minimum value of the resistance range
R_val_max = 10000;   %Maximum value of the resistance ragne
Rout_vals = linspace(R_val_min, R_val_max, 10);     %Resistance range setup


%% === SYMBOLIC SETUP ===
syms w Vin Rin Rout Rl C1 C2 L real
syms Vout Va Vb Vc

A = [ 0,                   (1/Rin)+1i*w*C1,               -1i*w*C1,                             0;
     -1i*w*C2,            -1i*w*C1,                       (1/(1i*w*L))+1i*w*(C1+C2),            -1/(1i*w*L);
      0,                   0,                              -1/(1i*w*L),                         (1/Rl)+(1/(1i*w*L));
     (1/Rout)+1i*w*C2,     0,                              -1i*w*C2,                            0];

b = [Vin/Rin; 0; 0; 0];

n = length(b);
x = sym('x', [n, 1]);
detA = det(A);
for i = 1:n
    Ai = A; Ai(:, i) = b;
    x(i) = simplify(det(Ai) / detA);
end

abs_Vout   = simplify(sqrt(x(1)*conj(x(1))));
abs_Va     = simplify(sqrt(x(2)*conj(x(2))));
phase_Vout = simplify((180/pi)*angle(x(1)));
phase_Va   = simplify((180/pi)*angle(x(2)));

gain_dB_Vout_expr = 20*log10(abs_Vout/abs_Va);

% Function handles (respect this argument order everywhere)
% [L, C1, C2, Rout, Rin, Rl, Vin, w]
gain_dB_func   = matlabFunction(real(gain_dB_Vout_expr), 'Vars', [L, C1, C2, Rout, Rin, Rl, Vin, w]);
phase_Vout_fun = matlabFunction(phase_Vout,         'Vars', [L, C1, C2, Rout, Rin, Rl, Vin, w]);
phase_Va_fun   = matlabFunction(phase_Va,           'Vars', [L, C1, C2, Rout, Rin, Rl, Vin, w]);

%% === TRACKING BEST CASES ===
max_phase_change = -inf;
max_gain_change  = -inf;

best_params_phase = struct('C1',NaN,'C2',NaN,'Rout',NaN,'freq',NaN,'raw',NaN,'unwrapped',NaN);
best_params_gain  = struct('C1',NaN,'C2',NaN,'Rout',NaN,'freq',NaN);

%% === PLOT: PHASE SEARCH ===
fig1 = figure('Name','Phase of Vout vs L for Different C,R and w');
ax1  = axes('Parent',fig1); hold(ax1,'on'); grid(ax1,'on');
xlabel(ax1,'L (\muH)'); ylabel(ax1,'Phase of V_{out} (deg)');
title(ax1,'Phase of V_{out} vs L (Unwrapped)');

for w_val = w_vals
    for c2 = C_vals
        c1 = c2;
        for r_out = Rout_vals

            % Raw wrapped phases at two L values:
            phi1_raw_Vout = phase_Vout_fun(L_measured_min, c1, c2, r_out, Rin_val, Rl_val, Vin_val, w_val);
            phi2_raw_Vout = phase_Vout_fun(L_measured_max, c1, c2, r_out, Rin_val, Rl_val, Vin_val, w_val);

            phi1_raw_Va   = phase_Va_fun(L_measured_min, c1, c2, r_out, Rin_val, Rl_val, Vin_val, w_val);
            phi2_raw_Va   = phase_Va_fun(L_measured_max, c1, c2, r_out, Rin_val, Rl_val, Vin_val, w_val);

            diff1_raw = abs(phi1_raw_Vout - phi1_raw_Va);
            diff2_raw = abs(phi2_raw_Vout - phi2_raw_Va);
            delta_raw = abs(diff2_raw - diff1_raw);

            % Unwrap for robustness
            phi1_rad_Vout = deg2rad(phi1_raw_Vout);
            phi2_rad_Vout = deg2rad(phi2_raw_Vout);
            phi1_rad_Va   = deg2rad(phi1_raw_Va);
            phi2_rad_Va   = deg2rad(phi2_raw_Va);

            phi_unwrapped1  = unwrap([phi1_rad_Vout, phi1_rad_Va]);
            phi_unwrapped2  = unwrap([phi2_rad_Vout,   phi2_rad_Va]);
            diff1_unwrapped = abs(phi_unwrapped1(1) - phi_unwrapped1(2));
            diff2_unwrapped = abs(phi_unwrapped2(1) - phi_unwrapped2(2));
            delta_unwrapped = abs(diff2_unwrapped - diff1_unwrapped);

            % Reject obvious wrap artifacts
            if abs(delta_unwrapped - deg2rad(delta_raw)) > pi
                continue;
            end

            % Update max and plot this curve immediately
            if delta_unwrapped > max_phase_change
                max_phase_change = delta_unwrapped;
                best_params_phase.C1   = c1;
                best_params_phase.C2   = c2;
                best_params_phase.Rout = r_out;
                best_params_phase.freq = w_val/(2*pi);
                best_params_phase.raw  = delta_raw;
                best_params_phase.unwrapped = delta_unwrapped;

                y_vals = arrayfun(@(LL) phase_Vout_fun(LL, c1, c2, r_out, Rin_val, Rl_val, Vin_val, w_val), L_vals);
                h = plot(ax1, L_vals*1e6, y_vals);
                h.UserData = struct('C1',c1,'C2',c2,'Rout',r_out,'freq',w_val/(2*pi), 'metric','Phase Vout');
            end
        end
    end
end

fprintf('Maximum phase change: %.2f deg (unwrapped)\n', rad2deg(max_phase_change));
fprintf('Best Phase Params:\n  C1 = %.2e F\n  C2 = %.2e F\n  Rout = %.0f Ohm\n  f = %.0f Hz\n', ...
    best_params_phase.C1, best_params_phase.C2, best_params_phase.Rout, best_params_phase.freq);

% Best Phase comparison
fig2 = figure('Name','Best Phase (Vout vs Va)');
ax2  = axes('Parent',fig2); hold(ax2,'on'); grid(ax2,'on');
xlabel(ax2,'L (\muH)'); ylabel(ax2,'Phase (deg)'); title(ax2,'Best Phase');

y_vout = arrayfun(@(LL) phase_Vout_fun(LL, best_params_phase.C1, best_params_phase.C2, ...
                                       best_params_phase.Rout, Rin_val, Rl_val, Vin_val, 2*pi*best_params_phase.freq), L_vals);
plot(ax2, L_vals*1e6, y_vout, 'LineWidth', 1.2);

y_va = arrayfun(@(LL) phase_Va_fun(LL, best_params_phase.C1, best_params_phase.C2, ...
                                   best_params_phase.Rout, Rin_val, Rl_val, Vin_val, 2*pi*best_params_phase.freq), L_vals);
plot(ax2, L_vals*1e6, y_va, 'LineWidth', 1.2);
legend(ax2, 'Phase Vout', 'Phase Va', 'Location','best');

%% === PLOT: GAIN SEARCH ===
fig3 = figure('Name','|Vout| vs L (Gain in dB)');
ax3  = axes('Parent',fig3); hold(ax3,'on'); grid(ax3,'on');
xlabel(ax3,'L (\muH)'); ylabel(ax3,'|V_{out}/V_{a}| (dB)');
title(ax3,'|V_{out}| vs L for different C,R and w');

for w_val = w_vals
    for c2 = C_vals
        c1 = c2;
        for r_out = Rout_vals
            g1 = gain_dB_func(L_measured_min, c1, c2, r_out, Rin_val, Rl_val, Vin_val, w_val);
            g2 = gain_dB_func(L_measured_max, c1, c2, r_out, Rin_val, Rl_val, Vin_val, w_val);

            if (abs(g2 - g1) > max_gain_change) && g1 > -30 && g2 > -30
                max_gain_change = abs(g2 - g1);
                best_params_gain.C1   = c1;
                best_params_gain.C2   = c2;
                best_params_gain.Rout = r_out;
                best_params_gain.freq = w_val/(2*pi);

                y_vals = arrayfun(@(LL) gain_dB_func(LL, c1, c2, r_out, Rin_val, Rl_val, Vin_val, w_val), L_vals);
                plot(ax3, L_vals*1e6, real(y_vals), 'LineWidth', 0.8); % plot inside loop
            end
        end
    end
end
    
fprintf('Maximum gain change: %.2f dB\n', max_gain_change);
fprintf('Best Gain Params:\n  C1 = %.2e F\n  C2 = %.2e F\n  Rout = %.0f Ohm\n  f = %.0f Hz\n', ...
    best_params_gain.C1, best_params_gain.C2, best_params_gain.Rout, best_params_gain.freq);

% Best Gain curve
fig4 = figure('Name','Best Gain');
ax4  = axes('Parent',fig4); hold(ax4,'on'); grid(ax4,'on');
xlabel(ax4,'L (\muH)'); ylabel(ax4,'|V_{out}/V_{a}| (dB)'); title(ax4,'Best Gain');

y_best_gain = arrayfun(@(LL) gain_dB_func(LL, best_params_gain.C1, best_params_gain.C2, ...
                                          best_params_gain.Rout, Rin_val, Rl_val, Vin_val, 2*pi*best_params_gain.freq), L_vals);
plot(ax4, L_vals*1e6, y_best_gain, 'LineWidth', 1.2);
legend(ax4, 'Best Gain (dB)', 'Location','best');
