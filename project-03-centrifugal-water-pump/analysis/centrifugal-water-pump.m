%% Centrifugal Pump Design & Performance Sizing (Imperial Units)
% Project: Centrifugal Water Pump
%
% Purpose:
%   Size an impeller and estimate inlet/outlet velocity triangles for a
%   centrifugal pump meeting a target operating point, then sanity-check
%   blade solidity and diffusion factor. Also includes a simple iterative
%   volute sizing loop.
%
% Notes:
%   - This script uses Imperial units throughout unless otherwise noted.
%   - Several coefficients are course/figure-based (Cordier diagram, etc.).
%   - "overrideBladeCount" lets you force blade count to satisfy solidity.

clc; clear;

%% ===========================
%  1) INPUTS (EDIT HERE)
%  ===========================

% Target operating point
flow_gpm       = 1500;   % [gpm]
speed_rpm      = 2500;   % [rpm]
head_ft        = 500;    % [ft]

% Material / sizing assumptions
allowableShear_psi = 10000;  % [psi] torsional shear stress limit
r1_in              = 3.2;    % [in] impeller inlet radius (eye radius)

% Blade angles (metal angles)
betaB1_deg = 17.0;            % [deg] inlet blade angle
betaB2_deg = 22.5;            % [deg] outlet blade angle (target band ~22.5–27.5)

% Cordier diagram selection (non-dimensional specific diameter)
Ds          = 8.2;            % [-] from Cordier diagram for selected family
etaHyd_est  = 0.78;           % [-] estimated hydraulic efficiency (Cordier)

% Controls
overrideBladeCount = 9;       % [-] set to 0 to use rounded computed blade count
bladeProfileSteps  = 10;      % [-] discretization steps from r1->r2 for blade profile

%% ===========================
%  2) CONSTANTS & CONVERSIONS
%  ===========================

g_ft_s2     = 32.2;            % [ft/s^2]
rhoWater_lb_ft3 = 62.4;        % [lb/ft^3] (weight density in Imperial)

inPerFt     = 12;
ft2PerIn2   = 1/144;

% Common design coefficients (course-based)
phiEye      = 0.27;            % [-] eye flow coefficient used in r_eye calc
velMargin   = 1.10;            % [-] V1 = 1.1 * Ve (inlet velocity margin)
epsilon2Nom = 0.90;            % [-] general blockage factor (legacy; not used directly)

% Solidity and diffusion acceptability bands
solidityMin = 2.5;
solidityMax = 3.0;
diffMin     = 0.7;
diffMax     = 1.0;

%% ===========================
%  3) BASIC DERIVED QUANTITIES
%  ===========================

flow_ft3_s = flow_gpm * (1/60) * 231 / (inPerFt^3);   % [ft^3/s]
omega_rad_s = 2*pi*speed_rpm/60;                     % [rad/s]

% Course: specific speed definitions
Ns_dim = (speed_rpm * sqrt(flow_gpm)) / (head_ft^(3/4));  % (dimensional form used in class)
omega_s = (omega_rad_s * sqrt(flow_ft3_s)) / ((g_ft_s2 * head_ft)^(3/4));  % [-]

%% ===========================
%  4) SHAFT POWER ESTIMATE
%  ===========================
% Pressure rise (as lb/ft^2) from head:
%   Δp = γ * H ; where γ ≈ 62.4 lb/ft^3
deltaP_lb_ft2 = rhoWater_lb_ft3 * head_ft;

% Shaft power from hydraulic power / efficiency:
%   P = Q*Δp/η
shaftPower_ftlbf_s = (flow_ft3_s * deltaP_lb_ft2) / etaHyd_est; % [ft·lbf/s]

% Torque:
%   τ = P/ω
torque_ftlbf = shaftPower_ftlbf_s / omega_rad_s;               % [ft·lbf]

%% ===========================
%  5) INLET (EYE) SIZING
%  ===========================

% Shaft diameter from torsion:
%   τ = T*r/J -> D = (16*T/(pi*τ_allow))^(1/3)
% Convert torque to [in·lbf] for psi usage:
torque_inlbf = torque_ftlbf * inPerFt;

shaftDiam_in = (16 * torque_inlbf / (pi * allowableShear_psi))^(1/3); % [in]

% Eye radius based on coefficient relationship from class (phiEye)
rEye_ft = (flow_ft3_s / (pi * phiEye * omega_rad_s))^(1/3);   % [ft]
rEye_in = rEye_ft * inPerFt;

dEye_ft = 2 * rEye_ft;
dEye_in = 2 * rEye_in;

% Mean inlet velocity at eye
Ve_ft_s = flow_ft3_s / (pi * (dEye_ft^2) / 4);     % [ft/s]
V1_ft_s = velMargin * Ve_ft_s;                     % [ft/s] design inlet velocity

% Inlet geometry at r1 (given)
d1_in = 2 * r1_in;
r1_ft = r1_in / inPerFt;
d1_ft = d1_in / inPerFt;

A1_ft2 = flow_ft3_s / V1_ft_s;                     % [ft^2]
A1_in2 = A1_ft2 / ft2PerIn2;                       % [in^2]

b1_ft = A1_ft2 / (2*pi*r1_ft);                     % [ft]
b1_in = b1_ft * inPerFt;                           % [in]

% Inlet velocity triangle
U1_ft_s = omega_rad_s * r1_ft;                     % [ft/s]

betaF1_rad = atan(V1_ft_s / U1_ft_s);              % [rad]
betaF1_deg = rad2deg(betaF1_rad);                  % [deg]
W1_ft_s = hypot(U1_ft_s, V1_ft_s);                 % [ft/s] relative magnitude

%% ===========================
%  6) OUTLET SIZING (CORDIER)
%  ===========================

% Sanity check: outlet blade angle band (course)
if betaB2_deg >= 22.5 && betaB2_deg <= 27.5
    outletAngleMsg = "Outlet blade angle OK (22.5–27.5 deg).";
else
    outletAngleMsg = "Outlet blade angle OUT OF RANGE; target 22.5–27.5 deg.";
end

% Cordier-based impeller diameter
d2_ft = (Ds * sqrt(flow_ft3_s)) / ((g_ft_s2 * head_ft)^(1/4)); % [ft]
d2_in = d2_ft * inPerFt;

r2_ft = d2_ft / 2;
r2_in = d2_in / 2;

U2_ft_s = omega_rad_s * r2_ft;                                   % [ft/s]

% Blade count estimate (course heuristic)
nu = d1_ft / d2_ft;
Zb_est = 6.5 * ((1 + nu)/(1 - nu)) * sind((betaB1_deg + betaB2_deg)/2);
Zb_round = round(Zb_est);

% Head coefficient / hydraulic efficiency correlation (as used in your original)
etaHyd = 1 - (0.80 / (flow_gpm^(0.25)));

% Slip factor approximation
muSlip = 1 - pi * sind(betaB2_deg) / Zb_round;

% Euler: V_u2 = gH/(η_h U2)
Vu2_ft_s = (g_ft_s2 * head_ft) / (etaHyd * U2_ft_s);

% Radial component from outlet triangle using blade angle and slip
Vr2_ft_s = ((Vu2_ft_s / muSlip) - U2_ft_s) * -tand(betaB2_deg);

W2_ft_s = hypot(Vr2_ft_s, (U2_ft_s - Vu2_ft_s));
diffusionFactor = W2_ft_s / W1_ft_s;

A2_ft2 = flow_ft3_s / Vr2_ft_s;
A2_in2 = A2_ft2 / ft2PerIn2;

% Passage blockage factor at outlet
epsilon2 = 1 - Zb_round * 0.25 / (pi * d2_in * sind(betaB2_deg));
b2_in = A2_in2 / (2 * epsilon2 * pi * (d2_in/2));              % [in]

%% ===========================
%  7) BLADE PROFILE (r-theta)
%  ===========================
% Discretize from r1 -> r2, linearly varying blade angle from betaB1 to betaB2.

dr_in = (r2_in - r1_in) / bladeProfileSteps;

idx = 1;
theta_deg = 0;
arcLen_in = 0;

r_in   = zeros(1, bladeProfileSteps+1);
betaB_deg = zeros(1, bladeProfileSteps+1);
dTheta_deg = zeros(1, bladeProfileSteps+1);
thetaVec_deg = zeros(1, bladeProfileSteps+1);
L_in   = zeros(1, bladeProfileSteps+1);

for rStep = r1_in : dr_in : r2_in
    r_in(idx) = rStep;

    % Linear interpolation of blade angle from inlet to outlet
    betaB_deg(idx) = betaB1_deg + (betaB2_deg - betaB1_deg) * ((rStep - r1_in) / (r2_in - r1_in));

    % From geometry: dθ = dr/(r*tan(beta))  [radians], convert to degrees
    dTheta_deg(idx) = rad2deg( dr_in / (r_in(idx) * tand(betaB_deg(idx))) );

    thetaVec_deg(idx) = theta_deg;
    theta_deg = theta_deg + dTheta_deg(idx);

    % Incremental blade length
    dL_in = hypot(r_in(idx) * deg2rad(dTheta_deg(idx)), dr_in);
    L_in(idx) = arcLen_in;
    arcLen_in = arcLen_in + dL_in;

    idx = idx + 1;
end

maxBladeLen_in = max(L_in);

%% ===========================
%  8) SOLIDITY CHECK & RE-CALC
%  ===========================
% Solidity (as used in your original):
%   σ = Z * (L / (π D2))

if overrideBladeCount ~= 0
    bladeCount = overrideBladeCount;
else
    bladeCount = Zb_round;
end

solidity = bladeCount * (maxBladeLen_in / (pi * d2_in));

solidityOK = (solidity >= solidityMin) && (solidity <= solidityMax);
diffusionOK = (diffusionFactor >= diffMin) && (diffusionFactor <= diffMax);

% If we override blade count, re-compute slip and dependent outlet terms using max(betaB)
if overrideBladeCount ~= 0
    betaB2_eff_deg = max(betaB_deg);

    muSlip = 1 - pi * sind(betaB2_eff_deg) / bladeCount;
    Vr2_ft_s = ((Vu2_ft_s / muSlip) - U2_ft_s) * -tand(betaB2_eff_deg);

    A2_ft2 = flow_ft3_s / Vr2_ft_s;
    A2_in2 = A2_ft2 / ft2PerIn2;

    epsilon2 = 1 - bladeCount * 0.25 / (pi * d2_in * sind(betaB2_eff_deg));
    b2_in = A2_in2 / (2 * epsilon2 * pi * (d2_in/2));

    W2_ft_s = hypot(Vr2_ft_s, (U2_ft_s - Vu2_ft_s));
    diffusionFactor = W2_ft_s / W1_ft_s;

    diffusionOK = (diffusionFactor >= diffMin) && (diffusionFactor <= diffMax);
end

%% ===========================
%  9) DIFFUSER / VOLUTE SIZING
%  ===========================
% This section preserves your original approach:
% - Angular momentum at r2
% - Use K3 from a figure for throat velocity
% - Iterate volute centroid radius r_c until convergence

angMom_ft2_s = r2_ft * Vu2_ft_s;            % [ft^2/s]
angMom_in2_s = angMom_ft2_s * inPerFt;      % [in^2/s] legacy scale as in original

K3 = 0.46;                                  % from course figure (example)
Vt_ft_s = K3 * sqrt(2 * g_ft_s2 * head_ft); % [ft/s]
Vt_in_s = Vt_ft_s * inPerFt;

At_ft2 = flow_ft3_s / Vt_ft_s;
At_in2 = At_ft2 / ft2PerIn2;

% d3 set as +8% over d2 per your original
d3_in = 1.08 * d2_in;
r3_in = d3_in / 2;

Ccoef = 0.95;
bVolute_in = 1.8 * b2_in;
phiTrap_deg = 30; %#ok<NASGU> (kept for traceability)

Vtheta_in_s = Vt_in_s;
rc_prev_in = 0;
rc_in = (Ccoef * angMom_in2_s) / Vtheta_in_s;

flow_in3_s = flow_ft3_s * (inPerFt^3);

% Convergence loop (kept logically consistent with your original)
while abs((rc_in - rc_prev_in) / Vtheta_in_s) >= 0.01

    Atheta_in2 = (flow_in3_s * max(thetaVec_deg)) / (2 * pi * Vtheta_in_s);

    % Solve for trapezoid dimensions (original algebra preserved)
    a_in = sqrt((2 * Atheta_in2 * (2 * tand(15))) + bVolute_in^2);
    h_in = (2 * Atheta_in2) / (a_in + bVolute_in);

    xc_in = (h_in * ((2 * a_in) + bVolute_in)) / (3 * (a_in + bVolute_in));

    rc_prev_in = rc_in;
    rc_in = r3_in + xc_in;

    % Update tangential velocity estimate
    Vtheta_in_s = (Ccoef * angMom_in2_s) / rc_in;
end

%% ===========================
%  10) PRINT SUMMARY
%  ===========================

hp = shaftPower_ftlbf_s / 550;

fprintf("\n=== Centrifugal Pump Sizing Summary ===\n");
fprintf("Target: Q = %.0f gpm, H = %.0f ft, N = %.0f rpm\n", flow_gpm, head_ft, speed_rpm);
fprintf("Specific speed (Ns): %.2f | omega_s: %.4f\n", Ns_dim, omega_s);
fprintf("%s\n\n", outletAngleMsg);

fprintf("Shaft power: %.0f ft·lbf/s (%.1f HP)\n", shaftPower_ftlbf_s, hp);
fprintf("Torque: %.0f ft·lbf | Shaft diameter: %.2f in\n\n", torque_ftlbf, shaftDiam_in);

fprintf("Inlet:\n");
fprintf("  Eye diameter: %.2f in | D1: %.2f in | b1: %.2f in\n", dEye_in, d1_in, b1_in);
fprintf("  V1: %.2f ft/s | U1: %.2f ft/s | W1: %.2f ft/s | beta_f1: %.2f deg\n\n", ...
    V1_ft_s, U1_ft_s, W1_ft_s, betaF1_deg);

fprintf("Outlet:\n");
fprintf("  D2: %.2f in | U2: %.2f ft/s\n", d2_in, U2_ft_s);
fprintf("  Blade count: %d (est %.2f -> round %d)\n", bladeCount, Zb_est, Zb_round);
fprintf("  Vu2: %.2f ft/s | Vr2: %.2f ft/s | b2: %.2f in\n", Vu2_ft_s, Vr2_ft_s, b2_in);
fprintf("  Diffusion factor: %.3f (%s)\n", diffusionFactor, ternary(diffusionOK,"OK","OUT OF RANGE"));
fprintf("  Solidity: %.3f (%s)\n\n", solidity, ternary(solidityOK,"OK","OUT OF RANGE"));

fprintf("Volute iteration result:\n");
fprintf("  Throat velocity Vt: %.2f ft/s | Throat area At: %.2f in^2\n", Vt_ft_s, At_in2);
fprintf("  Volute centroid radius rc: %.2f in\n", rc_in);
fprintf("======================================\n");

%% ===========================
%  Local helper (MATLAB allows local functions at end of script)
%  ===========================
function out = ternary(cond, a, b)
%TERNARY Simple inline conditional string selection
    if cond, out = a; else, out = b; end
end
