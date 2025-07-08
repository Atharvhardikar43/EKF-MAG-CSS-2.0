% SIMULATE_SATELLITE Simulates satellite attitude dynamics and EKF-based estimation
clear
% Simulation setup
dt = 0.1;
t_end = 100;
steps = t_end / dt;

epoch_ref = datetime(2000, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
now_utc = datetime('now', 'TimeZone', 'UTC');
utc0 = seconds(now_utc - epoch_ref);

% Initial Conditions
r = [6771e3; 0; 0];
v = [0; 7.67e3; 0];
q = [1; 0; 0; 0];
w = deg2rad([1; 0.1; 0]);
J = [10, 11, 12];
mu = 3.986e14;
q_est = [0.99; 0.1; 0; 0];
q_est = q_est / norm(q_est);
w_est = w;

% Sensor paramters
sigma_gyro = deg2rad(0.04)^2;
sigma_css = 0.004^2;
sigma_mag = 0.01^2;
sigma_q_init = 0.01^2;
noise_css = 0.003;  % Sun sensor noise (approx. 1%)
noise_mag = 0.02;   % Magnetometer noise (5%)

% EKF parameters
P = sigma_q_init * eye(4);
Q = sigma_gyro * eye(4);
R = [sigma_mag * eye(3), zeros(3); zeros(3), sigma_css * eye(3)];

q_true_log = zeros(4, steps);
q_est_log = zeros(4, steps);
t_log = now_utc;
ang_error = zeros(steps, 1);

for k = 1:steps
    t_utc = now_utc + seconds((k - 1) * dt);

    q_true_log(:, k) = q;
    q_est_log(:, k) = q_est;

    dot_product = dot(q, q_est);
    angle_rad = 2 * acos(dot_product);
    ang_error(k) = rad2deg(angle_rad);

    [r, v] = propagate_orbit(r, v, dt, mu, w);
    [q, w] = propagate_attitude(q, w, dt);
    [q_hat_prev, w_est] = propagate_attitude(q_est, w_est, dt);

    [sun_true_body, sun_true_eci] = true_sun_vector(t_utc, q);
    [mag_true_body, mag_true_eci] = true_mag_vector(r, t_utc, q);

    [sun_meas_body, ~, mag_meas_body, ~] = actual_measurements(sun_true_body, sun_true_eci, mag_true_body, mag_true_eci, noise_css, noise_mag);

    [q_est, w_est, P] = ekf_update(q_hat_prev, w_est, P, J, sun_meas_body, mag_meas_body, sun_true_eci, mag_true_eci, Q, R, dt);

    t_log(k) = t_utc;
end

hold on;

assignin('base', 'q_true', q_true_log);
assignin('base', 'q_est', q_est_log);

plot_results(t_log, q_true_log, q_est_log, ang_error);
