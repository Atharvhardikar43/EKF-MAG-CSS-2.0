function [q_est, w_t, P_next] = ekf_update(q_hat_next, w_t, P_prev, J, ...
    sun_meas_body, mag_meas_body, ...
    sun_true_eci, mag_true_eci, Q, R, dt)

    % Compute Jacobians
    [F_t, H_t] = calc_Jacobians(q_hat_next, w_t, J, mag_true_eci, sun_true_eci, dt);

    % Time update (prediction)
    P_hat_next = F_t * P_prev * F_t' + Q;

    % Measurement vectors
    z_t = [mag_meas_body; sun_meas_body];

    % Innovation
    y_t = z_t - H_t * q_hat_next;

    % Innovation covariance
    S_t = H_t * P_hat_next * H_t' + R;

    % Kalman gain
    K_t = P_prev * H_t' / S_t;

    % State update
    update = K_t * y_t;
    q_est = q_hat_next + update;
    q_est = q_est / norm(q_est);

    % Covariance update
    P_next = (eye(4) - K_t * H_t) * P_hat_next;
end
