function [F_t, H_t] = calc_Jacobians(q_hat_next, w_t, J, b, s, dt)

% Dynamics Jacobian F
wx = w_t(1); wy = w_t(2); wz = w_t(3);

F_t = [1, -0.5*dt*wx, -0.5*dt*wy, -0.5*dt*wz;
       0.5*dt*wx, 1, 0.5*dt*wz, -0.5*dt*wy;
       0.5*dt*wy, -0.5*dt*wz, 1, 0.5*dt*wx;
       0.5*dt*wz, 0.5*dt*wy, -0.5*dt*wx, 1];

% Normalize inertial vectors
b = b / norm(b);
s = s / norm(s);

q = q_hat_next;
q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);

% Helper function to compute Jacobian blocks
compute_Jacobian = @(q_vec, vec) 2 * [ ...
    q_vec(1),  q_vec(4), -q_vec(3);
   -q_vec(4),  q_vec(1),  q_vec(2);
    q_vec(3), -q_vec(2),  q_vec(1)
] * vec;

compute_Jacobian_alt = @(q_vec, vec) 2 * [ ...
    q_vec(2),  q_vec(3),  q_vec(4);
    q_vec(3), -q_vec(2),  q_vec(1);
    q_vec(4), -q_vec(1), -q_vec(2)
] * vec;

compute_Jacobian_alt2 = @(q_vec, vec) 2 * [ ...
   -q_vec(3),  q_vec(2), -q_vec(1);
    q_vec(2),  q_vec(3),  q_vec(4);
    q_vec(1),  q_vec(4), -q_vec(3)
] * vec;

compute_Jacobian_alt3 = @(q_vec, vec) 2 * [ ...
   -q_vec(4),  q_vec(1),  q_vec(2);
   -q_vec(1), -q_vec(4),  q_vec(3);
    q_vec(2),  q_vec(3),  q_vec(4)
] * vec;

% Magnetic field Jacobians
B1 = compute_Jacobian(q, b);
B2 = compute_Jacobian_alt(q, b);
B3 = compute_Jacobian_alt2(q, b);
B4 = compute_Jacobian_alt3(q, b);

% Sun vector Jacobians
S1 = compute_Jacobian(q, s);
S2 = compute_Jacobian_alt(q, s);
S3 = compute_Jacobian_alt2(q, s);
S4 = compute_Jacobian_alt3(q, s);

% Assemble Measurement Jacobian H
H_t = [B1, B2, B3, B4;
       S1, S2, S3, S4];

end
