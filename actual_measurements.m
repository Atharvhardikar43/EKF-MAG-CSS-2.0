function [sun_meas_body, sun_meas_eci, mag_meas_body, mag_meas_eci] = ...
    actual_measurements(sun_true_body, sun_true_eci, mag_true_body, mag_true_eci, noise_css, noise_mag)

sun_meas_body = sun_true_body + noise_css * randn(size(sun_true_body));
sun_meas_body = sun_meas_body / norm(sun_meas_body);

mag_meas_body = mag_true_body + noise_mag * randn(size(mag_true_body));
mag_meas_body = mag_meas_body / norm(mag_meas_body);

sun_meas_eci = sun_true_eci + noise_css * randn(size(sun_true_eci));
sun_meas_eci = sun_meas_eci / norm(sun_meas_eci);

mag_meas_eci = mag_true_eci + noise_mag * randn(size(mag_true_eci));
mag_meas_eci = mag_meas_eci / norm(mag_meas_eci);

end
