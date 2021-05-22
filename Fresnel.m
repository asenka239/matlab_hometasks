% To have 1x3992 vector
temp_var = 0:0.20:399;
temp_var = sort([-sqrt(temp_var), sqrt(temp_var)]);

% Create x, y
x = fresnelc(temp_var);
y = fresnels(temp_var);

% Plot figure
figure; hold on; grid on; plot(x, y);
