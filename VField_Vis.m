%input here R2 or S2 for your desired manifold to form a vector space over
space_type = 'R2';

switch space_type
    case 'R2'
        %input vector field here (first argument the d/dx term and second
        %the d/dy term)
        %if you dont know what to do try [-y,x] i guess, or anything its
        %your choice
        vector_field = @(x,y) [y^2, x^2];
        
        
    case 'S2'
     %input vector field here (same as above)
     %if you dont ... try 0.1sin(phi),0.2*cos(theta) or even 
        vector_field = @(theta,phi) [0.1*sin(phi), 0.2*cos(theta)];
        
    otherwise
        error('Invalid space type. Choose R2, or S2');
end
%% 

visualizeVectorField(space_type, vector_field)

function visualizeVectorField(space_type, vector_field)

    figure('Position', [100, 100, 800, 600]);
    
    switch space_type
        case 'R2'
            visualizeR2(vector_field);
        case 'S2'
            visualizeS2(vector_field);
    end

    title(sprintf('Vector Field on %s', space_type), 'FontSize', 14);
    grid on;
    axis equal;

    caption = sprintf('Generated: %s | Space: %s', datetime("now"), space_type);
    annotation('textbox', [0.1, 0.02, 0.8, 0.05], ...
        'String', caption, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end
%% 

function visualizeR2(vec_field)
  
    [X, Y] = meshgrid(-3:0.5:3, -3:0.5:3);
    
    U = zeros(size(X));
    V = zeros(size(Y));
    
    for i = 1:numel(X)
        vec = vec_field(X(i), Y(i));
        U(i) = vec(1);
        V(i) = vec(2);
    end
    
    %store magnitudes to allow colouring for later
    magnitude = sqrt(U.^2 + V.^2);
    
    % norm
    U_norm = U ./ (magnitude + eps);
    V_norm = V ./ (magnitude + eps);
    
    
    h = quiver(X, Y, U_norm, V_norm, 0.5, 'LineWidth', 1.0);
    hold on;
    
    
  
    %prep colour map
    cmap = parula(256);
    mag_flat = magnitude(:);
    mag_min = min(mag_flat);
    mag_max = max(mag_flat);
    if mag_max == mag_min
        mag_max = mag_min + 1; % avoid division by zero
    end
    mag_idx = round( ( (mag_flat - mag_min) / (mag_max - mag_min) ) * 255 ) + 1;
    mag_idx = reshape(mag_idx, size(magnitude));
    
    %because quiver wont allow colouring vectors by a rule afaia so redraw
    %the vectors using colours
    delete(h);
    
    scale = 0.4; 
    head_frac = 0.2;
    for i = 1:numel(X)
        xi = X(i);
        yi = Y(i);
        ui = U_norm(i) * scale;
        vi = V_norm(i) * scale;
        
        x_shaft = [xi, xi + ui*(1-head_frac)];
        y_shaft = [yi, yi + vi*(1-head_frac)];
        color = cmap(mag_idx(i), :);
        plot(x_shaft, y_shaft, 'Color', color, 'LineWidth', 1.5);
    
        dir = [ui, vi];
        if norm(dir) < eps
            continue;
        end
        dir_perp = [-dir(2), dir(1)];
        dir_perp = dir_perp / norm(dir_perp);
        tip = [xi + ui, yi + vi];
        base = [xi + ui*(1-head_frac), yi + vi*(1-head_frac)];
        head_width = 0.12; % constant width
        p1 = tip;
        p2 = base + dir_perp * head_width;
        p3 = base - dir_perp * head_width;
        patch('XData', [p1(1), p2(1), p3(1)], 'YData', [p1(2), p2(2), p3(2)], ...
            'FaceColor', color, 'EdgeColor', color);
    end
    
    colormap(cmap);
    c = colorbar;
    c.Label.String = 'Magnitude';
    clim([mag_min, mag_max]);
    
    xlabel('x', 'FontSize', 12);
    ylabel('y', 'FontSize', 12);
    xlim([-3.5, 3.5]);
    ylim([-3.5, 3.5]);
    hold off;
end

%%

 function visualizeS2(vector_field)
    
    n_theta = 28;   %longitude 
    n_phi = 20;    %latitude
    
    phi = linspace(0.15, pi-0.15, n_phi);
    theta = linspace(0, 2*pi, n_theta);
    [Theta, Phi] = meshgrid(theta, phi);
    
    
    X = sin(Phi) .* cos(Theta);
    Y = sin(Phi) .* sin(Theta);
    Z = cos(Phi);
    
    
    dtheta = zeros(size(Theta));
    dphi = zeros(size(Phi));
    
    for i = 1:size(Theta, 1)
        for j = 1:size(Theta, 2)
            vec = vector_field(Theta(i,j), Phi(i,j));
            dtheta(i,j) = vec(1);
            dphi(i,j) = vec(2);
        end
    end
    
    %create cartesian tangent vectors for surface of sphere
    e_theta_x = -sin(Theta);
    e_theta_y = cos(Theta);
    e_theta_z = zeros(size(Theta));
    
    e_phi_x = cos(Theta) .* cos(Phi);
    e_phi_y = sin(Theta) .* cos(Phi);
    e_phi_z = -sin(Phi);
    
    
    Vx = dtheta .* e_theta_x + dphi .* e_phi_x;
    Vy = dtheta .* e_theta_y + dphi .* e_phi_y;
    Vz = dtheta .* e_theta_z + dphi .* e_phi_z;
    
    
    magnitude = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
    
    %colour mapping by magnitude (pop pop) ...
    cmap = parula(256);
    mag_flat = magnitude(:);
    mag_min = min(mag_flat);
    mag_max = max(mag_flat);
    
    if mag_max == mag_min
        mag_max = mag_min + 1; % so no 0 div again
    end
    
    mag_idx = round(((mag_flat - mag_min) / (mag_max - mag_min)) * 255) + 1;
    mag_idx = reshape(mag_idx, size(magnitude));
    
    %normalise
    eps_val = 1e-2;
    Vx_norm = Vx ./ (magnitude + eps_val);
    Vy_norm = Vy ./ (magnitude + eps_val);
    Vz_norm = Vz ./ (magnitude + eps_val);
    
 
    scale = 0.22;
    Vx_plot = Vx_norm * scale;
    Vy_plot = Vy_norm * scale;
    Vz_plot = Vz_norm * scale;
    
 
  
    
    surf(X, Y, Z, ...
        'FaceColor', '#9266a9', ...
        'FaceAlpha', 1, ...
        'EdgeColor', 'none', ...
        'SpecularStrength', 0.1, ...
        'AmbientStrength', 0.7);
    hold on;
    
   %plot small individual quivers for each amgnitude grouping (better for
   %large amount of vectors with similar magnitude which is probably likely
   %based on regular fields although depends)
    for i = 1:size(X, 1)
        for j = 1:size(X, 2)
            if magnitude(i,j) > eps_val
                quiver3(X(i,j), Y(i,j), Z(i,j), ...
                    Vx_plot(i,j), Vy_plot(i,j), Vz_plot(i,j), ...
                    0, ...
                    'LineWidth', 1.4, ...
                    'Color', cmap(mag_idx(i,j), :), ...
                    'MaxHeadSize', 2, ...
                    'AutoScale', 'off');
            end
        end
    end
    
    
    n_grid = 40;
    [theta_grid, phi_grid] = meshgrid(linspace(0, 2*pi, n_grid), ...
                                       linspace(0.1, pi-0.1, n_grid/2));
    x_grid = sin(phi_grid) .* cos(theta_grid);
    y_grid = sin(phi_grid) .* sin(theta_grid);
    z_grid = cos(phi_grid);
    
    plot3(x_grid, y_grid, z_grid, 'Color', [0.65, 0.65, 0.65], 'LineWidth', 0.2);
    plot3(x_grid', y_grid', z_grid', 'Color', [0.65, 0.65, 0.65], 'LineWidth', 0.2);
    
    
    axis equal tight;
    axis off;
    box off;
    axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]);
    
    
    lighting gouraud;
    light('Position', [1, 1, 1], 'Style', 'infinite');
    light('Position', [-1, -1, 0.5], 'Style', 'infinite');
    
  
    view(60, 25);
    
  
    rotate3d on;
    
  
    colormap(gca, cmap);
    c = colorbar('Location', 'eastoutside', 'FontSize', 11);
    c.Label.String = 'Vector Magnitude';
    c.Label.FontSize = 12;
    c.Label.FontWeight = 'bold';
    
    %show min/max
    c.Ticks = linspace(0, 1, 5);
    c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), ...
        linspace(mag_min, mag_max, 5), 'UniformOutput', false);
    
    
    title('S^2 Vector Field (Parula = Magnitude)', ...
        'FontSize', 14, 'FontWeight', 'normal', 'Color', [0.2, 0.2, 0.2]);
    
    hold off;
end


