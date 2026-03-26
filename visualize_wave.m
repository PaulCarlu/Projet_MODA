function visualize_wave(wavebasis,mesh, phi, varargin)
% VISUALIZE     Plot a mode shape on 2D triangular mesh
%               Supports real and complex eigenvectors.
%
%   visualize(mesh, phi)
%   visualize(mesh, phi, 'animate', true)
%   visualize(mesh, phi, 'scale', 0.2, 'animate', true, 'nframes', 60)
%
% Inputs:
%   mesh      struct with .nodes (N×2) and .conn (M×3)
%   phi       (2N)×1 vector — can be real or complex
%   Name-Value pairs (optional):
%     'scale'     exaggeration factor (default: auto)
%     'animate'   logical — create animation of Re(φ e^{iθ})  (default: false)
%     'nframes'   number of frames in animation (default: 60)
%     'fps'       animation speed (default: 30)

    nodes = mesh.nodes;
    conn  = mesh.conn;
    N     = size(nodes,1);

    % ─────────────── Parse optional arguments ───────────────
    p = inputParser;
    addParameter(p, 'scale',   []);           % empty → auto
    addParameter(p, 'animate', false, @islogical);
    addParameter(p, 'nframes', 60,   @isscalar);
    addParameter(p, 'fps',     40,   @isscalar);
    parse(p, varargin{:});

    animate = p.Results.animate;
    nframes = p.Results.nframes;
    fps     = p.Results.fps;

    % ─────────────── Extract real & imaginary parts ───────────────
    if isreal(phi)
        ux = real(phi(1:2:end));
        uy = real(phi(2:2:end));
        title_suffix = ' (real part)';
    else
        ux = real(phi(1:2:end));
        uy = real(phi(2:2:end));
        title_suffix = ' – Re(φ)';
    end

    % ─────────────── Scaling ───────────────
    if isempty(p.Results.scale)
        amp = max(sqrt(ux.^2 + uy.^2));
        if amp > 1e-12
            scale = 0.1 * max(max(abs(nodes))) / amp;
        else
            scale = 0.1;
        end
    else
        scale = p.Results.scale;
    end

    % ─────────────── Static plot (always shown) ───────────────
    figure('Color','w', 'Name','Mode shape');

    % Original (reference) mesh
    triplot(conn, nodes(:,1), nodes(:,2), ...
        'Color',[0.78 0.78 0.78], 'LineWidth',0.6);

    hold on;

    % Deformed – real part
    x_def = nodes(:,1) + scale * ux;
    y_def = nodes(:,2) + scale * uy;

    disp_mag = sqrt(ux.^2 + uy.^2);

    patch('Faces',       conn, ...
          'Vertices',    [x_def y_def], ...
          'FaceVertexCData', disp_mag, ...
          'FaceColor',   'interp', ...
          'EdgeColor',   'interp', ...
          'LineWidth',   1.1);

    colormap jet;
    cb = colorbar('eastoutside');
    cb.Label.String = 'Displacement magnitude (real part)';
    cb.FontSize = 10;

    axis equal tight; box on; grid on;
    xlabel('x'); ylabel('y');
    title_str = sprintf('Mode shape – scale = %.2f%s', scale, title_suffix);
    title(title_str, 'Interpreter','none');
    set(gca, 'FontSize',11);

    text(0.02, 0.98, sprintf('Scale = %.2f', scale), ...
         'Units','normalized', 'VerticalAlignment','top', ...
         'FontSize',9, 'Color',[0.4 0.4 0.4]);

    hold off;

    % ─────────────── Optional animation of complex mode ───────────────
    if ~animate || isreal(phi)
        return;
    end

    fprintf('Creating animation of Re(φ ⋅ e^{iθ}) ...\n');

    % Prepare figure for animation
    fig_anim = figure('Color','w', 'Name','Complex mode animation');
    h_patch = patch('Faces',conn, 'Vertices',zeros(N,2), ...
                    'FaceVertexCData',zeros(size(conn,1),1), ...
                    'FaceColor','interp', 'EdgeColor','interp');

    colormap jet;
    cb = colorbar('eastoutside');
    cb.Label.String = 'Displacement magnitude';
    axis equal tight; box on; grid on;
    xlabel('x'); ylabel('y');
    title('Re(φ e^{iθ})  –  animation');
    set(gca,'FontSize',11);

    % Animation loop
    theta_vec = linspace(0, 2*pi, nframes+1);
    theta_vec(end) = [];   % avoid duplicating frame 0

    for k = 1:nframes
        theta = theta_vec(k);

        % Complex rotation: Re( φ * exp(iθ) )
        phi_rot = phi * exp(1i * theta);
        ux_rot = real(phi_rot(1:2:end));
        uy_rot = real(phi_rot(2:2:end));

        x_def_rot = nodes(:,1) + scale * ux_rot;
        y_def_rot = nodes(:,2) + scale * uy_rot;

        mag_rot = sqrt(ux_rot.^2 + uy_rot.^2);

        set(h_patch, ...
            'Vertices',    [x_def_rot y_def_rot], ...
            'FaceVertexCData', mag_rot);

        title(sprintf('Re(φ e^{iθ})   θ = %.1f°   scale = %.2f', ...
                      rad2deg(theta), scale), 'Interpreter','none');

        drawnow;

        if k == 1
            set(gca, 'XLimMode','manual', 'YLimMode','manual');
        end

        pause(1/fps);
    end

    fprintf('Animation complete.\n');
end

end

