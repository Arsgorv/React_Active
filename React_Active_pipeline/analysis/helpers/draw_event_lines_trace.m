function draw_event_lines_trace(ax, yl, t_stimoff, t_arr, t_rew)
axes(ax);
hold on
if isfinite(t_stimoff) && t_stimoff > 0
    patch([0 t_stimoff t_stimoff 0], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.7 0.7 0.7], 'EdgeColor','none', 'FaceAlpha',0.18);
end
plot([0 0], yl, 'k-', 'LineWidth',1.2);
if isfinite(t_arr)
    plot([t_arr t_arr], yl, 'k--', 'LineWidth',1.2);
end
if isfinite(t_rew)
    plot([t_rew t_rew], yl, 'r--', 'LineWidth',1.2);
end
hold off
end
