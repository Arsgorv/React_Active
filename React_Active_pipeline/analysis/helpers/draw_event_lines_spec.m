function draw_event_lines_spec(ax, ylims_use, t_stimoff, t_arr, t_rew)
axes(ax);
hold on
plot([0 0], ylims_use, 'k-', 'LineWidth',1.2);
if isfinite(t_stimoff) && t_stimoff > 0
    plot([t_stimoff t_stimoff], ylims_use, 'k-', 'LineWidth',1.2);
end
if isfinite(t_arr)
    plot([t_arr t_arr], ylims_use, 'k--', 'LineWidth',1.2);
end
if isfinite(t_rew)
    plot([t_rew t_rew], ylims_use, 'r--', 'LineWidth',1.2);
end
hold off
end
