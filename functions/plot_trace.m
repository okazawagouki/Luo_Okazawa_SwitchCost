function lh = plot_trace(x, y, se, color)
hold on;


if isempty(se)
    se = zeros(size(y));
end

palecol = color * .5 + [1 1 1] * .5;
x2 = [x(:); flipud(x(:))];
y2 = [y(:)-se(:); flipud(y(:)+se(:))];
h = fill(x2,y2,palecol);
set(h, 'EdgeColor', palecol);

    %draw main trace
lh = plot(x,y,'Color', color, 'LineWidth', 1);

