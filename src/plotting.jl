export zz_heatmap_figure


function zz_heatmap_figure(
    Z;
    x,
    y,
    xlabel = L"J_{12}\ (\mathrm{MHz})",
    ylabel = L"J_c\ (\mathrm{MHz})",
    title = "",
    colorbar_label = L"\zeta\ (\mathrm{kHz})",
    colormap = Reverse(:Blues),
    clims = (1e-2, 1e2),          # kHz range
    dashed_curve = nothing,      # (x_curve, y_curve)
    device_point = nothing,      # (x0, y0)
    device_label = "Our device",
)

    fig = Figure(resolution = (700, 550))

    ax = Axis(
        fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        limits = (extrema(x)..., extrema(y)...),
    )

    # --- Log-scaled heatmap ---
    hm = heatmap!(
        ax,
        x,
        y,
        Z;
        colormap = colormap,
        colorscale = log10,
    )

    # --- Colorbar ---
    Colorbar(
        fig[1, 2],
        hm;
        label = colorbar_label,
        ticks = ([1e-2, 1e-1, 1e0, 1e1, 1e2],
                 ["10⁻²", "10⁻¹", "10⁰", "10¹", "10²"]),
    )

    # --- Dashed minimum-ZZ curve ---
    if dashed_curve !== nothing
        xcurve, ycurve = dashed_curve
        lines!(
            ax,
            xcurve,
            ycurve;
            linestyle = :dash,
            linewidth = 2,
            color = :black,
        )
        text!(
            ax,
            "Minimum\nZZ";
            position = (xcurve[end], ycurve[end]),
            align = (:right, :bottom),
        )
    end

    # --- Device marker ---
    if device_point !== nothing
        x0, y0 = device_point
        scatter!(
            ax,
            [x0],
            [y0];
            markersize = 14,
            color = :dodgerblue,
        )
        text!(
            ax,
            device_label;
            position = (x0, y0),
            offset = (10, -10),
            color = :dodgerblue,
        )
    end

    # --- Extra annotation ---
    text!(
        ax,
        "No coupler";
        position = (maximum(x) * 0.6, minimum(y) + 40),
        color = :gray,
    )

    return fig, ax
end
