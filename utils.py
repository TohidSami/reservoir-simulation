import plotly.graph_objects as go
import numpy as np
import os
import webbrowser

def run_animation(history_sw, grid, title_text="Simulation Result"):
    """
    Displays Interactive 3D Animation.
    FIXED: Explicit coordinate generation to match matrix indexing [z][y][x].
    """
    nz, ny, nx = grid.NZ, grid.NY, grid.NX
    
    print(f"Generating Animation ({nx}x{ny}x{nz})... Please wait.")

    # --- 1. ساخت مختصات دقیق (Explicit Coordinates) ---
    # این بخش حیاتی است: ترتیب حلقه‌ها باید دقیقا مثل ترتیب ذخیره سازی در run_model باشد
    # Order: Z -> Y -> X
    x_coords = []
    y_coords = []
    z_coords = []
    
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                x_coords.append(i)
                y_coords.append(j)
                z_coords.append(k)

    frames = []
    for step, sw_data_3d in enumerate(history_sw):
        flat_values = []
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    val = sw_data_3d[k, j, i]
                    flat_values.append(val)
        
        frames.append(
            go.Frame(
                data=[go.Scatter3d(
                    marker=dict(color=flat_values)
                )],
                name=str(step)
            )
        )

    initial_flat = []
    initial_data = history_sw[0]
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                initial_flat.append(initial_data[k, j, i])


    fig = go.Figure(
        data=[go.Scatter3d(
            x=x_coords,
            y=y_coords,
            z=z_coords,
            mode='markers',
            marker=dict(
                size=20,     
                color=initial_flat,
                colorscale='Jet',
                cmin=0.1, cmax=0.8,
                opacity=0.8,   
                symbol='square',
                colorbar=dict(title='Sw')
            ),
            hovertemplate='X: %{x}<br>Y: %{y}<br>Layer: %{z}<br>Sw: %{marker.color:.3f}<extra></extra>'
        )],
        frames=frames
    )

    # تنظیمات ظاهری و محورها
    fig.update_layout(
        title=title_text,
        width=900, height=700,
        scene=dict(
            xaxis=dict(title='X Block', tickmode='linear', dtick=1, showgrid=True),
            yaxis=dict(title='Y Block', tickmode='linear', dtick=1, showgrid=True),
            zaxis=dict(title='Layer (Z)', tickmode='linear', dtick=1, autorange="reversed"),
            aspectmode='data' # تناسب ۱:۱:۱
        ),
        updatemenus=[{
            "type": "buttons",
            "buttons": [
                {"label": "Play", "method": "animate", "args": [None, {"frame": {"duration": 100, "redraw": True}, "fromcurrent": True}]},
                {"label": "Pause", "method": "animate", "args": [[None], {"frame": {"duration": 0, "redraw": False}, "mode": "immediate"}]}
            ],
            "direction": "left", "pad": {"r": 10, "t": 87}, "showactive": False, "x": 0.1, "xanchor": "right", "y": 0, "yanchor": "top"
        }],
        sliders=[{
            "active": 0, "yanchor": "top", "xanchor": "left",
            "currentvalue": {"font": {"size": 20}, "prefix": "Step: ", "visible": True, "xanchor": "right"},
            "pad": {"b": 10, "t": 50}, "len": 0.9, "x": 0.1, "y": 0,
            "steps": [{"args": [[str(k)], {"frame": {"duration": 300, "redraw": True}, "mode": "immediate"}], "label": str(k), "method": "animate"} for k in range(len(history_sw))]
        }]
    )

    # ذخیره و نمایش
    file_path = os.path.abspath("simulation_result_3d.html")
    fig.write_html(file_path)
    print(f"Animation saved to: {file_path}")
    webbrowser.open(f"file://{file_path}")