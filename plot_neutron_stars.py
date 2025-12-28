import gzip
import json
import argparse
from math import isnan, log10
import time
import numpy as np
import signal
import os

# Global flag for graceful shutdown
shutdown_requested = False


def signal_handler(sig, frame):
    global shutdown_requested
    shutdown_requested = True
    print("\n[SIGINT] Graceful shutdown requested...", flush=True)


def get_first(*keys, src=None):
    if not src:
        return None
    for k in keys:
        if k in src and src[k] is not None:
            return src[k]
    return None


def try_float(v):
    try:
        return float(v)
    except Exception:
        return None


def _clean_line_for_json(line):
    line = line.lstrip()
    if not line:
        return line
    if line == "[" or line == "]":
        return ""
    if line.startswith("["):
        line = line[1:].lstrip()
    while line.endswith(",") or line.endswith("]"):
        line = line[:-1].rstrip()
    return line


def extract_neutron_stars(gz_path, max_points=None, progress_interval=200000):
    global shutdown_requested
    points = []
    lines_read = 0
    systems_parsed = 0
    malformed_lines = 0
    stars_found = 0
    start = time.time()

    with gzip.open(gz_path, "rt", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            if shutdown_requested:
                elapsed = time.time() - start
                rate = lines_read / elapsed if elapsed > 0 else 0
                print(
                    f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] Shutdown requested — stopping read.",
                    flush=True,
                )
                print(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] lines={lines_read} systems={systems_parsed} neutron_stars={stars_found} malformed={malformed_lines} elapsed={elapsed:.0f}s rate={rate:.0f} L/s",
                    flush=True,
                )
                break

            lines_read += 1
            line = raw.strip()
            if not line:
                continue

            line = _clean_line_for_json(line)
            if not line:
                continue

            try:
                sysobj = json.loads(line)
            except Exception:
                malformed_lines += 1
            else:
                systems_parsed += 1
                sys_name = sysobj.get("name") or sysobj.get("systemName") or ""
                bodies = sysobj.get("bodies") or []
                for b in bodies:
                    if (
                        b.get("subType") != "Neutron Star"
                        and b.get("subtype") != "Neutron Star"
                    ):
                        continue
                    sr = try_float(
                        get_first("solarRadius", "radius", "stellarRadius", src=b)
                    )
                    rp = try_float(
                        get_first(
                            "rotationalPeriod", "rotationPeriod", "rotPeriod", src=b
                        )
                    )
                    age = try_float(
                        get_first("age", "age_Myr", "age_million_years", src=b)
                    )
                    mass = try_float(
                        get_first("solarMasses", "mass", "stellarMass", src=b)
                    )
                    temperature = try_float(
                        get_first("surfaceTemperature", "temperature", "temp", src=b)
                    )
                    if sr is None or rp is None or age is None:
                        continue
                    if isnan(sr) or isnan(rp) or isnan(age):
                        continue
                    if mass is not None and isnan(mass):
                        mass = None
                    if temperature is not None and isnan(temperature):
                        temperature = None
                    points.append(
                        {
                            "system": sys_name,
                            "body": b.get("name") or b.get("bodyName") or "",
                            "solarRadius": sr,
                            "rotationalPeriod": rp,
                            "age": age,
                            "mass": mass,
                            "temperature": temperature,
                        }
                    )
                    stars_found += 1
                    if max_points and len(points) >= max_points:
                        if progress_interval:
                            elapsed = time.time() - start
                            rate = lines_read / elapsed if elapsed > 0 else 0
                            print(
                                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] lines={lines_read} systems={systems_parsed} neutron_stars={stars_found} malformed={malformed_lines} elapsed={elapsed:.0f}s rate={rate:.0f} L/s",
                                flush=True,
                            )
                        return points

            if progress_interval and lines_read % progress_interval == 0:
                elapsed = time.time() - start
                rate = lines_read / elapsed if elapsed > 0 else 0
                print(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] lines={lines_read} systems={systems_parsed} neutron_stars={stars_found} malformed={malformed_lines} elapsed={elapsed:.0f}s rate={rate:.0f} L/s",
                    flush=True,
                )

    if progress_interval:
        elapsed = time.time() - start
        rate = lines_read / elapsed if elapsed > 0 else 0
        print(
            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] DONE lines={lines_read} systems={systems_parsed} neutron_stars={stars_found} malformed={malformed_lines} elapsed={elapsed:.0f}s avg_rate={rate:.0f} L/s",
            flush=True,
        )

    return points


def _range_for(vals, percentile_low=1, percentile_high=99):
    if len(vals) == 0:
        return (0.0, 1.0)
    arr = np.array(vals)
    mn = np.percentile(arr, percentile_low)
    mx = np.percentile(arr, percentile_high)
    if mn == mx:
        mn = arr.min()
        mx = arr.max()
        if mn == mx:
            mx = mn + (abs(mn) if mn != 0 else 1.0) * 1e-6
    return (mn, mx)


def value_to_color(val, vmin, vmax, colorscale="viridis"):
    """Convert value to RGB color using colorscale"""
    # Normalize value to 0-1
    if vmax == vmin:
        t = 0.5
    else:
        t = (val - vmin) / (vmax - vmin)
    t = max(0, min(1, t))

    # Colorscales
    if colorscale == "viridis":
        # Simple viridis approximation
        r = int(255 * (0.267 + t * (0.329 - 0.267 + t * (0.993 - 0.329))))
        g = int(255 * (0.005 + t * (0.569 + t * (0.906 - 0.569))))
        b = int(255 * (0.329 + t * (0.527 + t * (0.144 - 0.527))))
    elif colorscale == "plasma":
        # Simple plasma approximation
        r = int(255 * (0.050 + t * (0.940 + t * (0.988 - 0.940))))
        g = int(255 * (0.030 + t * (0.384 + t * (0.645 - 0.384))))
        b = int(255 * (0.528 + t * (0.885 - t * (0.775))))
    else:
        # Default grayscale
        v = int(255 * t)
        r, g, b = v, v, v

    return f"rgb({r},{g},{b})"


def create_svg_plot(
    x_vals,
    y_vals,
    colors,
    color_range,
    colorscale,
    systems,
    bodies,
    ages,
    radii,
    rps,
    xlabel,
    ylabel,
    title,
    colorbar_title,
    width=1000,
    height=700,
    margin=80,
    reference_line=None,
    reference_line_label=None,
):
    """Create Plotly WebGL scatter plot with log scales"""
    import random

    # Filter out invalid values for log scale
    valid_indices = []
    for i in range(len(x_vals)):
        if x_vals[i] > 0 and y_vals[i] > 0:
            valid_indices.append(i)

    if not valid_indices:
        return "<div><p>No valid data points</p></div>"

    x_vals = [x_vals[i] for i in valid_indices]
    y_vals = [y_vals[i] for i in valid_indices]
    colors = [colors[i] for i in valid_indices]
    systems = [systems[i] for i in valid_indices]
    bodies = [bodies[i] for i in valid_indices]
    ages = [ages[i] for i in valid_indices]
    radii = [radii[i] for i in valid_indices]
    rps = [rps[i] for i in valid_indices]

    # Generate unique ID
    uid = str(random.randint(100000, 999999))

    # Create hover text
    hover_text = []
    for i in range(len(x_vals)):
        hover_text.append(
            f"{bodies[i]}<br>"
            + f"System: {systems[i]}<br>"
            + f"Age: {ages[i]:.2e}<br>"
            + f"Radius: {radii[i]:.4f}<br>"
            + f"Rot Period: {rps[i]:.4f}"
        )

    # Create the HTML with Plotly
    html = f"""
    <div id="plot-{uid}" style="width:100%;height:{height}px;"></div>
    <script type="application/json" id="plotData-{uid}">
    {{
        "x": {json.dumps(x_vals)},
        "y": {json.dumps(y_vals)},
        "colors": {json.dumps(colors)},
        "hover_text": {json.dumps(hover_text)},
        "systems": {json.dumps(systems)},
        "color_range": {json.dumps(color_range)},
        "colorscale": "{colorscale}",
        "xlabel": "{xlabel}",
        "ylabel": "{ylabel}",
        "title": "{title}",
        "colorbar_title": "{colorbar_title}",
        "reference_line": {json.dumps(reference_line)},
        "reference_line_label": "{reference_line_label if reference_line_label else ''}"
    }}
    </script>
    <script>
    (function() {{
        const data = JSON.parse(document.getElementById('plotData-{uid}').textContent);
        
        const trace = {{
            x: data.x,
            y: data.y,
            mode: 'markers',
            type: 'scattergl',  // WebGL for performance
            marker: {{
                size: 4,
                color: data.colors,
                colorscale: data.colorscale === 'viridis' ? 'Viridis' : 'Plasma',
                showscale: true,
                colorbar: {{
                    title: data.colorbar_title,
                    thickness: 20,
                    len: 0.7
                }},
                cmin: data.color_range[0],
                cmax: data.color_range[1]
            }},
            text: data.hover_text,
            hoverinfo: 'text',
            customdata: data.systems
        }};
        
        const layout = {{
            title: data.title,
            xaxis: {{
                title: data.xlabel,
                type: 'log',
                gridcolor: '#444',
                zerolinecolor: '#666'
            }},
            yaxis: {{
                title: data.ylabel,
                type: 'log',
                gridcolor: '#444',
                zerolinecolor: '#666',
                range: data.reference_line ? [Math.log10(Math.min(...data.y) * 0.1), Math.log10(Math.max(...data.y, data.reference_line) * 10)] : undefined
            }},
            paper_bgcolor: '#1a1a1a',
            plot_bgcolor: '#1a1a1a',
            font: {{ color: '#fff' }},
            hovermode: 'closest',
            autosize: true,
            height: {height},
            margin: {{ l: 80, r: 80, t: 80, b: 80 }},
            shapes: data.reference_line ? [{{
                type: 'line',
                x0: Math.min(...data.x),
                x1: Math.max(...data.x),
                y0: data.reference_line,
                y1: data.reference_line,
                line: {{
                    color: '#00ff00',
                    width: 2,
                    dash: 'dash'
                }},
                label: {{
                    text: data.reference_line_label,
                    textposition: 'top right',
                    textangle: 0
                }}
            }}] : []
        }};
        
        const config = {{
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: ['select2d', 'lasso2d'],
            displaylogo: false
        }};
        
        Plotly.newPlot('plot-{uid}', [trace], layout, config);
        
        // Add click handler to open system in Canonn
        document.getElementById('plot-{uid}').on('plotly_click', function(data) {{
            if (data.points && data.points[0]) {{
                const system = data.points[0].customdata;
                const url = 'https://signals.canonn.tech/?system=' + encodeURIComponent(system);
                window.open(url, '_blank');
            }}
        }});
    }})();
    </script>
    """

    return html


def create_html_with_tabs(svg_files, output_file, point_count):
    """Create HTML page with tabbed interface and embedded SVG plots"""

    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Neutron Stars Analysis - {point_count} stars</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            background: #1a1a1a;
            color: #e0e0e0;
        }}
        .header {{
            background: #2d2d2d;
            padding: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.3);
        }}
        h1 {{
            font-size: 24px;
            font-weight: 600;
        }}
        .subtitle {{
            color: #888;
            font-size: 14px;
            margin-top: 5px;
        }}
        .tabs {{
            display: flex;
            background: #2d2d2d;
            padding: 0 20px;
            border-bottom: 1px solid #444;
            overflow-x: auto;
        }}
        .tab {{
            padding: 15px 25px;
            cursor: pointer;
            border: none;
            background: none;
            color: #888;
            font-size: 14px;
            font-weight: 500;
            white-space: nowrap;
            transition: all 0.2s;
            border-bottom: 3px solid transparent;
        }}
        .tab:hover {{
            color: #fff;
            background: rgba(255,255,255,0.05);
        }}
        .tab.active {{
            color: #4a9eff;
            border-bottom-color: #4a9eff;
        }}
        .tab-content {{
            display: none;
            padding: 20px;
            animation: fadeIn 0.3s;
        }}
        .tab-content.active {{
            display: block;
        }}
        @keyframes fadeIn {{
            from {{ opacity: 0; }}
            to {{ opacity: 1; }}
        }}
        .plot-container {{
            background: #2d2d2d;
            border-radius: 8px;
            padding: 20px;
            margin-bottom: 20px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.2);
            overflow-x: auto;
        }}
        .loading {{
            text-align: center;
            padding: 60px;
            color: #888;
        }}
        .spinner {{
            border: 3px solid #444;
            border-top: 3px solid #4a9eff;
            border-radius: 50%;
            width: 40px;
            height: 40px;
            animation: spin 1s linear infinite;
            margin: 0 auto 15px;
        }}
        @keyframes spin {{
            0% {{ transform: rotate(0deg); }}
            100% {{ transform: rotate(360deg); }}
        }}
        svg {{
            max-width: 100%;
            height: auto;
        }}
        .help {{
            background: #2d2d2d;
            padding: 15px;
            margin-bottom: 20px;
            border-radius: 8px;
            font-size: 13px;
            color: #aaa;
        }}
        .help strong {{
            color: #4a9eff;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Neutron Stars Analysis</h1>
        <div class="subtitle">{point_count} neutron stars analyzed | Interactive controls: zoom, pan, box select</div>
    </div>
    
    <div class="tabs">
        <button class="tab active" onclick="showTab(0)">Radius vs Rotational Period</button>
        <button class="tab" onclick="showTab(1)">Rotational Period vs Age</button>
        <button class="tab" onclick="showTab(2)">Radius vs Age</button>
        <button class="tab" onclick="showTab(3)">Tangential Velocity vs Age</button>
        <button class="tab" onclick="showTab(4)">Temperature vs Age</button>
        <button class="tab" onclick="showTab(5)">Mass vs Radius</button>
    </div>
    
    <div id="tab0" class="tab-content active">
        <div class="help">
            <strong>Controls:</strong> Zoom In/Out buttons, Reset view, Pan mode (drag to move), Box Select (drag to zoom to region), Mouse wheel to zoom at cursor
        </div>
        <div class="plot-container">
            <div class="loading">
                <div class="spinner"></div>
                Loading plot...
            </div>
        </div>
    </div>
    
    <div id="tab1" class="tab-content">
        <div class="help">
            <strong>Controls:</strong> Click and drag to pan, Mouse wheel to zoom, Click points to open in Canonn, Hover for details | WebGL rendering for millions of points
        </div>
        <div class="plot-container">
            <div class="loading">
                <div class="spinner"></div>
                Loading plot...
            </div>
        </div>
    </div>
    
    <div id="tab2" class="tab-content">
        <div class="help">
            <strong>Controls:</strong> Click and drag to pan, Mouse wheel to zoom, Click points to open in Canonn, Hover for details | WebGL rendering for millions of points
        </div>
        <div class="plot-container">
            <div class="loading">
                <div class="spinner"></div>
                Loading plot...
            </div>
        </div>
    </div>
    
    <div id="tab3" class="tab-content">
        <div class="help">
            <strong>Controls:</strong> Click and drag to pan, Mouse wheel to zoom, Click points to open in Canonn, Hover for details | Tangential velocity = 2πr/T
        </div>
        <div class="plot-container">
            <div class="loading">
                <div class="spinner"></div>
                Loading plot...
            </div>
        </div>
    </div>

    <div id="tab4" class="tab-content">
        <div class="help">
            <strong>Controls:</strong> Click and drag to pan, Mouse wheel to zoom, Click points to open in Canonn, Hover for details | Shows how neutron star temperature correlates with age
        </div>
        <div class="plot-container">
            <div class="loading">
                <div class="spinner"></div>
                Loading plot...
            </div>
        </div>
    </div>

    <div id="tab5" class="tab-content">
        <div class="help">
            <strong>Controls:</strong> Click and drag to pan, Mouse wheel to zoom, Click points to open in Canonn, Hover for details | Displays mass-radius relationship with temperature coloring
        </div>
        <div class="plot-container">
            <div class="loading">
                <div class="spinner"></div>
                Loading plot...
            </div>
        </div>
    </div>

    <script>
        const svgFiles = {svg_files_json};
        const loadedTabs = new Set();
        
        function showTab(index) {{
            const tabs = document.querySelectorAll('.tab');
            const contents = document.querySelectorAll('.tab-content');
            
            tabs.forEach((tab, i) => {{
                tab.classList.toggle('active', i === index);
            }});
            
            contents.forEach((content, i) => {{
                content.classList.toggle('active', i === index);
            }});
            
            if (!loadedTabs.has(index)) {{
                loadSVG(index);
            }}
        }}
        
        function loadSVG(index) {{
            const tabContent = document.getElementById('tab' + index);
            const container = tabContent.querySelector('.plot-container');
            const svgPath = svgFiles[index];
            
            fetch(svgPath)
                .then(response => response.text())
                .then(content => {{
                    container.innerHTML = content;
                    loadedTabs.add(index);
                    
                    // Re-execute scripts after a brief delay to ensure DOM is ready
                    setTimeout(() => {{
                        const scripts = container.querySelectorAll('script');
                        scripts.forEach(oldScript => {{
                            if (oldScript.type === 'application/json') {{
                                // Skip JSON data scripts
                                return;
                            }}
                            const newScript = document.createElement('script');
                            if (oldScript.src) {{
                                newScript.src = oldScript.src;
                            }} else {{
                                newScript.textContent = oldScript.textContent;
                            }}
                            oldScript.parentNode.replaceChild(newScript, oldScript);
                        }});
                    }}, 100);
                }})
                .catch(error => {{
                    container.innerHTML = '<div style="color: #f44336; padding: 40px; text-align: center;">Error loading plot: ' + error.message + '</div>';
                }});
        }}
        
        window.addEventListener('DOMContentLoaded', () => {{
            loadSVG(0);
        }});
    </script>
</body>
</html>"""

    svg_files_json = json.dumps(svg_files)
    html = html_template.format(point_count=point_count, svg_files_json=svg_files_json)

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html)


def main():
    signal.signal(signal.SIGINT, signal_handler)

    p = argparse.ArgumentParser(description="Neutron star plots")
    p.add_argument(
        "--input",
        "-i",
        default="/home/meddler/spansh/galaxy.json.gz",
        help="path to galaxy.json.gz",
    )
    p.add_argument(
        "--output", "-o", default="neutron_stars_plots.html", help="output HTML file"
    )
    p.add_argument("--max", type=int, default=0, help="max points to read (0 = all)")
    p.add_argument(
        "--progress-interval",
        type=int,
        default=200000,
        help="lines between progress messages (0 to disable)",
    )
    args = p.parse_args()

    print("Extracting neutron stars...")
    pts = extract_neutron_stars(
        args.input,
        max_points=(args.max if args.max > 0 else None),
        progress_interval=args.progress_interval,
    )
    if not pts:
        print("No neutron star data found.")
        return

    print(f"Building plots for {len(pts)} neutron stars...")

    ages = np.array([p["age"] for p in pts])
    radii = np.array([p["solarRadius"] for p in pts])
    rps = np.array([p["rotationalPeriod"] for p in pts])
    masses = np.array(
        [p["mass"] if p["mass"] is not None else 1.4 for p in pts]
    )  # Default to 1.4 solar masses
    temperatures = np.array(
        [p["temperature"] if p["temperature"] is not None else 6000 for p in pts]
    )  # Default to 6000 K
    systems = [p["system"] for p in pts]
    bodies = [p["body"] for p in pts]

    # Calculate tangential velocity: v = 2πr/T at equator
    # radius in solar radii, convert to km (1 solar radius = 696000 km)
    # period in DAYS, convert to seconds
    # Result in km/s
    radius_km = radii * 696000

    # Period is in days, convert to seconds
    period_seconds = rps * 86400

    tangential_velocity = (2 * np.pi * radius_km) / period_seconds

    age_range = _range_for(ages)
    radius_range = _range_for(radii)
    rp_range = _range_for(rps)
    mass_range = _range_for(masses)
    tang_vel_range = _range_for(tangential_velocity)
    temp_range = _range_for(temperatures)

    print(
        f"Color ranges: Age={age_range}, Radius={radius_range}, RotPeriod={rp_range}, Mass={mass_range}, TangVel={tang_vel_range}, Temp={temp_range}"
    )

    output_dir = os.path.dirname(args.output) or "."
    base_name = os.path.splitext(os.path.basename(args.output))[0]

    svg_files = []

    plots = [
        {
            "x": radii,
            "y": rps,
            "color": ages,
            "color_range": age_range,
            "colorscale": "viridis",
            "colorbar_title": "Age",
            "xlabel": "Solar Radius (R☉)",
            "ylabel": "Rotational Period",
            "title": "Radius vs Rotational Period (color=age)",
            "filename": f"{base_name}_plot1.html",
        },
        {
            "x": ages,
            "y": rps,
            "color": radii,
            "color_range": radius_range,
            "colorscale": "plasma",
            "colorbar_title": "Radius (R☉)",
            "xlabel": "Age",
            "ylabel": "Rotational Period",
            "title": "Rotational Period vs Age (color=radius)",
            "filename": f"{base_name}_plot2.html",
        },
        {
            "x": ages,
            "y": radii,
            "color": rps,
            "color_range": rp_range,
            "colorscale": "viridis",
            "colorbar_title": "Rotational Period",
            "xlabel": "Age",
            "ylabel": "Solar Radius (R☉)",
            "title": "Radius vs Age (color=rotational period)",
            "filename": f"{base_name}_plot3.html",
        },
        {
            "x": ages,
            "y": tangential_velocity,
            "color": masses,
            "color_range": mass_range,
            "colorscale": "plasma",
            "colorbar_title": "Mass (M☉)",
            "xlabel": "Age",
            "ylabel": "Tangential Velocity (km/s)",
            "title": "Tangential Velocity vs Age (color=mass)",
            "filename": f"{base_name}_plot4.html",
            "reference_line": 299792,
            "reference_line_label": "Speed of Light",
        },
        {
            "x": ages,
            "y": temperatures,
            "color": radii,
            "color_range": radius_range,
            "colorscale": "plasma",
            "colorbar_title": "Radius (R☉)",
            "xlabel": "Age",
            "ylabel": "Temperature (K)",
            "title": "Temperature vs Age (color=radius)",
            "filename": f"{base_name}_plot5.html",
        },
        {
            "x": masses,
            "y": radii,
            "color": rps,
            "color_range": rp_range,
            "colorscale": "viridis",
            "colorbar_title": "Rotational Period",
            "xlabel": "Mass (M☉)",
            "ylabel": "Radius (R☉)",
            "title": "Radius vs Mass (color=rotational period)",
            "filename": f"{base_name}_plot6.html",
        },
    ]

    for i, plot_config in enumerate(plots):
        print(f"Generating plot {i+1}/{len(plots)}: {plot_config['title']}'")

        svg_content = create_svg_plot(
            plot_config["x"],
            plot_config["y"],
            plot_config["color"],
            plot_config["color_range"],
            plot_config["colorscale"],
            systems,
            bodies,
            ages,
            radii,
            rps,
            plot_config["xlabel"],
            plot_config["ylabel"],
            plot_config["title"],
            plot_config["colorbar_title"],
            reference_line=plot_config.get("reference_line"),
            reference_line_label=plot_config.get("reference_line_label"),
        )

        svg_path = os.path.join(output_dir, plot_config["filename"])
        # Write atomically: write to a temporary file, flush+fsync, then replace
        tmp_path = svg_path + ".tmp"
        try:
            with open(tmp_path, "w", encoding="utf-8") as f:
                f.write(svg_content)
                f.flush()
                try:
                    os.fsync(f.fileno())
                except Exception:
                    # Not fatal on platforms without fsync support for this descriptor
                    pass
            os.replace(tmp_path, svg_path)
        except Exception as e:
            # If atomic write fails (permissions, fs issues, interrupted), cleanup and fallback
            try:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
            except Exception:
                pass
            try:
                with open(svg_path, "w", encoding="utf-8") as f:
                    f.write(svg_content)
            except Exception as e2:
                print(f"Failed to write SVG to {svg_path}: {e2}", flush=True)
        svg_files.append(plot_config["filename"])
        print(f"  Saved {svg_path}")

    print("Generating HTML with tabs...")
    create_html_with_tabs(svg_files, args.output, len(pts))

    print(f"Saved {len(pts)} points and plots to {args.output}")
    print(f"SVG files: {', '.join(svg_files)}")


if __name__ == "__main__":
    main()
