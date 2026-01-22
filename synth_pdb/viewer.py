"""
3D Molecular Viewer for synth-pdb.

Opens generated PDB structures in browser-based 3D viewer using 3Dmol.js.
Based on pdbstat's molecular viewer implementation.
"""

import logging
import tempfile
import webbrowser
from pathlib import Path

logger = logging.getLogger(__name__)


def view_structure_in_browser(
    pdb_content: str,
    filename: str = "structure.pdb",
    style: str = "cartoon",
    color: str = "spectrum",
) -> None:
    """
    Open 3D structure viewer in browser.
    
    EDUCATIONAL NOTE - Why Browser-Based Visualization:
    Browser-based viewers are ideal for quick structure inspection because:
    1. No installation required (works on any system with a browser)
    2. Interactive (rotate, zoom, change styles)
    3. Lightweight (uses 3Dmol.js JavaScript library)
    4. Shareable (can save HTML file and share with others)
    
    Args:
        pdb_content: PDB file contents as string
        filename: Name to display in viewer title
        style: Initial representation style (cartoon/stick/sphere/line)
        color: Initial color scheme (spectrum/chain/ss/white)
        
    Raises:
        Exception: If viewer fails to open
        
    Example:
        >>> pdb = generate_pdb_content(length=20)
        >>> view_structure_in_browser(pdb, "my_peptide.pdb")
    """
    try:
        logger.info(f"Opening 3D viewer for {filename}")
        
        # Generate HTML with embedded 3Dmol.js viewer
        html = _create_3dmol_html(pdb_content, filename, style, color)
        
        # Save to temporary file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".html", delete=False, encoding="utf-8"
        ) as f:
            f.write(html)
            temp_path = f.name
        
        # Open in default browser
        webbrowser.open(f"file://{temp_path}")
        
        logger.info(f"3D viewer opened in browser: {temp_path}")
        
    except Exception as e:
        logger.error(f"Failed to open 3D viewer: {e}")
        raise


def _create_3dmol_html(
    pdb_data: str, filename: str, style: str, color: str
) -> str:
    """
    Generate HTML with embedded 3Dmol.js viewer.
    
    EDUCATIONAL NOTE - 3Dmol.js:
    3Dmol.js is a JavaScript library for molecular visualization that:
    - Runs entirely in the browser (no server needed)
    - Supports PDB, SDF, MOL2, and other formats
    - Provides interactive controls (rotate, zoom, style changes)
    - Uses WebGL for hardware-accelerated 3D graphics
    
    Args:
        pdb_data: PDB file contents
        filename: Name of PDB file for display
        style: Representation style (cartoon/stick/sphere/line)
        color: Color scheme (spectrum/chain/ss/white)
        
    Returns:
        Complete HTML document as string
    """
    # Escape PDB data for JavaScript (prevent injection attacks)
    pdb_escaped = pdb_data.replace("\\", "\\\\").replace("`", "\\`")
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>3D Viewer - {filename}</title>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <style>
        body {{
            margin: 0;
            padding: 0;
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            background: #f5f5f5;
        }}
        #header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            text-align: center;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        #header h1 {{
            margin: 0 0 10px 0;
            font-size: 28px;
            font-weight: 600;
        }}
        #header p {{
            margin: 0;
            opacity: 0.9;
            font-size: 14px;
        }}
        #controls {{
            background: white;
            padding: 20px;
            border-bottom: 1px solid #e0e0e0;
            display: flex;
            gap: 30px;
            align-items: center;
            flex-wrap: wrap;
            box-shadow: 0 2px 5px rgba(0,0,0,0.05);
        }}
        .control-group {{
            display: flex;
            gap: 10px;
            align-items: center;
        }}
        label {{
            font-weight: 600;
            color: #333;
            font-size: 14px;
        }}
        button {{
            padding: 10px 18px;
            background: #667eea;
            color: white;
            border: none;
            border-radius: 6px;
            cursor: pointer;
            font-size: 14px;
            font-weight: 500;
            transition: all 0.2s;
            box-shadow: 0 2px 4px rgba(102, 126, 234, 0.2);
        }}
        button:hover {{
            background: #5568d3;
            transform: translateY(-1px);
            box-shadow: 0 4px 8px rgba(102, 126, 234, 0.3);
        }}
        button.active {{
            background: #10b981;
            box-shadow: 0 2px 4px rgba(16, 185, 129, 0.2);
        }}
        button.active:hover {{
            background: #059669;
        }}
        #viewer {{
            width: 100%;
            height: calc(100vh - 200px);
            position: relative;
            background: white;
        }}
        #instructions {{
            background: #f9fafb;
            padding: 12px;
            text-align: center;
            color: #6b7280;
            font-size: 13px;
            border-top: 1px solid #e0e0e0;
        }}
        .emoji {{
            font-size: 16px;
        }}
    </style>
</head>
<body>
    <div id="header">
        <h1>🧬 synth-pdb 3D Molecular Viewer</h1>
        <p>{filename}</p>
    </div>

    <div id="controls">
        <div class="control-group">
            <label>Style:</label>
            <button id="btn-cartoon" onclick="setStyle('cartoon')">Cartoon</button>
            <button id="btn-stick" onclick="setStyle('stick')">Stick</button>
            <button id="btn-sphere" onclick="setStyle('sphere')">Sphere</button>
            <button id="btn-line" onclick="setStyle('line')">Line</button>
        </div>

        <div class="control-group">
            <label>Color:</label>
            <button id="color-spectrum" onclick="setColor('spectrum')">Spectrum</button>
            <button id="color-chain" onclick="setColor('chain')">Chain</button>
            <button id="color-ss" onclick="setColor('ss')">Secondary Structure</button>
            <button id="color-white" onclick="setColor('white')">White</button>
        </div>

        <div class="control-group">
            <button onclick="resetView()">🔄 Reset View</button>
            <button onclick="toggleSpin()">🔄 Toggle Spin</button>
        </div>
    </div>

    <div id="viewer"></div>

    <div id="instructions">
        <span class="emoji">🖱️</span> Left-click + drag to rotate | Scroll to zoom | Right-click + drag to pan
    </div>

    <script>
        let viewer;
        let currentStyle = '{style}';
        let currentColor = '{color}';
        let spinning = false;

        // Initialize viewer when page loads
        window.addEventListener('load', function() {{
            let element = document.getElementById('viewer');
            let config = {{ backgroundColor: 'white' }};
            viewer = $3Dmol.createViewer(element, config);

            // Load PDB data
            let pdbData = `{pdb_escaped}`;
            viewer.addModel(pdbData, "pdb");

            // Set initial style and render
            applyStyle();
            viewer.zoomTo();
            viewer.render();

            // Highlight active buttons
            updateActiveButtons();
        }});

        function applyStyle() {{
            // Clear existing styles
            viewer.setStyle({{}}, {{}});

            // Apply new style based on current selection
            let styleObj = {{}};
            if (currentStyle === 'cartoon') {{
                styleObj.cartoon = {{ color: currentColor }};
            }} else if (currentStyle === 'stick') {{
                styleObj.stick = {{ colorscheme: currentColor }};
            }} else if (currentStyle === 'sphere') {{
                styleObj.sphere = {{ colorscheme: currentColor }};
            }} else if (currentStyle === 'line') {{
                styleObj.line = {{ colorscheme: currentColor }};
            }}

            viewer.setStyle({{}}, styleObj);
            viewer.render();
        }}

        function setStyle(style) {{
            currentStyle = style;
            applyStyle();
            updateActiveButtons();
        }}

        function setColor(color) {{
            currentColor = color;
            applyStyle();
            updateActiveButtons();
        }}

        function resetView() {{
            viewer.zoomTo();
            viewer.render();
        }}

        function toggleSpin() {{
            if (spinning) {{
                viewer.spin(false);
                spinning = false;
            }} else {{
                viewer.spin(true);
                spinning = true;
            }}
        }}

        function updateActiveButtons() {{
            // Remove all active classes
            document.querySelectorAll('button').forEach(btn => btn.classList.remove('active'));

            // Add active class to current style and color buttons
            let styleBtn = document.getElementById('btn-' + currentStyle);
            if (styleBtn) styleBtn.classList.add('active');

            let colorBtn = document.getElementById('color-' + currentColor);
            if (colorBtn) colorBtn.classList.add('active');
        }}
    </script>
</body>
</html>
"""
    return html
