"""
    generate_ground_track_svg(lats, lons, filename)

Generates a standalone SVG file visualizing the satellite ground track.
Handles the International Date Line crossing by creating discontinuous paths.

Security:
- Enforces strict type checking (must be Real numbers) to prevent injection.
- Enforces a maximum number of points to prevent SVG generation Denial of Service (DoS).
- Enforces input validation (finite numbers) to ensure valid SVG output.
"""
function generate_ground_track_svg(lats::AbstractVector{<:Real}, lons::AbstractVector{<:Real}, filename::AbstractString)
    if length(lats) != length(lons)
        throw(ArgumentError("Vectors lats and lons must have the same length. Got $(length(lats)) and $(length(lons))."))
    end

    # Security: Enforce .svg extension to prevent arbitrary file write vulnerabilities
    if !endswith(lowercase(filename), ".svg")
        throw(ArgumentError("Output filename must end with .svg to ensure security and correct file type. Got: $filename"))
    end

    # Security: Limit maximum points to prevent huge SVG generation (DoS)
    # 1,000,000 points creates roughly 30-50MB SVG text, which is already heavy for browsers but manageable.
    MAX_SVG_POINTS = 1_000_000
    if length(lats) > MAX_SVG_POINTS
        throw(ArgumentError("Too many points for SVG generation ($(length(lats)) > $MAX_SVG_POINTS). Please downsample the data to prevent memory exhaustion."))
    end

    # Security: Ensure valid coordinates (no NaNs or Infs)
    if !all(isfinite, lats) || !all(isfinite, lons)
        throw(ArgumentError("Coordinates must be finite real numbers (no NaNs or Infs)."))
    end

    # Security: Enforce .svg extension to prevent arbitrary file writes (e.g., overriding system files)
    if !endswith(lowercase(filename), ".svg")
        throw(ArgumentError("Filename must have .svg extension. Got: $filename"))
    end

    # Security: Prevent path traversal
    if any(x -> x == "..", splitpath(filename)) || isabspath(filename)
        throw(ArgumentError("Security violation: Path traversal and absolute paths are not allowed. Filename must be relative and strictly within the current directory tree."))
    end

    # Simple Equirectangular projection
    # Lon: -180 to 180 -> X: 0 to width
    # Lat: -90 to 90   -> Y: height to 0 (SVG Y is down)
    
    width = 1000
    height = 500
    
    function scale_x(lon)
        return (lon + 180.0) / 360.0 * width
    end
    
    function scale_y(lat)
        return (1.0 - (lat + 90.0) / 180.0) * height
    end
    
    # Build SVG Path(s)
    path_d = IOBuffer()
    markers = IOBuffer()

    if length(lats) > 0
        x0 = scale_x(lons[1])
        y0 = scale_y(lats[1])
        # OPTIMIZATION: Use multi-argument print instead of string interpolation
        # "M $x0 $y0 " allocates a new string for every point.
        # print(io, "M ", x0, " ", y0, " ") writes directly to the buffer.
        print(path_d, "M ", x0, " ", y0, " ")
        
        # Start Marker (Green Circle)
        start_lat = round(lats[1], digits=2)
        start_lon = round(lons[1], digits=2)
        print(markers, """<circle cx="$x0" cy="$y0" r="4" fill="#2ECC40" stroke="#fff" stroke-width="1" class="marker" tabindex="0"><title>Start Point (Lat: $start_lat, Lon: $start_lon)</title></circle>""")

        for i in 2:length(lats)
            lat = lats[i]
            lon = lons[i]
            prev_lon = lons[i-1]
            
            # Check for Date Line crossing (large jump in Longitude)
            if abs(lon - prev_lon) > 300.0
                # Start new path segment
                x = scale_x(lon)
                y = scale_y(lat)
                # OPTIMIZATION: Avoid string interpolation
                print(path_d, "M ", x, " ", y, " ")
            else
                x = scale_x(lon)
                y = scale_y(lat)
                # OPTIMIZATION: Avoid string interpolation
                print(path_d, "L ", x, " ", y, " ")
            end
        end

        # End Marker (Red Square)
        xe = scale_x(lons[end])
        ye = scale_y(lats[end])
        # Center the square (8x8) at (xe, ye)
        rect_x = xe - 4
        rect_y = ye - 4
        end_lat = round(lats[end], digits=2)
        end_lon = round(lons[end], digits=2)
        print(markers, """<rect x="$rect_x" y="$rect_y" width="8" height="8" fill="#FF4136" stroke="#fff" stroke-width="1" class="marker" tabindex="0"><title>End Point (Lat: $end_lat, Lon: $end_lon)</title></rect>""")
    end
    
    path_d_str = String(take!(path_d))
    markers_str = String(take!(markers))

    # Standalone SVG Content
    svg_content = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
    <svg width="$width" height="$height" viewBox="0 0 $width $height" xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="svg-title svg-desc">
        <title id="svg-title">Satellite Ground Track</title>
        <desc id="svg-desc">A map displaying the satellite's path over the Earth. The path starts at the green circle and ends at the red square.</desc>
        <style>
            .background { fill: #001f3f; } /* Deep Ocean Blue */
            .track { fill: none; stroke: #FFDC00; stroke-width: 2; stroke-opacity: 0.8; transition: stroke-width 0.3s, stroke-opacity 0.3s; }
            .track:hover, .track:focus { stroke-width: 4; stroke-opacity: 1.0; outline: none; }
            .marker { transition: stroke-width 0.3s; }
            .marker:focus { stroke-width: 3; outline: none; }
            .grid { stroke: #333; stroke-width: 1; stroke-dasharray: 4; }
            .axis-label { fill: #DDD; font-family: sans-serif; font-size: 12px; }
            .title { fill: #eee; font-family: sans-serif; font-size: 16px; text-anchor: middle; }
        </style>
        
        <!-- Background -->
        <rect width="100%" height="100%" class="background" aria-hidden="true" />
        
        <!-- Grid Lines -->
        <line x1="0" y1="$(height/2)" x2="$width" y2="$(height/2)" class="grid" aria-hidden="true" /> <!-- Equator -->
        <line x1="$(width/2)" y1="0" x2="$(width/2)" y2="$height" class="grid" aria-hidden="true" /> <!-- Prime Meridian -->
        
        <!-- Ground Track -->
        <path d="$path_d_str" class="track" tabindex="0">
            <title>Satellite Path ($(length(lats)) points)</title>
        </path>
        
        <!-- Markers -->
        $markers_str

        <!-- Labels & Title -->
        <text x="$(width/2)" y="25" class="title">DeepOrbit Ground Track</text>
        <text x="5" y="$(height/2 - 5)" class="axis-label">Equator</text>
        <text x="$(width/2 + 5)" y="$(height - 10)" class="axis-label">Prime Meridian</text>

        <!-- Legend (Palette Improvement) -->
        <g transform="translate(10, $(height - 70))">
            <rect width="80" height="60" fill="#000" fill-opacity="0.3" rx="5" />

            <!-- Start Legend -->
            <circle cx="20" cy="20" r="4" fill="#2ECC40" stroke="#fff" stroke-width="1" />
            <text x="35" y="24" class="axis-label">Start</text>

            <!-- End Legend -->
            <rect x="16" y="36" width="8" height="8" fill="#FF4136" stroke="#fff" stroke-width="1" />
            <text x="35" y="44" class="axis-label">End</text>
        </g>
    </svg>
    """
    
    open(filename, "w") do f
        write(f, svg_content)
    end
end
