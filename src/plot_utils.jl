default(fontfamily="Computer Modern", titlefont = (10, "Computer Modern"), legendfontsize = 10, guidefont = (10, :black), tickfont = (10, :black),  yminorgrid = false, grid=false, legend = :none, framestyle = :box)

function plot_kwargs(;kwargs...)
    return Dict(
        :fontfamily=>"Computer Modern", 
        :titlefont => (10, "Computer Modern"), 
        :legendfontsize => 10,
        :guidefont => (10, :black), 
        :tickfont => (10, :black),  
        :yminorgrid => false, 
        :grid=>false, 
        :legend => :none, 
        :framestyle => :box)
end

function nice_points(color::String) 
    return Dict(
        :msc => color, 
        :mc => get_alpha_hex(color,0.65))
end

function hex_to_rgb(value ;full=false) 
    """Convert a hex color to rgb tuple."""
    value = lstrip(value,'#')
    lv = length(value)
    step = lv ÷ 3 
    scale = 1.0 / 255.0
    if full 
        scale = 1.0
    end
    col = Tuple(scale*parse(Int,value[i+begin:i+step],base=16) for i = 0:step:lv-1)
 
    return col 
end

function hex_to_rgb(value, transmit; full=false) 
    """Convert a hex color to rgb tuple."""
    value = lstrip(value,'#')
    lv = len(value)
    step = lv ÷ 3 
    scale = 1.0 / 255.0
    if full 
        scale = 1.0
    end
    col = Tuple(scale*parse(Int,value[i+begin:i+step],base=16) for i = 0:step:lv-1)
 
    return col + (transmit,) 
end

function rgb_to_hex(value) 
    """Convert a rgb tuple to a hex color string."""
    if value[1] < 1 
        scale = 255.0
    else 
        scale = 1.0
    end
    rgb = [round(Int64,scale*k) for k in value]
    return @sprintf "#%02x%02x%02x" rgb[1] rgb[2] rgb[3]  
end

function get_alpha_hex(value,alpha) 
    """Convert a hex color to an equivalent non-transparent version."""

    #first we get the rgb
    rgb = hex_to_rgb(value)

    # apply the transparency
    target = [alpha*k + (0.999-alpha) for k in rgb] 

    return rgb_to_hex(target)
end

function get_colors() 

    cols = [
    "#B39EB5",  # Pastel Purple
    "#C4C3E9",  # Pastel Lavender
    "#AEC6CF",  # Pastel Teal 
    "#89CFF0",  # Pastel Cyan
    "#779ECB",  # Pastel Blue  
    "#FF6961",  # Pastel Red
    "#FF9AA2",  # Pastel Coral
    "#F4A460",  # Pastel Sandy Brown
    "#FFB347",  # Pastel Orange 
    "#D1A677",  # Pastel Brown
    "#D3B3A6",  # Pastel Taupe 
    "#CFCFC4",  # Pastel Magenta 
    "#B2B982",  # Pastel Olive
    "#77DD77",  # Pastel Green
    "#DA70D6",  # Pastel Orchid 
    "#FFB6C1",   # Pastel Light Pink 
]
 
    cols = reverse(cols)
    #circshift!(cols,9)
    return cols 
end
function get_color_lookup(Us) 
    colors = get_colors()
    color_lookup = Dict{Float64,String}()
    for (U,col) in zip(Us,colors) 
        color_lookup[U] = col
    end
    if 3.225 in keys(color_lookup)
        color_lookup[3.2255] = color_lookup[3.225]
    end
    return color_lookup
end
function get_marker_size_lookup() 
    return Dict(
			16 => 2,
			32 => 3.3,
			48 => 3.9,
			64 => 4.4,
			96 => 4.9,
			128 => 5.3
		) 
end