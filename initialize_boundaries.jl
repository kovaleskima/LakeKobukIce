using CSV, DataFrames
using Subzero  # bring Subzero into scope
using CairoMakie, GeoInterfaceMakie # bring plotting packages into scope
import GeoInterface as GI

# Open coastline geometry file as a DataFrame
vertices = CSV.read("kobuk_vertices.csv", DataFrame; header=false)
x0 = -4e4
xf = 4e4
y0 = -4e4
yf = 4e4 
grid = RegRectilinearGrid(; x0 = -4e4, xf = 4e4, y0 = -4e4, yf = 4e4, Nx = 100, Ny = 100)

fig = Figure();
ax1 = Axis(fig[1, 1];  # set up axis tick marks to match grid cells
    title = "Grid Setup",
    xticks = range(grid.x0, grid.xf, 5), xminorticks = IntervalsBetween(5),
    xminorticksvisible = true, xminorgridvisible = true, xticklabelrotation = pi/4,
    yticks = range(grid.y0, grid.yf, 5), yminorticks = IntervalsBetween(5),
    yminorticksvisible = true, yminorgridvisible = true,
)
lines!(  # plot boundary of grid with a dashed line
    [grid.x0, grid.x0, grid.xf, grid.xf, grid.x0],  # point x-values
    [grid.y0, grid.yf, grid.yf, grid.y0, grid.y0];  # point y-values
    linestyle = :dash, linewidth = 3.0)

#colsize!(fig.layout, 1, Aspect(1, 1))
#resize_to_layout!(fig)
fig  # display the figure

north_bound = OpenBoundary(North; grid)
south_bound = CollisionBoundary(South; grid)
east_bound = OpenBoundary(East; grid)
west_bound = CollisionBoundary(West; grid)

#=
poly!(   # plot each of the boundaries with a 50% transparent color so we can see the overlap
    [north_bound.poly, south_bound.poly, east_bound.poly, west_bound.poly];
    color = [(:purple, 0.5), (:purple, 0.5), (:teal, 0.5), (:teal, 0.5)],
)
=#
ax1.title = "Grid and Boundary Setup"
fig  # display the figure

x = vertices[:,1]
y = vertices[:,2]

# outer ring (closed)
outer_ring = GI.LineString([(x0,y0), (x0,yf), (xf,yf), (xf,y0), (x0,y0)])

# coastline as inner ring (must also be a closed LineString)
coastline_points = collect(zip(x, y))
if coastline_points[begin] != coastline_points[end]
    push!(coastline_points, coastline_points[begin])  # close the loop
end
hole_ring = GI.LineString(coastline_points)

topo_features = GI.Polygon([outer_ring, hole_ring])

# Get rings from polygon
rings = GI.getring(topo_features)

# Outer ring
outer_coords = GI.coordinates(rings[1])
outer_x, outer_y = map(x -> getindex.(outer_coords, x), (1, 2))

hole_coords = GI.coordinates(rings[2])
hole_x, hole_y = map(x -> getindex.(hole_coords, x), (1, 2))

# Plot
fig = Figure()
ax = Axis(fig[1, 1], title="Lake Geometry", xlabel="x", ylabel="y")

# Plot outer polygon filled
poly!(ax, outer_x, outer_y; color=:green, strokecolor=:black)

# Plot hole by drawing it in white to simulate cut-out
poly!(ax, hole_x, hole_y; color=:white, strokecolor=:black)

fig

#topo_field = initialize_topography_field([; topo_features])

#topo_color = RGBf(115/255, 93/255, 55,355)
#poly!(topo_field.poly; color = topo_color)
#ax1.title = "Grid and Domain Setup"

