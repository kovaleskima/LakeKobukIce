using CSV, DataFrames
using Subzero  # bring Subzero into scope
using CairoMakie, GeoInterfaceMakie # bring plotting packages into scope
import GeoInterface as GI

# Open coastline geometry file as a DataFrame
vertices = CSV.read("kobuk_vertices.csv", DataFrame; header=false)

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
multipoint = GI.MultiPoint(GI.Point.(zip(x,y)));
plot!(ax1, multipoint; marker = '.', markersize = 20)
fig


