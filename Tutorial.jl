using Subzero  # bring Subzero into scope
using CairoMakie, GeoInterfaceMakie # bring plotting packages into scope
import GeoInterface as GI

grid = RegRectilinearGrid(; x0 = -1e5, xf = 1e5, y0 = 0.0, yf = 1e5, Nx = 20, Ny = 10)

fig = Figure();
ax1 = Axis(fig[1, 1];  # set up axis tick marks to match grid cells
    title = "Grid Setup",
    xticks = range(grid.x0, grid.xf, 5), xminorticks = IntervalsBetween(5),
    xminorticksvisible = true, xminorgridvisible = true, xticklabelrotation = pi/4,
    yticks = range(grid.y0, grid.yf, 3), yminorticks = IntervalsBetween(5),
    yminorticksvisible = true, yminorgridvisible = true,
)#using vim now
lines!(  # plot boundary of grid with a dashed line
    [grid.x0, grid.x0, grid.xf, grid.xf, grid.x0],  # point x-values
    [grid.y0, grid.yf, grid.yf, grid.y0, grid.y0];  # point y-values
    linestyle = :dash, linewidth = 3.0)

colsize!(fig.layout, 1, Aspect(1, 2))
resize_to_layout!(fig)
fig  # display the figure

north_bound = CollisionBoundary(North; grid)
south_bound = CollisionBoundary(South; grid)
east_bound = PeriodicBoundary(East; grid)
west_bound = PeriodicBoundary(West; grid)

poly!(   # plot each of the boundaries with a 50% transparent color so we can see the overlap
    [north_bound.poly, south_bound.poly, east_bound.poly, west_bound.poly];
    color = [(:purple, 0.5), (:purple, 0.5), (:teal, 0.5), (:teal, 0.5)],
)
ax1.title = "Grid and Boundary Setup"
fig  # display the figure

island1 = GI.Polygon([[(-6e4, 7.5e4), (-4e4, 5e4), (-2.5e4, 7e4), (-6e4, 7.5e4)]])
island2 = GI.Polygon([[(5e4, 2.5e4), (5e4, 5.5e4), (7.5e4, 3e4), (5e4, 2.5e4)]])
topo_list = [island1, island2]

topo_field = initialize_topography_field(; polys = topo_list)

topo_color = RGBf(115/255, 93/255, 55/255)  # brown color for topography
poly!(topo_field.poly; color = topo_color) # plot the topography
ax1.title = "Grid and Domain Setup"
fig  # display the figure

domain = Domain(; north = north_bound, south = south_bound, east = east_bound, west = west_bound, topography = topo_field)

function shear_flow(Nx, Ny, min_u, max_u)
    increasing = true
    curr_u = min_u
    Δu = fld(Nx, 2)  # divide Nx by 2 and round down
    u_vals = zeros(Nx, Ny)
    for (i, col) in enumerate(eachcol(u_vals))
        (i == 1 || i == Ny) && continue  # edges already set to 0m/s
        col .= curr_u
        if increasing && curr_u == max_u  # reach max value and start decreasing velocity
            increasing = false
        end
        if increasing  # update velocity for next column
            curr_u += Δu
        else
            curr_u -= Δu
        end
    end
    return curr_u
end

u_vals = shear_flow(grid.Nx + 1, grid.Ny + 1, 0.0, 0.25)
ocean = Ocean(; u = u_vals, v = 0.0, temp = -1.0, grid)

ax2 = Axis(fig[2, 1]; title = "Ocean U-Velocities [m/s]", xticklabelrotation = pi/4)
xs = grid.x0:grid.Δx:grid.xf
ys = grid.y0:grid.Δy:grid.yf
u_hm = heatmap!(ax2, xs, ys, ocean.u)
Colorbar(fig[2, end+1], u_hm)
resize_to_layout!(fig)
fig

