using CSV, DataFrames
using Subzero  # bring Subzero into scope
using CairoMakie, GeoInterfaceMakie # bring plotting packages into scope
import GeoInterface as GI
import GeometryOps as GO
using Random, Statistics

const FT = Float64
const hmean = 0.25
const Δh = 0.0
const Δt = 20

# Open coastline geometry file as a DataFrame
vertices = CSV.read("kobuk_vertices.csv", DataFrame; header=false)
grid = RegRectilinearGrid(FT, x0 = -2e4, xf = 2e4, y0 = -1e4, yf = 4e4, Nx = 100, Ny = 100)

fig = Figure()
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

colsize!(fig.layout, 1, Aspect(1, 1))
resize_to_layout!(fig)
fig  # display the figure

north_bound = CollisionBoundary(North; grid)
south_bound = CollisionBoundary(South; grid)
east_bound = OpenBoundary(East; grid)
west_bound = OpenBoundary(West; grid)

x = vertices[:,1]
y = vertices[:,2]

# outer ring (closed)
outer_ring = GI.LinearRing([(grid.x0, grid.y0), (grid.x0, grid.yf), (grid.xf , grid.yf), (grid.xf , grid.y0), (grid.x0 , grid.y0)])

# coastline as inner ring (must also be a closed LineString)
coastline_points = collect(zip(x, y))
if coastline_points[begin] != coastline_points[end]
    push!(coastline_points, coastline_points[begin])  # close the loop
end

lake_ring = GI.LinearRing(coastline_points)

lake_poly = GI.Polygon(lake_ring)
coastline_poly = GI.Polygon(outer_ring)

diff_poly = GO.difference(coastline_poly, lake_poly; target = GI.PolygonTrait())
println(diff_poly)

poly!(ax1, diff_poly, color=:lightblue, strokecolor=:black, strokewidth=1)

#####################
# CREATE TOPO FIELD #
#####################

coastline_coords = GI.coordinates.(diff_poly)
topo = initialize_topography_field(FT; coords = coastline_coords)

topo_color = RGBf(115/255, 93/255, 55/255)  # brown color for topography
poly!(topo.poly; color = topo_color) # plot the topography
ax1.title = "Grid and Domain Setup"
fig  # display the figure
save("coastline.png", fig)

domain = Domain(; north = north_bound, south = south_bound, east = east_bound, west = west_bound, topography = topo)

########################
# OCEAN & ATMOST SETUP #
########################

ocean = Ocean(;
    u = 0,
    grid,
    v = 0,
    temp = 0,
)
atmos = Atmos(; grid, u = -1.0, v = 0.5, temp = -1.0)

########################
# FLOES & MODEL PARAMS #
########################

coupling_settings = CouplingSettings(
    two_way_coupling_on = true,
)
floe_settings = FloeSettings(
    subfloe_point_generator = SubGridPointsGenerator(grid, 2),
)

floe_arr = initialize_floe_field(
    FT,
    500,
    [0.8],
    domain,
    hmean,
    Δh;
    rng = Random.Xoshiro(1),
    floe_settings = floe_settings
)

model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)

# Run simulation
run_time!(simulation) =  @time run!(simulation)
dir = "output/kobuk"

# Output setup
initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
writers = OutputWriters(initwriter, floewriter)
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 5000,
    verbose = true,
    writers = writers,
    rng = Xoshiro(1),
    coupling_settings = coupling_settings,
)
run_time!(simulation)
plot_sim(
    "output/kobuk/floes.jld2",
    "output/kobuk/initial_state.jld2",
    20,
    "output/kobuk/wind_forcing.mp4",
)