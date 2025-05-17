using Subzero  # bring Subzero into scope
using CairoMakie, GeoInterfaceMakie # bring plotting packages into scope
import GeoInterface as GI
import Random.Xoshiro

const FT = Float64
grid = RegRectilinearGrid(; x0 = -1e5, xf = 1e5, y0 = 0.0, yf = 1e5, Nx = 20, Ny = 10)

fig = Figure();
ax1 = Axis(fig[1, 1];  # set up axis tick marks to match grid cells
    title = "Grid Setup",
    xticks = range(grid.x0, grid.xf, 5), xminorticks = IntervalsBetween(5),
    xminorticksvisible = true, xminorgridvisible = true, xticklabelrotation = pi/4,
    yticks = range(grid.y0, grid.yf, 3), yminorticks = IntervalsBetween(5),
    yminorticksvisible = true, yminorgridvisible = true,
)
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

atmos = Atmos(; u=-0.3, v=0.2, temp=0.0, grid)

floe_settings = FloeSettings(
  min_floe_area = 1e5,
  max_floe_height = 5,
  subfloe_point_generator = SubGridPointsGenerator(grid, 2),
  stress_calculator = DecayAreaScaledCalculator(0.0, 0.5),
)

coords = [[[3e4, 1e4], [3e4, 1.5e4], [3.5e4, 1.5e4], [3.5e4, 1e4], [3e4, 1e4]]]
hmean = 0.25
Δh = 0.1
floe = Floe(
    FT,
    coords,
    hmean,  # Floe height will be between 0.15 - 0.35
    Δh;  # Δh is the maximum difference between hmean and actual floe height
    floe_settings = floe_settings,
    u = 0.0,
    v = 0.0,
    ξ = 0.0,
    rng = Xoshiro(1), # seed of 1
)

floe1 = [[[6e4, 2e4], [6e4, 5e4], [9e4, 5e4], [9e4, 2e4], [6e4, 2e4]]]
floe2 = [[[5.5e4, 2e4], [5.25e4, 4e4], [5.75e4, 4e4], [5.5e4, 2e4]]]
floe_field = initialize_floe_field(
  FT,
  [floe1, floe2],
  domain,
  0.25,  # mean height of 0.25
  0.0;  # all floes will be the same height
  rng = Xoshiro(1),
  floe_settings = floe_settings,
)

floe_arr = initialize_floe_field(
    FT,
    100,  # attempt to initialize 100 floes
    [1.0; 0.0],  # the top half of the domain is fully packed and the bottom has no floes
    domain,
    0.25,  # mean height of 0.25
    0.10;  # floe heights will range from 0.15-0.35
    floe_settings = floe_settings,
    rng = Xoshiro(1),
)

model = Model(grid, ocean, atmos, domain, floe_arr)

#####################
# CREATE SIMULATION #
#####################

consts = Constants(μ = 0.0) #redefine constants that need to be changed

couple_settings = CouplingSettings(
  coupling_on = true,
  Δt = 10,
  Δd = 1,
  two_way_coupling_on = false,
)

collision_settings = CollisionSettings(
  FT,
  collisions_on = true,
  floe_floe_max_overlap = 0.55,
  floe_domain_max_overlap = 0.75,
)

fracture_settings = FractureSettings(
  fractures_on = false,
  criteria = NoFracture(),
  Δt = 0,
  deform_on = false,
  npieces = 3,
)

pstar = 2.25e5 #Thermodynamic Ice Model stuff -> can be tuned for fracture criteria
c = 20.0
criteria = HiblerYieldCurve(FT, floe_arr, pstar, c)
Δtout = 5

simp_settings = SimplificationSettings(; smooth_vertices_on = true, max_vertices = 30, Δt_smooth = 20)

ridgeraft_settings = RidgeRaftSettings(
  ridge_raft_on = false,
  Δt = 0,
  ridge_probability = 0.95,
  raft_probability = 0.95,
  min_overlap_frac = 0.01,
  min_ridge_height = 0.2,
  max_floe_ridge_height = 5.0,
  max_domain_ridge_height = 1.25,
  max_floe_raft_height = 0.25,
  max_domain_raft_height = 0.25,
  domain_gain_probability = 1.0,
)

weld_settings = WeldSettings(
    weld_on = false,
    Δts = Vector{Int}(),
    Nxs = Vector{Int}(),
    Nys = Vector{Int}(),
    min_weld_area = 1e6,
    max_weld_area = 2e9,
    welding_coeff = 150,
)

initwriter1 = InitialStateOutputWriter(
    dir = ".",
    filename = "initial_state.jld2",
    overwrite = true,
    jld2_kw = Dict{Symbol, Any}(),
)

checkpointer1 = CheckpointOutputWriter(
    Δtout;
    dir = ".",
    filename = "checkpoint.jld2",
    overwrite = true,
    jld2_kw = Dict{Symbol, Any}(),
)

# This thing is what contains info for plotting
floewriter1 = FloeOutputWriter(
    Δtout;
    outputs = collect(fieldnames(Floe)),
    dir = ".",
    filename = "floes.jld2",
    overwrite = true,
    jld2_kw = Dict{Symbol, Any}(),
)

using StructArrays, Subzero

outputwriters2 = OutputWriters(
   initwriter1,
   checkpointer1,
   floewriter1
)

my_simulation = Simulation(
    model = model,
    consts = consts,
    Δt = 5,  # Timestep of 5 seconds
    coupling_settings = couple_settings,
    simp_settings = simp_settings,
    writers = outputwriters2,
)

run!(my_simulation)

timestep_sim!(
   my_simulation,
   tstep,
)

plot_sim(
    floe_output_writer_file_path = ".",
    initial_state_output_writer_file_path = ".",
    Δt = 5,
    output_file_path = "."
)