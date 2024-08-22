# # Forced-dissipative barotropic QG beta-plane turbulence
#
#md # This example can be run online via [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/barotropicqg_betaforced.ipynb).
#md # Also, it can be viewed as a Jupyter notebook via [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/barotropicqg_betaforced.ipynb).
#
# A simulation of forced-dissipative barotropic quasi-geostrophic turbulence on
# a beta plane. The dynamics include linear drag and stochastic excitation.

using FourierFlows
using Plots
using Statistics
using Printf
using Random

using FourierFlows: parsevalsum
using FFTW: irfft
using Statistics: mean
using Random: seed!
using MAT

import GeophysicalFlows.SurfaceQG
import GeophysicalFlows.SurfaceQG: kinetic_energy, buoyancy_variance
import GeophysicalFlows.SurfaceQG: buoyancy_dissipation_hyperviscosity, buoyancy_work
import GeophysicalFlows.SurfaceQG: buoyancy_dissipation_hypoviscosity

# ## Choosing a device: CPU or GPU

dev = CPU()    # Device (CPU/GPU)ENV["GRDIR"]=""
#Pkg.build("GR")
#nothing # hide


# ## Numerical parameters and time-stepping parameters

     nx = 512            # 2D resolution = nx^2   [512 for IC, 1024 for simulations]
stepper = "FilteredRK4"  # timestepper
     dt = 0.016*((1024/nx))         # timestep [0.02 for large beta, 0.01 for small/no beta]
 nsteps = round( Int, 30000/dt ) #600000           # total number of time-steps [300000 for large beta, 600000 for small/no beta]
 nsubs  = round( Int, 50/dt ) #1000            # number of time-steps for intermediate logging/plotting [1000]
nothing # hide

end_time = nsteps*dt


# ## Physical parameters

Lx = 2π        # domain size   [2π]
 #βʸ = 0.0      # planetary PV gradient in y-direction [10.0, 10/sqrt(2), 1.0, 1/sqrt(2), 0.0]
 #βˣ = 0.0    # planetary PV gradient in x-direction  [0.0, 10/sqrt(2), 0.0, 1/sqrt(2), 0.0]
 nμ = -1       # Hypo-viscosity order [-1]
 ε = 1e-5              # energy input rate by the forcing [1e-5]
 μ = 0.005      # bottom drag [0.044]
 ν = 1e-17 #1e-19      # Hyper-viscous coefficient  [1e-21]
 nν = 4        # Power of hyper-viscosity   [4]
 β = 0.5
nothing # hide


# ## Forcing
#
# We force the buoyancy equation with stochastic excitation that is delta-correlated
# in time and while spatially homogeneously and isotropically correlated. The forcing
# has a spectrum with power in a ring in wavenumber space of radious $k_f$ and
# width $\delta k_f$, and it injects energy per unit area and per unit time equal
# to $\varepsilon$.
#

forcing_wavenumber = 7.0    # the central forcing wavenumber for a spectrum that is a ring in wavenumber space
forcing_bandwidth  = 1.5     # the width of the forcing spectrum

gr  = TwoDGrid(nx, Lx)

k = [ gr.kr[i] for i=1:gr.nkr, j=1:gr.nl] # a 2D grid with the zonal wavenumber

forcing_spectrum = @. exp( -(sqrt(gr.Krsq)-forcing_wavenumber)^2 / (2forcing_bandwidth^2) )
@. forcing_spectrum[ gr.Krsq < (2π/Lx*2)^2  ] = 0
@. forcing_spectrum[ gr.Krsq > (2π/Lx*200)^2 ] = 0
@. forcing_spectrum[ k .< 2π/Lx ] .= 0 # make sure forcing does not have power at k=0 component
# Unlike the other QG models which force vorticity but constrain energy,
# the SQG model forces buoyancy and constrains buoyancy variance.
ε0 = parsevalsum(forcing_spectrum, gr)/(gr.Lx*gr.Ly)
@. forcing_spectrum = ε/ε0 * forcing_spectrum  # normalization so that forcing injects energy ε per domain area per unit time

seed!(1234) # reset of the random number generator for reproducibility
nothing # hide

# Next we construct function `calcF!` that computes a forcing realization every timestep
function calcF!(Fh, sol, t, clock, vars, params, grid)
  ξ = ArrayType(dev)(exp.(2π*im*rand(eltype(grid), size(sol)))/sqrt(clock.dt))
  @. Fh = ξ*sqrt.(forcing_spectrum)
  Fh[abs.(grid.Krsq).==0] .= 0
  nothing
end
nothing # hide


# ## Problem setup]
# We initialize a `Problem` by providing a set of keyword arguments. Not providing
# a viscosity coefficient ν leads to the module's default value: ν=0. In this
# example numerical instability due to accumulation of enstrophy in high wavenumbers
# is taken care with the `FilteredTimestepper` we picked.
prob = SurfaceQG.Problem(dev; nx=nx, Lx=Lx, dt=dt, stepper=stepper, β=β,
                            ν=ν, nν=nν, μ=μ, nμ=nμ, calcF=calcF!, stochastic=true)
nothing # hide

# Let's define some shortcuts.
sol, cl, vs, pr, gr = prob.sol, prob.clock, prob.vars, prob.params, prob.grid
x, y = gr.x, gr.y
nothing # hide


# First let's see how a forcing realization looks like.
calcF!(vs.Fh, sol, 0.0, cl, vs, pr, gr)

heatmap(x, y, irfft(vs.Fh, gr.nx),
     aspectratio = 1,
               c = :balance,
            clim = (-8, 8),
           xlims = (-gr.Lx/2, gr.Lx/2),
           ylims = (-gr.Ly/2, gr.Ly/2),
          xticks = -3:3,
          yticks = -3:3,
          xlabel = "x",
          ylabel = "y",
           title = "a buoyancy forcing realization",
      framestyle = :box)


# ## Setting initial conditions

# ## Output

# We choose folder for outputing `.jld2` files and snapshots (`.png` files).
# Define base filename so saved data can be distinguished from other runs
base_filename = string("SurfaceQG_n_", nx, "_visc_", round(ν, sigdigits=1),
"_order_", 2*nν, "_hypovisc_", round(μ, sigdigits=1), "_order_", 2*nμ, "_kf_",
forcing_wavenumber, "_F_", ε, "_endtime_", end_time, "_beta_", β)
# We choose folder for outputing `.jld2` files and snapshots (`.png` files).
datapath = string("../output/data/", base_filename, "/")
plotpath = string("../output/media/", base_filename, "/")
plotname = "snapshots_hypoviscosity"
#if !isdir(plotpath); mkdir(plotpath);
#if !isdir(datapath); mkdir(datapath);

dataname = joinpath(datapath, base_filename)
plotname = joinpath(plotpath, base_filename)
jld2_filename = joinpath(datapath, string(base_filename, ".jld2"))
matfilename= joinpath(datapath, base_filename)
nothing # hide

# Do some basic file management,
if isfile(jld2_filename); rm(jld2_filename); end
if isfile(matfilename); rm(matfilename); end
if !isdir(plotpath); mkdir(plotpath); end
if !isdir(datapath); mkdir(datapath); end
nothing # hide

# Specify the appropriate initial condition file (lower resolution nx = 512)
initial_condition_filename = string("SurfaceQG_n_512_visc_", round(ν, sigdigits=1),
    "_order_", 2*nν, "_hypovisc_", round(μ, sigdigits=1), "_order_", 2*nμ, "_kf_",
    forcing_wavenumber, "_F_", ε, "_endtime_8000.0_8000.mat")

#file = matopen(string("simulations/Output/Data/Spin-up/", initial_condition_filename))
#buoyancy = read(file, "b")
#bᵢ = zeros(nx, nx)
#bᵢ[1:4:end-3,1:4:end-3] = buoyancy # Interpolate onto the higher resolution grid
#bᵢ[3:4:end-5,:] = (bᵢ[1:4:end-7,:]+bᵢ[5:4:end-3,:])/2
#bᵢ[end-1,:] = (bᵢ[1,:]+bᵢ[end-3,:])/2
#bᵢ[2:2:end-2,:] = (bᵢ[1:2:end-3,:]+bᵢ[3:2:end-1,:])/2
#bᵢ[end,:] = (bᵢ[1,:]+bᵢ[end-1,:])/2
#bᵢ[:,3:4:end-5] = (bᵢ[:,1:4:end-7]+bᵢ[:,5:4:end-3])/2
#bᵢ[:,end-1] = (bᵢ[:,1]+bᵢ[:,end-3])/2
#bᵢ[:,2:2:end-2] = (bᵢ[:,1:2:end-3]+bᵢ[:,3:2:end-1])/2
#bᵢ[:,end] = (bᵢ[:,1]+bᵢ[:,end-1])/2

# Our initial condition is simply fluid at rest.
SurfaceQG.set_b!(prob, zeros(gr.nx, gr.ny))
#SurfaceQG.set_b!(prob, bᵢ)


# ## Diagnostics

# Create Diagnostic -- `energy` and `enstrophy` are functions imported at the top.
#E = Diagnostic(energy, prob; nsteps=nsteps)
#R = Diagnostic(energy_drag,        prob; nsteps=nsteps) # dissipation by drag
#D = Diagnostic(energy_dissipation, prob; nsteps=nsteps) # dissipation by hyperviscosity
#W = Diagnostic(energy_work,        prob; nsteps=nsteps) # work input of energy
#Z = Diagnostic(enstrophy, prob; nsteps=nsteps)
#RZ = Diagnostic(enstrophy_drag,        prob, nsteps=nsteps) # enstrophy dissipation by drag
#WZ = Diagnostic(enstrophy_work,        prob, nsteps=nsteps) # work input of enstrophy
#DZ = Diagnostic(enstrophy_dissipation,        prob, nsteps=nsteps) # enstrophy dissipation by visc
#diags = [E, D, W, R, Z, RZ, WZ, DZ] # A list of Diagnostics types passed to "stepforward!" will  be updated every timestep.
#nothing # hidenothing # hide

bb = Diagnostic(buoyancy_variance, prob; nsteps=nsteps)
Dν = Diagnostic(buoyancy_dissipation_hyperviscosity, prob; nsteps=nsteps) # dissipation by hyperviscosity
Dμ = Diagnostic(buoyancy_dissipation_hypoviscosity, prob; nsteps=nsteps) # dissipation by hypoviscosity
W = Diagnostic(buoyancy_work,        prob; nsteps=nsteps) # work input of energy
KE = Diagnostic(kinetic_energy, prob; nsteps=nsteps)

diags = [bb, Dν, Dμ, W, KE] # A list of Diagnostics types passed to "stepforward!" will  be updated every timestep.
nothing # hidenothing # hide


# and then create Output.
get_sol(prob) = sol # extracts the Fourier-transformed solution
get_u(prob) = irfft(im*gr.l.*gr.invKrsq.*sol, gr.nx)
out = Output(prob, dataname, (:sol, get_sol), (:u, get_u))
nothing # hide

function write_matlab_file(writepath, vars, grid, time, diags)
    bb, Dν, Dμ, W, KE = diags
    step_path = string(writepath, "_", round(Int, cl.t), ".mat")
    matfile = matopen(step_path, "w")
    write(matfile, "b", vars.b)
    write(matfile, "u", vars.u)
    write(matfile, "v", vars.v)
    write(matfile, "nx", grid.nx)
    write(matfile, "ny", grid.ny)
    write(matfile, "dx", grid.dx)
    write(matfile, "dy", grid.dy)
    write(matfile, "Lx", grid.Lx)
    write(matfile, "Ly", grid.Ly)
    write(matfile, "kinetic_energy", KE.data[KE.i])
    write(matfile, "buoyancy_variance", bb.data[bb.i])
    write(matfile, "buoyancy_work", W.data[W.i])
    write(matfile, "buoyancy_dissipation_hyperviscosity", Dν.data[Dν.i])
    write(matfile, "buoyancy_dissipation_hypoviscosity", Dμ.data[Dμ.i])
    close(matfile)
end



# ## Visualizing the simulation

# We define a function that plots the vorticity and streamfunction fields, their
# corresponding zonal mean structure and timeseries of energy and enstrophy.

function plot_output(prob)
    bˢ = prob.vars.b
    b̅ˢ = mean(bˢ, dims=1)'
    ū = mean(prob.vars.u, dims=1)'
    Energy =  0.5*(prob.vars.u.^2 + prob.vars.v.^2)

  pζ = heatmap(x, y, bˢ,
       aspectratio = 1,
            legend = false,
                 c = :balance,
              clim = (-maximum(abs.(bˢ)), maximum(abs.(bˢ))),
             xlims = (-gr.Lx/2, gr.Lx/2),
             ylims = (-gr.Ly/2, gr.Ly/2),
            xticks = -3:3,
            yticks = -3:3,
            xlabel = "x",
            ylabel = "y",
             title = "buoyancy bˢ",
        framestyle = :box)

        pψ = heatmap(x, y, bˢ.*bˢ,
             aspectratio = 1,
                  legend = false,
                       c = :balance,
                    clim = (-0.1, 0.1),
                   xlims = (-gr.Lx/8, gr.Lx/8),
                   ylims = (-gr.Ly/8, gr.Ly/8),
                  xticks = -3:3,
                  yticks = -3:3,
                  xlabel = "x",
                  ylabel = "y",
                   title = "buoyancy variance",
              framestyle = :box)

   pζm = plot(b̅ˢ, y,
             legend = false,
          linewidth = 2,
              alpha = 0.7,
             yticks = -3:3,
              xlims = (-0.25, 0.25),
             xlabel = "zonal mean bˢ",
             ylabel = "y")
   plot!(pζm, 0*y, y, linestyle=:dash, linecolor=:black)

   pum = plot(ū, y,
             legend = false,
          linewidth = 2,
              alpha = 0.7,
             yticks = -3:3,
              xlims = (-0.25, 0.25),
             xlabel = "zonal mean u",
             ylabel = "y")
   plot!(pum, 0*y, y, linestyle=:dash, linecolor=:black)

   pE = plot(1,
              label = "kinetic energy",
          linewidth = 2,
              alpha = 0.7,
              xlims = (-0.1, end_time),
              ylims = (0, 2e-1),
             xlabel = "t")

             pZ = plot(1,
                        label = "buoyancy variance",
                    linecolor = :red,
                       legend = :bottomright,
                    linewidth = 2,
                        alpha = 0.7,
                        xlims = (-0.1, end_time),
                        ylims = (0, 2e-1),
                       xlabel = "t")

             l = @layout grid(2, 3)
             p = plot(pζ, pζm, pE, pψ, pum, pZ, layout=l, size = (1000, 600), dpi=150)

  return p
end
nothing # hide


# ## Time-stepping the `Problem` forward

# We time-step the `Problem` forward in time.

startwalltime = time()

p = plot_output(prob)

anim = @animate for j=0:round( Int, nsteps/nsubs )

  cfl = cl.dt*maximum([maximum(vs.u)/gr.dx, maximum(vs.v)/gr.dy])

  #log = @sprintf("step: %04d, t: %d, cfl: %.2f, KE: %.2e, Ens: %.2e, W: %.2e, Diss: %.2e, Drag: %.2e, walltime: %.2f min",
  #cl.step, cl.t, cfl, E.data[E.i], Z.data[Z.i], W.data[W.i], D.data[D.i], R.data[R.i],
  #(time()-startwalltime)/60)
  #log = @sprintf("step: %04d, t: %d, cfl: %.2f, KE: %.2e, Ens: %.2e, walltime: %.2f min",
  #cl.step, cl.t, cfl, E.data[E.i], Z.data[Z.i],
  #(time()-startwalltime)/60)

  log = @sprintf("step: %04d, t: %.3f, cfl: %.3f, walltime: %.2f min",
        cl.step, cl.t, cfl, (time()-startwalltime)/60)

  if j%(nsubs/nsubs)==0; println(log) end

  log = @sprintf("General Diagnostics - Twice KE: %.2e, Buoy. Variance: %.2e",
            2*KE.data[KE.i], bb.data[bb.i])

  if j%(nsubs/nsubs)==0; println(log) end

  log = @sprintf("b² budget diagnostics - Work: %.2e, Dissipation by ν (%.2e) and by μ (%.2e)",
            W.data[W.i], Dν.data[Dν.i], Dμ.data[Dμ.i])

  if j%(nsubs/nsubs)==0; println(log) end

  if j%(10*nsubs/nsubs)==0;  write_matlab_file(matfilename, vs, gr, cl.t, diags) end

  #p[1][1][:z] = Array(vs.b)
  #p[1][:title] = "buoyancy, t="*@sprintf("%.2f", cl.t)
  #p[4][1][:z] = Array(vs.b.*vs.b)
  #p[2][1][:x] = mean(vs.b, dims=1)'
  #p[5][1][:x] = mean(vs.u, dims=1)'
  #push!(p[3][1], KE.t[KE.i], KE.data[KE.i])
  #push!(p[6][1], bb.t[bb.i], bb.data[bb.i])

  stepforward!(prob, diags, nsubs)
  SurfaceQG.updatevars!(prob)

end

write_matlab_file(matfilename, vs, gr, cl.t, diags)

mp4(anim, string(plotname, ".mp4"), fps=10)

i₀ = 1
dbbdt_numerical = (bb[(i₀+1):bb.i] - bb[i₀:bb.i-1])/cl.dt #numerical first-order approximation of energy tendency
ii = (i₀):bb.i-1
ii2 = (i₀+1):bb.i

t = bb.t[ii]
#dEdt_computed = W[ii2] - D[ii] - R[ii]

# residual = dEdt_computed - dEdt_numerical

l = @layout grid(1, 3)

Eᵏ = 0.5*(vs.u.^2+vs.v.^2) # Update diagnosed Kinetic Energy
L = Lx

pb1 = heatmap(x, y, vs.b',
          aspectratio = 1,
          legend = true,
               c = :balance,
            clim = (-maximum(abs.(vs.b)), maximum(abs.(vs.b))),
           xlims = (-L/2, L/2),
           ylims = (-L/2, L/2),
          xticks = -3:3,
          yticks = -3:3,
          xlabel = "x",
          ylabel = "y",
           title = "∇²ψ(x, y, t="*@sprintf("%.2f", cl.t)*")",
      framestyle = :box)

pb2 = plot(t, [W[ii2] ε.+0*t -Dν[ii] -Dμ[ii]],
           label = ["work, W" "ensemble mean work, <W>" "hyper-viscous dissipation, Dν" "hypo-viscous dissipation, Dμ"],
       linestyle = [:solid :dash :solid :solid],
       linewidth = 2,
           alpha = 0.8,
           ylims = (-5e-5, 5e-5),
          xlabel = "t",
          ylabel = "buoyancy variance sources and sinks")

pb3 = plot(t, bb[ii],
                     label = ["Buoyancy variance"],
                 linestyle = [:solid :dash :solid :solid],
                 linewidth = 2,
                     alpha = 0.8,
                    xlabel = "t",
                    ylabel = "")


plot_budgets = plot(pb1, pb2, pb3, layout=l, size = (2400, 800))

png(plotname)
