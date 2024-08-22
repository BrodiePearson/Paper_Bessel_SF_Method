# # Forced-dissipative barotropic QG beta-plane turbulence
#
#md # This example can be run online via [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/barotropicqg_betaforced.ipynb).
#md # Also, it can be viewed as a Jupyter notebook via [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/barotropicqg_betaforced.ipynb).
#
# A simulation of forced-dissipative barotropic quasi-geostrophic turbulence on
# a beta plane. The dynamics include linear drag and stochastic excitation.

using FourierFlows, Plots, Statistics, Printf, Random

using FourierFlows: parsevalsum
using FFTW: irfft
using Statistics: mean
using Random: seed!
using MAT

import GeophysicalFlows.BarotropicQG
import GeophysicalFlows.BarotropicQG: energy, enstrophy
import GeophysicalFlows.BarotropicQG: energy_dissipation, energy_work, energy_drag
import GeophysicalFlows.BarotropicQG: enstrophy_dissipation, enstrophy_work, enstrophy_drag


# ## Choosing a device: CPU or GPU

dev = CPU()    # Device (CPU/GPU)ENV["GRDIR"]=""
#Pkg.build("GR")
#nothing # hide


# ## Numerical parameters and time-stepping parameters

     nx = 2048            # 2D resolution = nx^2   [512 for IC, 1024 for simulations]
stepper = "FilteredRK4"  # timestepper
     dt = 0.005           # timestep [0.02]
 nsteps = 300000*2           # total number of time-steps [300000]
 nsubs  = 1000            # number of time-steps for intermediate logging/plotting [1000]
nothing # hide


# ## Physical parameters

Lx = 2π        # domain size   [2π]
 βʸ = 10.0      # planetary PV gradient in y-direction [10.0, 10/sqrt(2), 1.0, 1/sqrt(2), 0.0]
 βˣ = 0.0    # planetary PV gradient in x-direction  [0.0, 10/sqrt(2), 0.0, 1/sqrt(2), 0.0]
 nμ = -1       # Hypo-viscosity order (not built in yet) [-1]
 ε = 1e-5              # energy input rate by the forcing [1e-5]
 μ = 0.044      # bottom drag [0.044]
 ν = 1e-21      # Hyper-viscous coefficient  [1e-21]
 nν = 4        # Power of hyper-viscosity   [4]
nothing # hide


# ## Forcing
#
# We force the vorticity equation with stochastic excitation that is delta-correlated
# in time and while spatially homogeneously and isotropically correlated. The forcing
# has a spectrum with power in a ring in wavenumber space of radious $k_f$ and
# width $\delta k_f$, and it injects energy per unit area and per unit time equal
# to $\varepsilon$.

forcing_wavenumber = 100.0    # the central forcing wavenumber for a spectrum that is a ring in wavenumber space
forcing_bandwidth  = 1.5     # the width of the forcing spectrum

gr  = TwoDGrid(nx, Lx)

k = [ gr.kr[i] for i=1:gr.nkr, j=1:gr.nl] # a 2D grid with the zonal wavenumber

forcing_spectrum = @. exp( -(sqrt(gr.Krsq)-forcing_wavenumber)^2 / (2forcing_bandwidth^2) )
@. forcing_spectrum[ gr.Krsq < (2π/Lx*2)^2  ] = 0
@. forcing_spectrum[ gr.Krsq > (2π/Lx*200)^2 ] = 0
@. forcing_spectrum[ k .< 2π/Lx ] .= 0 # make sure forcing does not have power at k=0 component
ε0 = parsevalsum(forcing_spectrum .* gr.invKrsq/2, gr)/(gr.Lx*gr.Ly)
@. forcing_spectrum = ε/ε0 * forcing_spectrum  # normalization so that forcing injects energy ε per domain area per unit time

seed!(1234) # reset of the random number generator for reproducibility
nothing # hide

# Next we construct function `calcF!` that computes a forcing realization every timestep
function calcFq!(Fh, sol, t, clock, vars, params, grid)
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
prob = BarotropicQG.Problem(dev; nx=nx, Lx=Lx, β=βʸ, βₓ=βˣ, μ=μ, nμ=nμ, dt=dt, stepper=stepper,
                            ν=ν, nν=nν, calcFq=calcFq!, stochastic=true)
nothing # hide

# Let's define some shortcuts.
sol, cl, vs, pr, gr = prob.sol, prob.clock, prob.vars, prob.params, prob.grid
x, y = gr.x, gr.y
nothing # hide


# First let's see how a forcing realization looks like.
calcFq!(vs.Fqh, sol, 0.0, cl, vs, pr, gr)

heatmap(x, y, irfft(vs.Fqh, gr.nx),
     aspectratio = 1,
               c = :balance,
            clim = (-8, 8),
           xlims = (-gr.Lx/2, gr.Lx/2),
           ylims = (-gr.Ly/2, gr.Ly/2),
          xticks = -3:3,
          yticks = -3:3,
          xlabel = "x",
          ylabel = "y",
           title = "a forcing realization",
      framestyle = :box)


# ## Setting initial conditions

# Use intial conditions simulation to set this up.
file = matopen("simulations/Output/Data/Spin-up/Anisotropic2D_Initial_Condition_LM.mat")
ω = read(file, "zeta")
ωᵢ = zeros(nx, nx)
ωᵢ[1:4:end-3,1:4:end-3] = ω # Interpolate onto the higher resolution grid
ωᵢ[3:4:end-5,:] = (ωᵢ[1:4:end-7,:]+ωᵢ[5:4:end-3,:])/2
ωᵢ[end-1,:] = (ωᵢ[1,:]+ωᵢ[end-3,:])/2
ωᵢ[2:2:end-2,:] = (ωᵢ[1:2:end-3,:]+ωᵢ[3:2:end-1,:])/2
ωᵢ[end,:] = (ωᵢ[1,:]+ωᵢ[end-1,:])/2
ωᵢ[:,3:4:end-5] = (ωᵢ[:,1:4:end-7]+ωᵢ[:,5:4:end-3])/2
ωᵢ[:,end-1] = (ωᵢ[:,1]+ωᵢ[:,end-3])/2
ωᵢ[:,2:2:end-2] = (ωᵢ[:,1:2:end-3]+ωᵢ[:,3:2:end-1])/2
ωᵢ[:,end] = (ωᵢ[:,1]+ωᵢ[:,end-1])/2

BarotropicQG.set_zeta!(prob, ωᵢ)


# ## Diagnostics

# Create Diagnostic -- `energy` and `enstrophy` are functions imported at the top.
E = Diagnostic(energy, prob; nsteps=nsteps)
R = Diagnostic(energy_drag,        prob; nsteps=nsteps) # dissipation by drag
D = Diagnostic(energy_dissipation, prob; nsteps=nsteps) # dissipation by hyperviscosity
W = Diagnostic(energy_work,        prob; nsteps=nsteps) # work input of energy
Z = Diagnostic(enstrophy, prob; nsteps=nsteps)
RZ = Diagnostic(enstrophy_drag,        prob, nsteps=nsteps) # enstrophy dissipation by drag
WZ = Diagnostic(enstrophy_work,        prob, nsteps=nsteps) # work input of enstrophy
DZ = Diagnostic(enstrophy_dissipation,        prob, nsteps=nsteps) # enstrophy dissipation by visc
diags = [E, D, W, R, Z, RZ, WZ, DZ] # A list of Diagnostics types passed to "stepforward!" will  be updated every timestep.
nothing # hidenothing # hide


# ## Output

# We choose folder for outputing `.jld2` files and snapshots (`.png` files).
# Define base filename so saved data can be distinguished from other runs
base_filename = string("Anisotropic2D_n_", nx, "_drag_", round(μ, sigdigits=1), "_order_", 2*nμ, "_visc_", round(ν, sigdigits=1), "_order_", 2*nν, "_kf_", forcing_wavenumber, "_F_", ε, "_betay_", βʸ, "_betax_", βˣ)
# We choose folder for outputing `.jld2` files and snapshots (`.png` files).
datapath = "./simulations/Output/Data/Equilibrated/"
plotpath = "./simulations/Output/Figures/Equilibrated/"
plotname = "snapshots"

dataname = joinpath(datapath, base_filename)
plotname = joinpath(plotpath, base_filename)
jld2_filename = joinpath(datapath, string(base_filename, ".jld2"))
matfilename= joinpath(datapath, base_filename)
nothing # hide

# Do some basic file management,
if isfile(jld2_filename); rm(jld2_filename); end
if isfile(matfilename); rm(matfilename); end
#if !isdir(plotpath); mkdir(plotpath); end
nothing # hide

# and then create Output.
get_sol(prob) = sol # extracts the Fourier-transformed solution
get_u(prob) = irfft(im*gr.l.*gr.invKrsq.*sol, gr.nx)
out = Output(prob, dataname, (:sol, get_sol), (:u, get_u))
nothing # hide

function write_matlab_file(writepath, vars, grid, time, diags)
    E, D, W, R, Z = diags
    step_path = string(writepath, "_", round(Int, cl.t), ".mat")
    matfile = matopen(step_path, "w")
    write(matfile, "zeta", vars.zeta)
    write(matfile, "u", vars.u)
    write(matfile, "v", vars.v)
    write(matfile, "nx", grid.nx)
    write(matfile, "ny", grid.ny)
    write(matfile, "dx", grid.dx)
    write(matfile, "dy", grid.dy)
    write(matfile, "Lx", grid.Lx)
    write(matfile, "Ly", grid.Ly)
    write(matfile, "KE", E.data[E.i])
    write(matfile, "ENS", Z.data[Z.i])
    write(matfile, "Work_Rate_KE", W.data[W.i])
    write(matfile, "Diss_Rate_KE", D.data[D.i])
    write(matfile, "Drag_Rate_KE", R.data[R.i])
    write(matfile, "Work_Rate_ENS", WZ.data[WZ.i])
    write(matfile, "Diss_Rate_ENS", DZ.data[DZ.i])
    write(matfile, "Drag_Rate_ENS", RZ.data[RZ.i])
    close(matfile)
end



# ## Visualizing the simulation

# We define a function that plots the vorticity and streamfunction fields, their
# corresponding zonal mean structure and timeseries of energy and enstrophy.

function plot_output(prob)
    ζ = prob.vars.zeta
    ζ̄ = mean(ζ, dims=1)'
    ū = mean(prob.vars.u, dims=1)'
    Energy =  0.5*(prob.vars.u.^2 + prob.vars.v.^2)

  pζ = heatmap(x, y, ζ,
       aspectratio = 1,
            legend = false,
                 c = :balance,
              clim = (-8, 8),
             xlims = (-gr.Lx/2, gr.Lx/2),
             ylims = (-gr.Ly/2, gr.Ly/2),
            xticks = -3:3,
            yticks = -3:3,
            xlabel = "x",
            ylabel = "y",
             title = "vorticity ζ=∂v/∂x-∂u/∂y",
        framestyle = :box)

  pψ = contourf(x, y, Energy,
            levels = -0.32:0.04:0.32,
       aspectratio = 1,
         linewidth = 1,
            legend = false,
              clim = (0, 1e-3),
                 c = :viridis,
             xlims = (-gr.Lx/2, gr.Lx/2),
             ylims = (-gr.Ly/2, gr.Ly/2),
            xticks = -3:3,
            yticks = -3:3,
            xlabel = "x",
            ylabel = "y",
             title = "Energy",
        framestyle = :box)

  pζm = plot(ζ̄, y,
            legend = false,
         linewidth = 2,
             alpha = 0.7,
            yticks = -3:3,
             xlims = (-3, 3),
            xlabel = "zonal mean ζ",
            ylabel = "y")
  plot!(pζm, 0*y, y, linestyle=:dash, linecolor=:black)

  pum = plot(ū, y,
            legend = false,
         linewidth = 2,
             alpha = 0.7,
            yticks = -3:3,
             xlims = (-0.5, 0.5),
            xlabel = "zonal mean u",
            ylabel = "y")
  plot!(pum, 0*y, y, linestyle=:dash, linecolor=:black)

  pE = plot(1,
             label = "energy",
         linewidth = 2,
             alpha = 0.7,
             xlims = (-0.1, 4.1),
             ylims = (0, 5e-1),
            xlabel = "μt")

  pZ = plot(1,
             label = "enstrophy",
         linecolor = :red,
            legend = :bottomright,
         linewidth = 2,
             alpha = 0.7,
             xlims = (-0.1, 4.1),
             ylims = (0, 3),
            xlabel = "μt")

  l = @layout grid(2, 3)
  p = plot(pζ, pζm, pE, pψ, pum, pZ, layout=l, size = (1000, 600), dpi=150)

  return p
end
nothing # hide


# ## Time-stepping the `Problem` forward

# We time-step the `Problem` forward in time.

startwalltime = time()

p = plot_output(prob)

anim = @animate for j=0:Int(nsteps/nsubs)

  cfl = cl.dt*maximum([maximum(vs.u)/gr.dx, maximum(vs.v)/gr.dy])

  #log = @sprintf("step: %04d, t: %d, cfl: %.2f, KE: %.2e, Ens: %.2e, W: %.2e, Diss: %.2e, Drag: %.2e, walltime: %.2f min",
  #cl.step, cl.t, cfl, E.data[E.i], Z.data[Z.i], W.data[W.i], D.data[D.i], R.data[R.i],
  #(time()-startwalltime)/60)
  #log = @sprintf("step: %04d, t: %d, cfl: %.2f, KE: %.2e, Ens: %.2e, walltime: %.2f min",
  #cl.step, cl.t, cfl, E.data[E.i], Z.data[Z.i],
  #(time()-startwalltime)/60)

  log = @sprintf("step: %04d, t: %.1f, cfl: %.3f, walltime: %.2f min",
        cl.step, cl.t, cfl, (time()-startwalltime)/60)

  if j%(1000/nsubs)==0; println(log) end

  log = @sprintf("Energy diagnostics - Energy: %.2e, Drag: %.2e, Diss: %.2e, Work Input: %.2e",
            E.data[E.i], R.data[R.i], D.data[D.i], W.data[W.i])

  if j%(1000/nsubs)==0; println(log) end

  log = @sprintf("Enstrophy diagnostics - Enstrophy: %.2e, Drag: %.2e, Diss: %.2e, Work Input: %.2e",
      Z.data[Z.i], RZ.data[RZ.i], DZ.data[DZ.i], WZ.data[WZ.i])

  if j%(1000/nsubs)==0; println(log) end

  if j%(10000/nsubs)==0;  write_matlab_file(matfilename, vs, gr, cl.t, diags) end

  p[1][1][:z] = Array(vs.zeta)
  p[1][:title] = "vorticity, μt="*@sprintf("%.2f", μ*cl.t)
  p[4][1][:z] = Array(vs.psi)
  p[2][1][:x] = mean(vs.zeta, dims=1)'
  p[5][1][:x] = mean(vs.u, dims=1)'
  push!(p[3][1], μ*E.t[E.i], E.data[E.i])
  push!(p[6][1], μ*Z.t[Z.i], Z.data[Z.i])

  stepforward!(prob, diags, nsubs)
  BarotropicQG.updatevars!(prob)

end

mp4(anim, string(plotname, ".mp4"), fps=10)


i₀ = 1
dEdt_numerical = (E[(i₀+1):E.i] - E[i₀:E.i-1])/cl.dt #numerical first-order approximation of energy tendency
ii = (i₀):E.i-1
ii2 = (i₀+1):E.i

t = E.t[ii]
#dEdt_computed = W[ii2] - D[ii] - R[ii]

# residual = dEdt_computed - dEdt_numerical

l = @layout grid(2, 3)

Eᵏ = 0.5*(vs.u.^2+vs.v.^2) # Update diagnosed Kinetic Energy
L = Lx

pb4 = heatmap(x, y, vs.zeta,
          aspectratio = 1,
          legend = false,
               c = :balance,
            clim = (-maximum(abs.(vs.zeta)), maximum(abs.(vs.zeta))),
           xlims = (-L/2, L/2),
           ylims = (-L/2, L/2),
          xticks = -3:3,
          yticks = -3:3,
          xlabel = "x",
          ylabel = "y",
           title = "∇²ψ(x, y, t="*@sprintf("%.2f", cl.t)*")",
      framestyle = :box)

pb2 = plot(μ*t, [W[ii2] ε.+0*t -D[ii] -R[ii]],
           label = ["work, W" "ensemble mean work, <W>" "dissipation, D" "drag, D=-2μE"],
       linestyle = [:solid :dash :solid :solid],
       linewidth = 2,
           alpha = 0.8,
          xlabel = "μt",
          ylabel = "energy sources and sinks")

          pb3 = plot(μ*t, E[ii],
                     label = ["Energy, E"],
                 linestyle = [:solid :dash :solid :solid],
                 linewidth = 2,
                     alpha = 0.8,
                    xlabel = "μt",
                    ylabel = "Energy evolution")

                    pb1 = heatmap(x, y, Eᵏ,
                              aspectratio = 1,
                              legend = false,
                                   c = :balance,
                                clim = (0, maximum(Eᵏ)),
                               xlims = (-L/2, L/2),
                               ylims = (-L/2, L/2),
                              xticks = -3:3,
                              yticks = -3:3,
                              xlabel = "x",
                              ylabel = "y",
                               title = "Energy(x, y, t="*@sprintf("%.2f", cl.t)*")",
                          framestyle = :box)

          pb5 = plot(μ*t, [WZ[ii2] -DZ[ii] -RZ[ii]],
                     label = ["work, WZ" "dissipation, DZ" "drag, D=-2μZ"],
                 linestyle = [:solid :dash :solid :solid],
                 linewidth = 2,
                     alpha = 0.8,
                    xlabel = "μt",
                    ylabel = "Enstrophy sources and sinks")

                    pb6 = plot(μ*t, Z[ii],
                               label = ["Enstrophy, Z"],
                           linestyle = [:solid :dash :solid :solid],
                           linewidth = 2,
                               alpha = 0.8,
                              xlabel = "μt",
                              ylabel = "Enstrophy evolution")

plot_budgets = plot(pb1, pb2, pb3, pb4, pb5, pb6, layout=l, size = (2400, 1200))

png(plotname)

# Last we save the output.
saveoutput(out)
