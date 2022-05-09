using Plots
using Random, Distributions
using LinearAlgebra
using StaticArrays
using ProgressMeter

G = 4.498502151469553e-6

mutable struct Molecule
    #Position
    pos::SVector{3,Float64}
    #Velocity
    vel::SVector{3,Float64}
    #Mass
    mass::Float64
    #Force
    f::SVector{3,Float64}
end
function clear_f(body::Molecule)
    body.f=SA_F64[0, 0, 0]
end
function move_r(body::Molecule,dt::Float64,cte::Float64=1)
    body.pos+=body.vel*dt*cte
end
function move_v(body::Molecule,dt::Float64,cte::Float64=1)
    body.vel+=body.f*dt*cte/body.mass
end
function calc_f(body::Array{Molecule,1}; ϵ::Float64=0.02)
    Threads.@threads for i in 1:length(body)
        for j in i+1:length(body)
            dr = body[j].pos-body[i].pos
            d = norm(dr)
            F = -G*body[i].mass*body[j].mass*dr/((d+ϵ)^3)
            body[i].f+=-F #3 Newton Law
            body[j].f+=F 
        end
    end
end
function step_PEFRL(body::Array{Molecule,1},Δt::Float64)
    const1 = 0.1786178958448091      #ζ
    const3 = -0.6626458266981849e-1  #χ
    const4 = -0.2123418310626054     #λ
    const2 = (1-2*const4)/2          #(1-2λ)/2
    const5 = 1-2*(const3+const1)     #1-2*(ζ+χ)
    Threads.@threads for particle in body
        move_r(particle,Δt,const1)
        clear_f(particle)
    end
    calc_f(body)
    Threads.@threads for particle in body
        move_v(particle,Δt,const2)
        move_r(particle,Δt,const3)
        clear_f(particle)
    end
    calc_f(body)
    Threads.@threads for particle in body
        move_v(particle,Δt,const4)
        move_r(particle,Δt,const5)
        clear_f(particle)
    end
    calc_f(body)
    Threads.@threads for particle in body
        move_v(particle,Δt,const4)
        move_r(particle,Δt,const3)
        clear_f(particle)
    end
    calc_f(body)
    Threads.@threads for particle in body
        move_v(particle,Δt,const2)
        move_r(particle,Δt,const1)
    end
end
function Propagate(particles::Array{Molecule,1}, steps::Int64; Δt::Float64 = 0.1, Nframes=200, show_progress=false)
    frames = Array{Molecule,1}[] 
    draw=round(Int32,steps/Nframes)
    push!(frames, deepcopy(particles)) #save initial condition
    p = Progress(steps, enabled=show_progress, showspeed=true)
    for i in 1:steps
        #Perform time step integration
        step_PEFRL(particles,Δt)
        
        #Save simulation data every draw steps 
        if i%draw==0
            push!(frames, deepcopy(particles))
        end
        next!(p)
    end
    return frames
end

function show_particles(particles, θ::Float64 = 30.0, ϕ::Float64 = 30.0; Lim=1)
    scatter([Tuple(p.pos) for p = particles],
        lims = (-abs(Lim), Lim),
        camera = (θ, ϕ),
        size = (500, 500),
        label = "",
        background_color = :black,
        marker = (:circle, 3, 0.8, :white, stroke(0)))
end
function show_quadview(particles; Lim=1)
    scatter([Tuple(p.pos) for p = particles],
        lims = (-abs(Lim), Lim),
        camera = (30, 30),
        size = (800, 800),
        label = "",
        background_color = :black,
        marker = (:circle, 3, 0.8, :white, stroke(0)),
        layout = 4,
        subplot = 1)
    scatter!([Tuple(p.pos) for p = particles],
        lims = (-abs(Lim), Lim),
        camera = (0, 90),
        label = "",
        marker = (:circle, 3, 0.8, :white, stroke(0)),
        subplot = 2)
    scatter!([Tuple(p.pos) for p = particles],
        lims = (-abs(Lim), Lim),
        camera = (90, 0),
        label = "",
        marker = (:circle, 3, 0.8, :white, stroke(0)),
        subplot = 3)
    scatter!([Tuple(p.pos) for p = particles],
        lims = (-abs(Lim), Lim),
        camera = (60, 30),
        label = "",
        marker = (:circle, 3, 0.8, :white, stroke(0)),
        subplot = 4)
end
function animate_frames(frames; frame_rate::Int64 = 1, file_name = "animation",
                        θ::Float64 = 30.0, ϕ::Float64 = 30.0, quadview::Bool = false, Lim=1)
    anim = Animation()
    for i = 1:frame_rate:length(frames)
        if quadview
            show_quadview(frames[i], Lim=Lim)
        else
            show_particles(frames[i], θ, ϕ, Lim=Lim)
        end
        frame(anim)
    end
    gif(anim, "$file_name.gif")
end