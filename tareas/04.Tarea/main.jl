using Random, Distributions
using SpecialFunctions
include("ParallelNormalAlgorith.jl")

Random.seed!(123)
d = Normal(0.0,150.0)
u1 = Uniform(0,1)
u11 = Uniform(-1,1)
u50 = Uniform(1,50)

function graph_vel(particles)
    data_r = zeros(length(particles))
    data_v = zeros(length(particles))
    for i in 1:length(particles)
        data_r[i]= norm(particles[i].pos)
        data_v[i]= norm(particles[i].vel)
    end
    plot = scatter(data_r,data_v)
    xlabel!("r")
    ylabel!("v")
    return plot
end

function rand_kepler_particles2(num_particles::Int64, Black_Hole_Mass, R) #2D galaxy
    u1 = Uniform(0,1)
    u50 = Uniform(1,50)
    particles = Molecule[]
    push!(particles, Molecule([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], Black_Hole_Mass, [0.0,0.0,0.0])) #super massive black hole
    for i = 2:num_particles
        θ = 2π*rand(u1)  #Uniform random numbers between [0,1]
        r = R*rand(u1)  #Uniform random numbers between [0,R]
        z = 0
        rr = [r*cos(θ), r*sin(θ), z]
        v = cross([0,0,sqrt(G*Black_Hole_Mass/(norm(rr)^3))], rr) #v=Ω×r
        mass = rand(u50) #Uniform random numbers between [1,50]
        push!(particles, Molecule(rr, v, mass, [0.0,0.0,0.0]))
            
    end
    particles
end

function main()
    t = 13     #Gaños
    dt = 0.005 #Gaños
    n = round(Int64, t/dt)
    N = 500
    particles = rand_kepler_particles2(N,4e6,30)
    display(show_particles(particles,Lim=35))
    frames = Propagate(particles, n, Δt=dt, Nframes=300, show_progress=true)
    animate_frames(frames, frame_rate=5, quadview = true, file_name = "PNanimation", Lim=40)
    savefig(graph_vel(particles),"rvplot.png")
    println("done")
end

@time main()