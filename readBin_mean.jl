using Plots;
# using PyPlot;

# plotly(size=(500,500))
gr(size=(600,300))
# pyplot(size=(500,500))

filename = "data_mean.bin";
stream = open(filename, "r");

total_step = 4760
softness_collect = Array{Float64}(undef, total_step);
strain_collect = Array{Float64}(undef, total_step);
println("new session")

for i=1:total_step
Nstep=read(stream,Float64);
softness_collect[i] = read(stream, Float64);
strain_collect[i] = sqrt(read(stream, Float64));
println("Step: ", Nstep, "  MeanSoftness: ", softness_collect[i],"  Total energy:", strain_collect[i])
end
close(stream)


ps=plot(softness_collect,legend=false)
pe=plot(strain_collect,legend=false)
hp=plot(ps,pe,layout=2)
fn = "./softness_strain.png";
png(hp,fn)
