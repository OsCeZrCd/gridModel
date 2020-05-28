using Plots;
# using PyPlot;

# plotly(size=(500,500))
gr(size=(800,800))
# pyplot(size=(500,500))

filename = "data_Softness.bin";
stream = open(filename, "r");

total_step = 200
Ntotal = 400*400
Gvals = Array{Float64}(undef, Ntotal);
println("new session")

for i=1:total_step

Ntotal=read(stream,Float64);
read!(stream, Gvals);
# pos_G = findall(x -> x > 0, Gvals);
# L_pos = length(pos_G);
println("Step: ", (i-1) * 50 +1, "  Total: ", Ntotal)

Gval_mat = reshape(Gvals,(400,400));

handle = heatmap(Gval_mat,clim=(-3,3))
# plot!(clims(-0.25,0.25))
fn = "./softness_$((i-1) * 20 +1).png";
png(handle,fn)

end
close(stream)
