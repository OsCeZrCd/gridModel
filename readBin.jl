using Plots;

#plotly(size=(500,500))
gr(size=(500,500))

filename = "data_hasRe.bin";
stream = open(filename, "r");

total_step = 64
Ntotal = 250000
Gvals = Array{Int32}(undef, Ntotal);
println("new session")

for i=1:total_step

Ntotal=read(stream,Int32);
read!(stream, Gvals);
pos_G = findall(x -> x > 0, Gvals);
L_pos = length(pos_G);
println("Total: ", Ntotal,". Number of rearrangers: ", L_pos)

Gval_mat = reshape(Gvals,(500,500));


handle = heatmap(Gval_mat,legend = :none)
fn = "./frame_$(i).png";
png(handle,fn)

end
close(stream)
