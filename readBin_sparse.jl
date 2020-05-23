using Plots;

#plotly(size=(500,500))
gr(size=(500,500))

filename = "data_hasRe.bin";
stream = open(filename, "r");

total_step = 5000
Ptotal = 500*500
Gvals = Array{Int32}(undef, Ptotal);
println("new session")
# Gvals = Gvals .* 0;

sum_step = 1200;

for i=1:total_step

 Ntotal=read(stream,Int32);
 println("Frame: ", i,". Number of rearrangers: ", Ntotal)

  if Ntotal>0
   Gvals_frame = Array{Int32}(undef, Ntotal);
   read!(stream, Gvals_frame);
   Gvals_frame = Gvals_frame .+ 1;
   global Gvals[Gvals_frame] = fill(1,Ntotal);
  end

  if (i%sum_step==0)
  global Gval_mat = reshape(Gvals,(500,500));
  handle = heatmap(Gval_mat,legend = :none)
  fn = "./frame_$(i).png";
  png(handle,fn)
  global Gvals = Gvals .* 0;
  end

end

close(stream)
