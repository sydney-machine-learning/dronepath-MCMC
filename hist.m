x = iterations;
hist3([betas(126:x), deltas(126:x)],[20,20], 'CDataMode','auto','FaceColor','interp')
% histogram2(betas, deltas,20)
% histogram2(betas(x/4:x), deltas(x/4:x),20, "DisplayStyle","bar3")
% histogram(betas(iterations/4:iterations), "Normalization","pdf")
% plot(X,betas(X),Color=[cl 1-cl cl*0.5])   
% histogram(postprob)