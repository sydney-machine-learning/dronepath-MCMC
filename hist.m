x = 100:2000;
% hist3([betas(x), deltas(x)],[20,20], 'CDataMode','auto','FaceColor','interp')
% hist3([betas(x), deltas(x)],[20,20], 'CDataMode','auto','FaceColor','interp') %histogram
% histogram2(betas, deltas,20)
% histogram2(betas(x/4:x), deltas(x/4:x),20, "DisplayStyle","bar3")
% histogram(betas(iterations/4:iterations), "Normalization","pdf")
plot3(x,betas(x), deltas(x))   %traceplots
% plot(x, ll(x)) %loglikelihood plot
% histogram(postprob)