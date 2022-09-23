range = 120:2000;
% hist3([betas(range), deltas(range)],[20,20], 'CDataMode','auto','FaceColor','interp') %histogram
% histogram2(betas, deltas,20)
% histogram2(betas(x/4:x), deltas(x/4:x),20, "DisplayStyle","bar3")
plot3(range,betas(range), deltas(range))   %traceplots
% plot(range, ll(range)) %loglikelihood plot
% plot3(betas,deltas,ll(2:501))