function VR2=pvl_WVM_compute_VR(dist,tmscales,CS)

%computes variablity reduction (VR), which is a funciton of the distances
%between sites (dist), the timescales, and the clouds speed (CS)

%% old function (uses all points in dist, was slow due to often excessive size of dist)
%  fn=@(x)abs((x.^2-x)./2-length(dist));
% [len_dist]=round(fminsearch(fn,sqrt(length(dist))));
% 
% for i=1:length(tmscales)
%     r1=exp(-dist./(CS/2*tmscales(i)));
%     rhosum=sum(r1);
%     VR(i)=1./((2*rhosum+len_dist)./(len_dist.^2));
% end

ld=10000; %sample only 1000 points for speed, minimal impact on VR (typically order 10^-3)

r=randi(length(dist),ld,1);
d1=dist(r);
 fn2=@(x)abs((x.^2-x)./2-length(d1));
[len_dist2]=round(fminsearch(fn2,sqrt(length(d1))));


for i=1:length(tmscales)
    r2=exp(-d1./(CS/2*tmscales(i)));
    rhosum2=sum(r2);
    VR2(i)=1./((2*rhosum2+len_dist2)./(len_dist2.^2));
end
