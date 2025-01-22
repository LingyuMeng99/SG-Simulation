function [pos,Nd]=find_droplet_3D(phi,phi_c,L1,L2,L3)

state = phi > phi_c;
statesq = reshape(state,L1,L2,L3);
ball = find(state(:)==1);
Nall = length(ball);
check = zeros(Nall,1);
pos = {};
Nd = 0;
for i = 1:Nall
    if check(i) == 1
        continue;
    else
        nb = ball(i);
        [cluster,cnum] = find_cluster_3D(nb,statesq,L1,L2,L3);
        check = check + ismember(ball,cluster);
        if cnum > 1
            
            pos{end+1} = cluster;
            Nd = Nd+1;            
        end
    end    
end