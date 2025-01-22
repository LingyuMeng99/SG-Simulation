%% We find the all the particles that are connected to the focused one
function [cluster,cluster_num] = find_cluster_3D(n,state,L1,L2,L3)

N = L1*L2*L3;

cluster = zeros(N,1);
cluster_num = 0;
cluster(cluster_num+1) = n;
cluster_num = cluster_num + 1;
unexplore = zeros(N,1);

newmember = zeros(N,1);
newmember_num = 0;
newmember(newmember_num+1) = n;
newmember_num = newmember_num + 1;


while(1)
    
    count_un = 0;
    
    for i = 1:newmember_num
        n_temp = newmember(i);
        [n1,n2,n3,n4,n5,n6] = find_neigh_3D(n_temp,L1,L2,L3);
        unexplore(count_un+1:count_un+6) = [n1 n2 n3 n4 n5 n6];
        count_un = count_un + 6;
    end
    
    %  unexplore_unique=unique(unexplore(1:count_un));
    %  count_unique=length(unexplore_unique);
    
    
    newmember_num = 0;
    for k = 1:count_un
        
        %ismember(unexplore(k),cluster)==0
        
        if  state(unexplore(k)) == 1 %&& ismember(unexplore(k),cluster)==0
            
            % check if unexplore(k) is in the cluster
            check = 0;
            m = 1;
            while(1)
                if unexplore(k) == cluster(m)
                    check = 1;
                    break;
                end
                
                if m == cluster_num
                    break;
                end
                
                m = m+1;
            end
            
            if check == 0
                
                cluster(cluster_num+1) = unexplore(k);
                newmember(newmember_num+1) = unexplore(k);
                cluster_num = cluster_num + 1;
                newmember_num = newmember_num + 1;
            end
        end
    end
    
    if newmember_num == 0
        break;
    end
    
    
end

cluster = cluster(1:cluster_num);
end




