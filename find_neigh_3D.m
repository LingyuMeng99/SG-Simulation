function [n1, n2, n3, n4, n5, n6]=find_neigh_3D(n,L1,L2,L3)
n1 = n - 1;
n2 = n + 1;

n3 = n - L1;
n4 = n + L1;

n5 = n - L1*L2;
n6 = n + L1*L2;

n_x = mod(n-1,L1)+1;
n_y = mod(floor((n-1)/L1),L2)+1;
n_z = mod(floor((n-1)/L1/L2),L3)+1;

if n_x == 1
    n1 = n1 + L1;
elseif n_x == L1
    n2 = n2 - L1;
end

if n_y == 1
    n3 = n3 + L1*L2;
elseif n_y == L2
    n4 = n4 - L1*L2;
end

if n_z == 1
    n5 = n5 + L1*L2*L3;
elseif n_z == L3
    n6 = n6 - L1*L2*L3;
end

end