%% Make video
v = VideoWriter('output','MPEG-4');
open(v);
for i = 1:N_rec
    cla;
    ti = i*dt_rec;
    phi_shift_vec = [6,16,8];
    phi_plot = circshift(reshape(phi_RNA_rec(:,i),L,L,L),phi_shift_vec);
    phi_RNA_ref = circshift(reshape(phi_RNA_rec(:,i),L,L,L),phi_shift_vec);

    % phi_idx_1 = 12:26;%7:21;
    % phi_idx_2 = 9:23;%11:25;
    % phi_idx_3 = 16:30;
    % phi_plot = phi_plot(phi_idx_1,phi_idx_2,phi_idx_3);
    % phi_RNA_ref = phi_RNA_ref(phi_idx_1,phi_idx_2,phi_idx_3);

    set(gcf,'position',[561 527 520 470])
    % Plot condensates
    N_isovalue = 5;
    phi_plot_max = 15;
    % color_base = lines(1); color_base = color_base(1,:);
    color_base = [197,149,220]./255;
    color_list = [linspace(1,color_base(1),2*N_isovalue+1)',linspace(1,color_base(2),2*N_isovalue+1)',linspace(1,color_base(3),2*N_isovalue+1)'];
    for isovalue_i = 1:N_isovalue
        isovalue = isovalue_i/N_isovalue*phi_plot_max;
        surf = isosurface(phi_plot,isovalue);
        p = patch(surf);
        isonormals(phi_plot,p);
        set(p,'FaceColor',color_list(isovalue_i+N_isovalue+1,:),'EdgeColor','none','FaceAlpha',isovalue_i/N_isovalue/2);
        daspect([1,1,1])
    end
    % Plot background
    color_list_plot = [linspace(1,color_base(1),200)',linspace(1,color_base(2),200)',linspace(1,color_base(3),200)'];

    phi_plot_bg = mean(phi_plot(phi_RNA_ref(:)<1.5));%disp(phi_plot_bg);
    vert = [0 0 0; L 0 0; L L 0; 0 L 0; 0 0 L; L 0 L; L L L; 0 L L];
    fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
    patch('Vertices',vert,'Faces',fac,...
        'FaceVertexCData',color_list_plot(100+round(phi_plot_bg/phi_plot_max*100),:),'FaceColor','flat','FaceAlpha',0.1)
    
    view(3); axis square
    view([1 0.2 0.2]);
    colorbar; colormap(color_list_plot(100:end,:)); clim([0 phi_plot_max])
    xlim([0 L]);ylim([0 L]);zlim([0 L]);
    set(gca,'Box','On')
    set(gca,'BoxStyle','Full')
    set(gca,'position',[0.1190,0.0800,0.7094,0.8106])
    set(gca,'FontSize',15)
    
    if ti == t_stress_on
        title(['t = ',num2str(i*dt_rec),newline,'Concentration of RNA, Stress On'],...
            'Fontsize',20)
        for count=1:60
        frame = getframe(gcf);
        tmp = frame2im(frame);
        writeVideo(v,tmp);
        end
    elseif ti==t_stress_on+t_stress_duration
        title(['t = ',num2str(i*dt_rec),newline,'Concentration of RNA, Stress Off'],...
            'Fontsize',20)
        for count=1:60
        frame = getframe(gcf);
        tmp = frame2im(frame);
        writeVideo(v,tmp);
        end
    elseif mod(i,2)==0
        title(['t = ',num2str(i*dt_rec),newline,'Concentration of RNA'],'Fontsize',20)
        frame = getframe(gcf);
        tmp = frame2im(frame); 
        writeVideo(v,tmp);
    end    
end
close(v)





%% 3D plot

    i = 200;
    % ti = i*dt_rec;
    phi_shift_vec = [2,6,-7];
    phi_plot = circshift(reshape(phi_RNA_rec(:,i),L,L,L),phi_shift_vec);
    phi_RNA_ref = circshift(reshape(phi_RNA_rec(:,i),L,L,L),phi_shift_vec);

    length_unit = 0.55;

    % phi_idx_1 = 12:26;%7:21;
    % phi_idx_2 = 9:23;%11:25;
    % phi_idx_3 = 16:30;
    % phi_plot = phi_plot(phi_idx_1,phi_idx_2,phi_idx_3);
    % phi_RNA_ref = phi_RNA_ref(phi_idx_1,phi_idx_2,phi_idx_3);

    set(gcf,'position',[561 527 520 470])
    % Plot condensates
    N_isovalue = 5;
    phi_plot_max = 15;
    color_base = lines(1); color_base = color_base(1,:);
    color_list = [linspace(1,color_base(1),2*N_isovalue+1)',linspace(1,color_base(2),2*N_isovalue+1)',linspace(1,color_base(3),2*N_isovalue+1)'];
    for isovalue_i = 1:N_isovalue
        isovalue = isovalue_i/N_isovalue*phi_plot_max;
        [x,y,z] = meshgrid((1:L).*length_unit);
        surf = isosurface(x,y,z,phi_plot,isovalue);
        p = patch(surf);
        isonormals(phi_plot,p);
        set(p,'FaceColor',color_list(isovalue_i+N_isovalue+1,:),'EdgeColor','none','FaceAlpha',isovalue_i/N_isovalue/2);
        daspect([1,1,1])
    end
    % Plot background
    color_list_plot = [linspace(1,color_base(1),200)',linspace(1,color_base(2),200)',linspace(1,color_base(3),200)'];

    phi_plot_bg = mean(phi_plot(phi_RNA_ref(:)<1.5));disp(phi_plot_bg);
    vert = [0 0 0; L 0 0; L L 0; 0 L 0; 0 0 L; L 0 L; L L L; 0 L L].*length_unit;
    fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
    patch('Vertices',vert,'Faces',fac,...
        'FaceVertexCData',color_list_plot(100+round(phi_plot_bg/phi_plot_max*100),:),'FaceColor','flat','FaceAlpha',0.1)
    
    view(3); axis square
    view([1 -0.8 0.5]);
    colorbar; colormap(color_list_plot(100:end,:)); clim([0 phi_plot_max])
    xlim([0 L*length_unit]);ylim([0 L*length_unit]);zlim([0 L*length_unit]);
    set(gca,'Box','On')
    set(gca,'BoxStyle','Full')
    set(gca,'position',[0.1190,0.0800,0.7094,0.8106])
    set(gca,'FontSize',15)
    % set(gca,'Projection','perspective');


%% 3D data 3D plot PART slice balls
phi_80splot_rec = phi_80s_rec+phi_80s1RNA_rec+2*phi_80s2RNA_rec+3*phi_80s3RNA_rec+4*phi_80s4RNA_rec+5*phi_80s5RNA_rec;

i = 200;
ti = i*dt_rec;
phi_idx_1 = (12:13)+2;
phi_idx_2 = (2:3)+8;
phi_idx_3 = (14:15)+16;
phi_shift_vec = [2,6,-7];
phi_40s_plot = circshift(reshape(phi_40s_rec(:,i),L,L,L),phi_shift_vec);
phi_40s_plot = phi_40s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_43s_plot = circshift(reshape(phi_43s_rec(:,i),L,L,L),phi_shift_vec);
phi_43s_plot = phi_43s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_60s_plot = circshift(reshape(phi_60s_rec(:,i),L,L,L),phi_shift_vec);
phi_60s_plot = phi_60s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_80s_plot = circshift(reshape(phi_80splot_rec(:,i),L,L,L),phi_shift_vec);
phi_80s_plot = phi_80s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_RNA_ref = circshift(reshape(phi_RNA_rec(:,i),L,L,L),phi_shift_vec);
phi_RNA_ref = phi_RNA_ref(phi_idx_1,phi_idx_2,phi_idx_3);
phi_RNA_plot = phi_RNA_ref;

color_list = [133 173 214; 255 158 158; 82 156 111; 255 214 92]./255;%[169 197 226; 255 184 184; 145 222 176; 255 239 189]./255;
size_list = [10,10,12,14,8];

% color_base = [0.7,0.7,0.7];
% color_list_plot = [linspace(1,color_base(1),200)',linspace(1,color_base(2),200)',linspace(1,color_base(3),200)'];
% imagesc(phi_RNA_ref,'Interpolation','bilinear'); colorbar; clim([0 30]); colormap(color_list_plot);
% set(gca,'YDir','normal')
% hold on

[L1,L2,L3] = size(phi_RNA_ref);
[pos_SG,N_SG] = find_droplet_3D(phi_RNA_ref,5,L1,L2,L3);

ball_num_mat = zeros(4,N_SG+1);
% inside
for SG_i = 1:N_SG
    ball_num_mat(1,SG_i) = sum(phi_40s_plot(pos_SG{SG_i}));
    ball_num_mat(2,SG_i) = sum(phi_43s_plot(pos_SG{SG_i}));
    ball_num_mat(3,SG_i) = sum(phi_60s_plot(pos_SG{SG_i}));
    ball_num_mat(4,SG_i) = sum(phi_80s_plot(pos_SG{SG_i}));
    % ball_num_mat(5,SG_i) = sum(phi_RNA_plot(pos_SG{SG_i}));
end
% outside
pos_Cyto = 1:length(phi_40s_plot(:));
for SG_i = 1:N_SG
    pos_Cyto = setdiff(pos_Cyto,pos_SG{SG_i});
end
ball_num_mat(1,N_SG+1) = sum(phi_40s_plot(pos_Cyto));
ball_num_mat(2,N_SG+1) = sum(phi_43s_plot(pos_Cyto));
ball_num_mat(3,N_SG+1) = sum(phi_60s_plot(pos_Cyto));
ball_num_mat(4,N_SG+1) = sum(phi_80s_plot(pos_Cyto));
% ball_num_mat(5,N_SG+1) = sum(phi_RNA_plot(pos_Cyto));

ball_num_mat = round(ball_num_mat.*100) % .*100 for real number


plot_interval = 0.002;
plot_interval_x = plot_interval; plot_interval_y = plot_interval; plot_interval_z = plot_interval;

pos_40s_plot = [];
pos_43s_plot = [];
pos_60s_plot = [];
pos_80s_plot = [];

% inside plot
for SG_i = 1:N_SG
    pos_list = pos_SG{SG_i};
    pos_y_list = mod(pos_list-1,L1)+1;
    pos_x_list = mod(floor((pos_list-1)/L1),L2)+1;
    pos_z_list = mod(floor((pos_list-1)/L1/L2),L3)+1;
    
    add_length = round((1/plot_interval_x)*(1/plot_interval_y)*(1/plot_interval_z));
    pos_x_mat = zeros(1,round(length(pos_list)*add_length));
    pos_y_mat = zeros(1,round(length(pos_list)*add_length));
    pos_z_mat = zeros(1,round(length(pos_list)*add_length));
    
    add_idx_first = 1;
    add_idx_next = add_length;
    for pos_list_i = 1:length(pos_list)
        pos_x_add = pos_x_list(pos_list_i)-1:plot_interval_x:pos_x_list(pos_list_i)+0-plot_interval_x;
        pos_y_add = pos_y_list(pos_list_i)-1:plot_interval_y:pos_y_list(pos_list_i)+0-plot_interval_y;
        pos_z_add = (pos_z_list(pos_list_i)-1:plot_interval_z:pos_z_list(pos_list_i)+0-plot_interval_z);
        [pos_x_add_mat,pos_y_add_mat,pos_z_add_mat] = meshgrid(pos_x_add,pos_y_add,pos_z_add);
        
        pos_x_mat(add_idx_first:add_idx_next) = pos_x_add_mat(:)';
        pos_y_mat(add_idx_first:add_idx_next) = pos_y_add_mat(:)';
        pos_z_mat(add_idx_first:add_idx_next) = pos_z_add_mat(:)';

        add_idx_first = add_idx_first + add_length;
        add_idx_next = add_idx_next + add_length;
    end
    plot_idx_rand = randperm(length(pos_x_mat(:)),sum(ball_num_mat(:,SG_i)));
    pos_x_plot = pos_x_mat(plot_idx_rand).*0.55;
    pos_y_plot = pos_y_mat(plot_idx_rand).*0.55;
    pos_z_plot = pos_z_mat(plot_idx_rand).*0.55;

    for plot_i = 1:4
        ball_num_plot = sum(ball_num_mat(1:plot_i-1,SG_i))+1:sum(ball_num_mat(1:plot_i,SG_i));
        plot3(pos_x_plot(ball_num_plot),pos_y_plot(ball_num_plot),pos_z_plot(ball_num_plot),'.','Markersize',size_list(plot_i), ...
            'Color',color_list(plot_i,:),'Linewidth',0.01);
        view(3)
        hold on

        if plot_i == 1
            pos_40s_plot = [pos_40s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
        elseif plot_i == 2
            pos_43s_plot = [pos_43s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
        elseif plot_i == 3
            pos_60s_plot = [pos_60s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
        elseif plot_i == 4
            pos_80s_plot = [pos_80s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
        end
    end

end

% outside plot
pos_list = pos_Cyto;
pos_y_list = mod(pos_list-1,L1)+1;
pos_x_list = mod(floor((pos_list-1)/L1),L2)+1;
pos_z_list = mod(floor((pos_list-1)/L1/L2),L3)+1;

add_length = (1/plot_interval_x)*(1/plot_interval_y)*(1/plot_interval_z);
pos_x_mat = zeros(1,round(length(pos_list)*add_length));
pos_y_mat = zeros(1,round(length(pos_list)*add_length));
pos_z_mat = zeros(1,round(length(pos_list)*add_length));

add_idx_first = 1;
add_idx_next = add_length;
for pos_list_i = 1:length(pos_list)
    pos_x_add = pos_x_list(pos_list_i)-1:plot_interval_x:pos_x_list(pos_list_i)+0-plot_interval_x;
    pos_y_add = pos_y_list(pos_list_i)-1:plot_interval_y:pos_y_list(pos_list_i)+0-plot_interval_y;
    pos_z_add = (pos_z_list(pos_list_i)-1:plot_interval_z:pos_z_list(pos_list_i)+0-plot_interval_z);
    [pos_x_add_mat,pos_y_add_mat,pos_z_add_mat] = meshgrid(pos_x_add,pos_y_add,pos_z_add);
    
    pos_x_mat(add_idx_first:add_idx_next) = pos_x_add_mat(:)';
    pos_y_mat(add_idx_first:add_idx_next) = pos_y_add_mat(:)';
    pos_z_mat(add_idx_first:add_idx_next) = pos_z_add_mat(:)';

    add_idx_first = add_idx_first + add_length;
    add_idx_next = add_idx_next + add_length;
end

plot_idx_rand = randperm(length(pos_x_mat(:)),sum(ball_num_mat(:,N_SG+1)));
pos_x_plot = pos_x_mat(plot_idx_rand).*0.55;
pos_y_plot = pos_y_mat(plot_idx_rand).*0.55;
pos_z_plot = pos_z_mat(plot_idx_rand).*0.55;

for plot_i = 1:4
    ball_num_plot = sum(ball_num_mat(1:plot_i-1,N_SG+1))+1:sum(ball_num_mat(1:plot_i,N_SG+1));
    plot3(pos_x_plot(ball_num_plot),pos_y_plot(ball_num_plot),pos_z_plot(ball_num_plot),'.','Markersize',size_list(plot_i), ...
        'Color',color_list(plot_i,:),'Linewidth',0.01);
    view(3)
    hold on

    if plot_i == 1
        pos_40s_plot = [pos_40s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
    elseif plot_i == 2
        pos_43s_plot = [pos_43s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
    elseif plot_i == 3
        pos_60s_plot = [pos_60s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
    elseif plot_i == 4
        pos_80s_plot = [pos_80s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
    end
end

% contourf(phi_RNA_ref>=1.5,1,'FaceColor','none','Linewidth',1.5)
axis('square')
% set(gcf,'position',[561 527 520 470])
% set(gca,'position',[0.1190,0.0800,0.7094,0.8106])
% set(gca,'position',[0.10,0.0800,0.8094,0.8106])
set(gca,'FontSize',15)




%% 3D data 3D plot All slice balls
phi_80splot_rec = phi_80s_rec+phi_80s1RNA_rec+2*phi_80s2RNA_rec+3*phi_80s3RNA_rec+4*phi_80s4RNA_rec+5*phi_80s5RNA_rec;

i = 200;
ti = i*dt_rec;
phi_idx_1 = 2:16;%12:26;%7:21;
phi_idx_2 = 8:22;%9:23;%11:25;
phi_idx_3 = 16:30;%16:30;
phi_shift_vec = [2,6,-7];%[-1,3,-1];
phi_40s_plot = circshift(reshape(phi_40s_rec(:,i),L,L,L),phi_shift_vec);
phi_40s_plot = phi_40s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_43s_plot = circshift(reshape(phi_43s_rec(:,i),L,L,L),phi_shift_vec);
phi_43s_plot = phi_43s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_60s_plot = circshift(reshape(phi_60s_rec(:,i),L,L,L),phi_shift_vec);
phi_60s_plot = phi_60s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_80s_plot = circshift(reshape(phi_80splot_rec(:,i),L,L,L),phi_shift_vec);
phi_80s_plot = phi_80s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_RNA_ref = circshift(reshape(phi_RNA_rec(:,i),L,L,L),phi_shift_vec);
phi_RNA_ref = phi_RNA_ref(phi_idx_1,phi_idx_2,phi_idx_3);
phi_RNA_plot = phi_RNA_ref;

color_list = [133 173 214; 255 158 158; 82 156 111; 255 214 92]./255;%[169 197 226; 255 184 184; 145 222 176; 255 239 189]./255;
size_list = [10,10,12,14,8];

% color_base = [0.7,0.7,0.7];
% color_list_plot = [linspace(1,color_base(1),200)',linspace(1,color_base(2),200)',linspace(1,color_base(3),200)'];
% imagesc(phi_RNA_ref,'Interpolation','bilinear'); colorbar; clim([0 30]); colormap(color_list_plot);
% set(gca,'YDir','normal')
% hold on

[L1,L2,L3] = size(phi_RNA_ref);
[pos_SG,N_SG] = find_droplet_3D(phi_RNA_ref,5,L1,L2,L3);

ball_num_mat = zeros(4,N_SG+1);
% inside
for SG_i = 1:N_SG
    ball_num_mat(1,SG_i) = sum(phi_40s_plot(pos_SG{SG_i}));
    ball_num_mat(2,SG_i) = sum(phi_43s_plot(pos_SG{SG_i}));
    ball_num_mat(3,SG_i) = sum(phi_60s_plot(pos_SG{SG_i}));
    ball_num_mat(4,SG_i) = sum(phi_80s_plot(pos_SG{SG_i}));
    % ball_num_mat(5,SG_i) = sum(phi_RNA_plot(pos_SG{SG_i}));
end
% outside
pos_Cyto = 1:length(phi_40s_plot(:));
for SG_i = 1:N_SG
    pos_Cyto = setdiff(pos_Cyto,pos_SG{SG_i});
end
ball_num_mat(1,N_SG+1) = sum(phi_40s_plot(pos_Cyto));
ball_num_mat(2,N_SG+1) = sum(phi_43s_plot(pos_Cyto));
ball_num_mat(3,N_SG+1) = sum(phi_60s_plot(pos_Cyto));
ball_num_mat(4,N_SG+1) = sum(phi_80s_plot(pos_Cyto));
% ball_num_mat(5,N_SG+1) = sum(phi_RNA_plot(pos_Cyto));

ball_num_mat = round(ball_num_mat.*100) % .*100 for real number


plot_interval = 0.05;
plot_interval_x = plot_interval; plot_interval_y = plot_interval; plot_interval_z = plot_interval;

pos_40s_plot = [];
pos_43s_plot = [];
pos_60s_plot = [];
pos_80s_plot = [];

% inside plot
for SG_i = 1:N_SG
    pos_list = pos_SG{SG_i};
    pos_y_list = mod(pos_list-1,L1)+1;
    pos_x_list = mod(floor((pos_list-1)/L1),L2)+1;
    pos_z_list = mod(floor((pos_list-1)/L1/L2),L3)+1;
    
    add_length = round((1/plot_interval_x)*(1/plot_interval_y)*(1/plot_interval_z));
    pos_x_mat = zeros(1,round(length(pos_list)*add_length));
    pos_y_mat = zeros(1,round(length(pos_list)*add_length));
    pos_z_mat = zeros(1,round(length(pos_list)*add_length));
    
    add_idx_first = 1;
    add_idx_next = add_length;
    for pos_list_i = 1:length(pos_list)
        pos_x_add = pos_x_list(pos_list_i)-1:plot_interval_x:pos_x_list(pos_list_i)+0-plot_interval_x;
        pos_y_add = pos_y_list(pos_list_i)-1:plot_interval_y:pos_y_list(pos_list_i)+0-plot_interval_y;
        pos_z_add = (pos_z_list(pos_list_i)-1:plot_interval_z:pos_z_list(pos_list_i)+0-plot_interval_z);
        [pos_x_add_mat,pos_y_add_mat,pos_z_add_mat] = meshgrid(pos_x_add,pos_y_add,pos_z_add);
        
        pos_x_mat(add_idx_first:add_idx_next) = pos_x_add_mat(:)';
        pos_y_mat(add_idx_first:add_idx_next) = pos_y_add_mat(:)';
        pos_z_mat(add_idx_first:add_idx_next) = pos_z_add_mat(:)';

        add_idx_first = add_idx_first + add_length;
        add_idx_next = add_idx_next + add_length;
    end
    plot_idx_rand = randperm(length(pos_x_mat(:)),sum(ball_num_mat(:,SG_i)));
    pos_x_plot = pos_x_mat(plot_idx_rand).*0.55;
    pos_y_plot = pos_y_mat(plot_idx_rand).*0.55;
    pos_z_plot = pos_z_mat(plot_idx_rand).*0.55;

    for plot_i = 1:4
        ball_num_plot = sum(ball_num_mat(1:plot_i-1,SG_i))+1:sum(ball_num_mat(1:plot_i,SG_i));
        plot3(pos_x_plot(ball_num_plot),pos_y_plot(ball_num_plot),pos_z_plot(ball_num_plot),'.','Markersize',size_list(plot_i), ...
            'Color',color_list(plot_i,:),'Linewidth',0.01);
        view(3)
        hold on

        if plot_i == 1
            pos_40s_plot = [pos_40s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
        elseif plot_i == 2
            pos_43s_plot = [pos_43s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
        elseif plot_i == 3
            pos_60s_plot = [pos_60s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
        elseif plot_i == 4
            pos_80s_plot = [pos_80s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
        end
    end

end

% outside plot
pos_list = pos_Cyto;
pos_y_list = mod(pos_list-1,L1)+1;
pos_x_list = mod(floor((pos_list-1)/L1),L2)+1;
pos_z_list = mod(floor((pos_list-1)/L1/L2),L3)+1;

add_length = (1/plot_interval_x)*(1/plot_interval_y)*(1/plot_interval_z);
pos_x_mat = zeros(1,round(length(pos_list)*add_length));
pos_y_mat = zeros(1,round(length(pos_list)*add_length));
pos_z_mat = zeros(1,round(length(pos_list)*add_length));

add_idx_first = 1;
add_idx_next = add_length;
for pos_list_i = 1:length(pos_list)
    pos_x_add = pos_x_list(pos_list_i)-1:plot_interval_x:pos_x_list(pos_list_i)+0-plot_interval_x;
    pos_y_add = pos_y_list(pos_list_i)-1:plot_interval_y:pos_y_list(pos_list_i)+0-plot_interval_y;
    pos_z_add = (pos_z_list(pos_list_i)-1:plot_interval_z:pos_z_list(pos_list_i)+0-plot_interval_z);
    [pos_x_add_mat,pos_y_add_mat,pos_z_add_mat] = meshgrid(pos_x_add,pos_y_add,pos_z_add);
    
    pos_x_mat(add_idx_first:add_idx_next) = pos_x_add_mat(:)';
    pos_y_mat(add_idx_first:add_idx_next) = pos_y_add_mat(:)';
    pos_z_mat(add_idx_first:add_idx_next) = pos_z_add_mat(:)';

    add_idx_first = add_idx_first + add_length;
    add_idx_next = add_idx_next + add_length;
end

plot_idx_rand = randperm(length(pos_x_mat(:)),sum(ball_num_mat(:,N_SG+1)));
pos_x_plot = pos_x_mat(plot_idx_rand).*0.55;
pos_y_plot = pos_y_mat(plot_idx_rand).*0.55;
pos_z_plot = pos_z_mat(plot_idx_rand).*0.55;

for plot_i = 1:4
    ball_num_plot = sum(ball_num_mat(1:plot_i-1,N_SG+1))+1:sum(ball_num_mat(1:plot_i,N_SG+1));
    plot3(pos_x_plot(ball_num_plot),pos_y_plot(ball_num_plot),pos_z_plot(ball_num_plot),'.','Markersize',size_list(plot_i), ...
        'Color',color_list(plot_i,:),'Linewidth',0.01);
    view(3)
    hold on

    if plot_i == 1
        pos_40s_plot = [pos_40s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
    elseif plot_i == 2
        pos_43s_plot = [pos_43s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
    elseif plot_i == 3
        pos_60s_plot = [pos_60s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
    elseif plot_i == 4
        pos_80s_plot = [pos_80s_plot,[pos_x_plot(ball_num_plot);pos_y_plot(ball_num_plot);pos_z_plot(ball_num_plot)]];
    end
end

% contourf(phi_RNA_ref>=1.5,1,'FaceColor','none','Linewidth',1.5)
axis('square')
% set(gcf,'position',[561 527 520 470])
% set(gca,'position',[0.1190,0.0800,0.7094,0.8106])
% set(gca,'position',[0.10,0.0800,0.8094,0.8106])
set(gca,'FontSize',15)