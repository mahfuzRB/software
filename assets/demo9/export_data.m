fclose all;
clear;

convert_Aq_Aq;

load n_bfs.dat;

N = n_bfs;
for in = 1:N
    in
    s = num2str((in-1)/1000,'%4.3f');
    fid = fopen(['bf.gmv.',s(3:5)],'r');
    while not(feof(fid))
        header = fscanf(fid,'%s',1);
        if strcmp(header,'nodes') && in == 1
            node_num = fscanf(fid,'%d',1);
            node(:,1) = fscanf(fid,'%f',node_num);
            node(:,2) = fscanf(fid,'%f',node_num);
            node(:,3) = fscanf(fid,'%f',node_num);
        elseif strcmp(header,'cells') && in == 1
            elem_num = fscanf(fid,'%d',1);            
            elem = fscanf(fid,['%*s %d %d %d %d %d'],[5,elem_num]);
            elem = elem';            
            elem = elem(:,2:5);   
        elseif strcmp(header,'variable')
            fscanf(fid,'%s',1); % u
            numf = fscanf(fid,'%d',1); % 1?
            for i = 1:numf
                Z(1:node_num,in) = fscanf(fid,'%f',node_num);
            end
        elseif strcmp(header,'material')
            fscanf(fid,'%d',2); % 1 0?
        elseif strcmp(header,'proc_0')

        end
    end
    fclose(fid);
end

save read_data;

clear;
load read_data;

mesh.p = node';
elem = [elem(:,[1,3,2]);...
        elem(:,[1,2,4]);...
        elem(:,[2,3,4]);...
        elem(:,[1,4,3])];

tmp = elem;
tmp = sort(tmp,2);
[ttmp,Is,Iv] = unique(tmp,'rows');
[sIv,I] = sort(Iv);
dsIv = diff(sIv);
I1 = all([dsIv(1:end-1),dsIv(2:end)],2);
I2 = find(I1);
elem = elem(I(I2+1),:);

% for i = 1:size(elem,1)
%     c = cross(node(elem(i,2),:)-node(elem(i,1),:),node(elem(i,3),:)-node(elem(i,1),:));
%     if dot(c,sum(node(elem(i,:),:))) < 0
%         elem(i,:) = elem(i,[1,3,2]);
%     end
% end

% patch('Vertices',node,'Faces',elem,'FaceColor','none'); axis equal;

t = elem';
t_all = unique(t(:));
node = node(t_all,:);
t_map = 1:length(t_all);
[t_all,I] = sort(t_all);
t_map = t_map(I);
tt_map = 1:max(t_all);
tt_map(t_all) = t_map;
elem = tt_map(t(1:3,:))';
t_map = tt_map;
    
mesh.p = node';
mesh.t = [elem';ones(1,size(elem,1))];
nodemap = [1];

save mod_data;

mesh.t(1:3,:) = mesh.t([1,3,2],:);
fid = fopen('geometry.dat','w');
fprintf(fid,'%d ',size(mesh.p,2));
fprintf(fid,'%3.2f ',mesh.p(:));
fprintf(fid,'%d ',size(nodemap,1));
fprintf(fid,'%d ',size(mesh.t,2));
tt = mesh.t(1:3,:);
fprintf(fid,'%d ',tt(:)-1);
node_reg = zeros(size(mesh.p,2),1);
for i = 1:size(nodemap,1)
    I = find(mesh.t(4,:) == i);
    I = mesh.t(1:3,I); I = unique(I(:));
    node_reg(I) = i;
end
fprintf(fid,'%d ',node_reg-1);
fprintf(fid,'%d ',mesh.t(4,:)-1);
fclose(fid);
% keyboard
% fid = fopen('Z_000.dat','w');
% fprintf(fid,'%15.14e ',Z(:));
% fclose(fid);

fid = fopen('calN.dat','w');
fprintf(fid,'%d',size(mesh.p,2));
fclose(fid);

for i = 1:size(Z,2)
    s = num2str((i-1)/1000,'%4.3f');
    fid = fopen(['Z_000_',s(3:5),'.bin'],'w');
    fwrite(fid,Z(t_all,i),'single');
    fclose(fid);
end

system('javac AffineFunctions.java');
CD = cd;
copyfile('AffineFunctions.class','d:\Eclipse\Windows\android\platforms\android-2.0\tools\AffineFunctions.class');
cd('d:\Eclipse\Windows\android\platforms\android-2.0\tools');
system('dx --dex --output=AffineFunctions.jar AffineFunctions.class');
delete('AffineFunctions.class');
cd(CD);
copyfile('d:\Eclipse\Windows\android\platforms\android-2.0\tools\AffineFunctions.jar','AffineFunctions.jar');
delete('d:\Eclipse\Windows\android\platforms\android-2.0\tools\AffineFunctions.jar');
