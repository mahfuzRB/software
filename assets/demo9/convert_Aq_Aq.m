clear;

fid = fopen('Aq_Aq_norms.dat','r');
data = fscanf(fid,'%f');
fclose(fid);

fid = fopen('n_bfs.dat','r');
N = fscanf(fid,'%d');
fclose(fid);

count = 0;
Qa = (-1 + sqrt(1+4*2*size(data,1)/(N^2)))/2;
for i = 1:Qa
    for j = i:Qa
        si = num2str((i-1)/1000,'%4.3f');
        sj = num2str((j-1)/1000,'%4.3f');
        fid = fopen(['Aq_Aq_',si(3:5),'_',sj(3:5),'_norms.bin'],'w');
        fwrite(fid,data(count*(N^2)+1:(count+1)*(N^2)),'double');
        fclose(fid);
        count = count + 1;
        [i,j]
    end
end

fid = fopen('Aq_Aq_norms.bin','w');
fwrite(fid,data,'double');
fclose(fid);