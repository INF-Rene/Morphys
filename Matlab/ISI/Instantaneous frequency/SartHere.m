function [cellnames,dir2analyse]=StartHere()
dir2analyse='/Users/annagalakhova/PhD INF CNCR VU/DATA/Instantaneous frequency/ipfx_output_json';
filelist=dir(dir2analyse);
k=1;
for i=1:size(filelist,1)
    if  ~strcmpi(filelist(i).name,'.')&& ~strcmpi(filelist(i).name,'.DS_Store') && ~strcmpi(filelist(i).name,'..')
    cellnames{k,1}=filelist(i).name;
    k=k+1;
    end
end
end