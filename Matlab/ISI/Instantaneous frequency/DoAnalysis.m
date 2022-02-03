% reads json files for every cell and saves the talbe with features 
%(including instantaneous freq) 
% each cell

for i=1:size(cellnames,1)
cell_id=cellnames{i};
file2analyse=fullfile(dir2analyse,cell_id,'output.json');
readJson(file2analyse, cell_id, path2saveSpikes);
end