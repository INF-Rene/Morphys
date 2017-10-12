mainfolder='D:\Test abfbatch\Metadata';
dir_info = 'D:\Test abfbatch';
fn_info  = 'AbfInfo4Import_Rene_Test.xlsx';
tt = readtable(fullfile(dir_info,fn_info),'Sheet','Channel');
tt(:,~ismember(tt.Properties.VariableNames,{'cellName','ABF_name','ampChannel'}))=[];

abf = ImportCSVstoStruct(mainfolder, 'testbatch.mat');
ccabf = getccstep(abf);

%find cells that have multiple cc step files
tt{:,2}=cellfun(@(x) strcat(x ,'.abf'), tt{:,2}, 'UniformOutput', false);
cells=tt{ismember(tt{:,2},{ccabf.filename}),1};
[C,ia,ic] = unique(cells);
[counts,~]=hist(ic,ia);
if any(counts>1))
    %determine which abf of the duplicates should be used for analysis
end




