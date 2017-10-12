mainfolder='D:\Test abfbatch\Metadata';
dir_info = 'D:\Test abfbatch';
fn_info  = 'AbfInfo4Import_Rene_Test.xlsx';
tt = readtable(fullfile(dir_info,fn_info),'Sheet','Channel');
tt(:,~ismember(tt.Properties.VariableNames,{'cellName','ABF_name','ampChannel'}))=[];

abf = ImportCSVstoStruct(mainfolder, 'testbatch.mat');
ccabf = getccstep(abf);

%find cells that have multiple cc step files
tt{:,2}=cellfun(@(x) strcat(x ,'.abf'), tt{:,2}, 'UniformOutput', false);
cells=tt(ismember(tt{:,2},{ccabf.filename}),:);
[C,ia,ic] = unique(cells(:,1), 'stable');
[counts,~]=hist(ic,ia);
if any(counts>1))
    %determine which abf of the duplicates should be used for analysis
    dupes={cells{ia(counts>1),1}};
    for i=1:numel(dupes)
        a=cells(ismember(cells{:,1},dupes{i}),2:3);
        idx=find(ismember({ccabf.filename}, a{:,1}));
        b={ccabf(idx).filetimestart};
        b=cellfun(@(x) strsplit(x, {':', '-', ' '}), b, 'UniformOutput', false);
        clear c
        for j=1:numel(b)
        c(j)=datetime(cellfun(@str2num, b{j}));
        end
        [~,d]=max(c);           %selects the abf that was recorded the latest
        idx(d)=[];
        ccabf(idx)=[];
    end
end




