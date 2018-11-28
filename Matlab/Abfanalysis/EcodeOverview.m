function EcodeOverview
% Prints out a pdf overview of the ecodes recorded in a cell
%   Morphys code is necessary for this to function
%   Bioinformatics toolbox is necessary for the suptitle function to work
% Instruction: insert the filename to be saved (use the cell name)
%              Select the abf files % setupsettings
%              Figures and PS file with overview will be saved in location
%              folder


[fn, path] = uiputfile('*.ps');
fullname=fullfile(path, fn);
if ~exist([path 'figs\'], 'dir'), mkdir([path 'figs\']); end

b = Abfbatch('gui');
pronames={b.getabf.proname};


% make a plot of the RAMP protocol
b.getabf(contains(pronames, 'RAMP')).plot;
suptitle([fn(1:end-3) '\_RAMP']);
savefig([path 'figs\' fn(1:end-3) '_RAMP.fig']);
print(fullname, '-dpsc', '-r300');
close all

% make a plot of the LSCOARSE protocol
b.getabf(contains(pronames, 'LSCOARSE')).plot;
suptitle([fn(1:end-3) '\_LSCOARSE']);
print(fn, '-dpsc', '-r300', '-append');
savefig([path 'figs\' fn(1:end-3) '_LSCOARSE.fig']);
close all

% make a plot of the LSFINEST protocol
aplist=b.getabf(contains(pronames, 'LSFINEST')).apsweeps;
%plot the first sweep that has more than 5 AP's
sweep=find(aplist>4,1);
if isempty(sweep), sweep=numel(aplist); end
b.getabf(contains(pronames, 'LSFINEST')).selectsweep(sweep).plot;
suptitle([fn(1:end-3) '\_LSFINEST']);
print(fullname, '-dpsc', '-r300', '-append');
savefig([path 'figs\' fn(1:end-3) '_LSFINEST.fig']);
close all

% make a plot of the SSFINEST protocol
b.getabf(contains(pronames, 'SSFINEST')).plot;
xlim([1050 1070])
suptitle([fn(1:end-3) '\_SSFINEST']);
print(fullname, '-dpsc', '-r300', '-append');
savefig([path 'figs\' fn(1:end-3) '_SSFINEST.fig']);
close all

% make a plot of the SSFINEST protocol
b.getabf(contains(pronames, 'SSFINEST')).plot;
xlim([1050 1070])
suptitle([fn(1:end-3) '\_SSFINEST']);
print(fullname, '-dpsc', '-r300', '-append');
savefig([path 'figs\' fn(1:end-3) '_SSFINEST.fig']);
close all

% make a plot of the CHIRP protocol
b.getabf(contains(pronames, 'CHIRP')).plot;
suptitle([fn(1:end-3) '\_CHIRP']);
print(fullname, '-dpsc', '-r300', '-append');
savefig([path 'figs\' fn(1:end-3) '_CHIRP.fig']);
close all

% make a plot of the TRIPLE protocol
aplist=b.getabf(contains(pronames, 'TRIPLE')).apsweeps;
sweep=find(aplist<3,1);
if isempty(sweep), sweep=numel(aplist); else, sweep=[sweep, sweep-1]; end
b.getabf(contains(pronames, 'TRIPLE')).selectsweep(sweep).plot;
suptitle([fn(1:end-3) '\_TRIPLE']);
xlim([400 500])
print(fullname, '-dpsc', '-r300', '-append');
savefig([path 'figs\' fn(1:end-3) '_TRIPLE.fig']);
close all

% make a plot of the C1NSD2SHORT protocol
b.getabf(contains(pronames, 'C1NSD2SHORT')).plot;
suptitle([fn(1:end-3) '\_C1NSD2SHORT']);
print(fullname, '-dpsc', '-r300', '-append');
savefig([path 'figs\' fn(1:end-3) '_C1NSD2SHORT.fig']);
close all


end
