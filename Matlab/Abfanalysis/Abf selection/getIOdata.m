function [t] = getIOdata(file)
    ss=load('D:\Morphys\Data\Electrophysiology\SetupSettings\Setupsettings_INF.mat');
    ss=ss.obj;
    load(file);
    
    if exist('a', 'var')
        sweep=a.getchannel.getin('signal', 'primary').getsweep;
    elseif exist('a1', 'var')
        sweep1 = a1.getchannel.getin('signal','primary').getsweep ;
        sweep2 = a2.getchannel.getin('signal','primary').getsweep ;
        sweep=[sweep1, sweep2([2 1 3:end])];
    elseif exist('c', 'var')
        sweep=c.getin('signal', 'primary').getsweep;
    end
    
    % find current injection epoch and assign aps to sweep
    for step = 1:length(sweep(1).getepoch)
        if sweep(1).getepoch(step).stepdiff < 0 && abs(sweep(1).getepoch(step).stepdiff + sweep(1).getepoch(step+1).stepdiff) < 5
            break
        end
    end
    
    sweep([sweep.nrofepochs]==1)=[];

    %get sweep data
    t=table;
    t.pA=[sweep.getepoch(step).amplitude]';
    t.nrofaps=[sweep.getepoch(step).nrofaps]';
    
    
    for i=1:length(sweep)
        if t.nrofaps(i)>3
            t.freq(i)=mean([sweep(i).getepoch(step).getap(4:end).freq]);
            t.isi(i)=mean([sweep(i).getepoch(step).getap(4:end).isi]);
        end
    end
    t.sag=[sweep.getepoch(step).sag]';
    for i=1:length(sweep)
        if t.pA(i)<0
            vmbase=sweep(i).getepoch(step-1).steadystate;
            vmresponse=sweep(i).getepoch(step).vstep;
            PkDeflect = vmbase - vmresponse;
            t.vmresponse(i)=vmresponse;
            t.sag_perc(i)=t.sag(i)/PkDeflect;
        else
            t.vmresponse(i)=NaN;
            t.sag_perc(i)=NaN;
        end
    end
end
