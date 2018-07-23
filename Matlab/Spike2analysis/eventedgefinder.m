function up_idx = eventedgefinder(data,threshold,fs,minStateLength,mingap,up,Num)

if size(data,2)==1;data=data';end

Wn = [(4000/(fs/2))];
[B,A] = butter(2,Wn);
data = filtfilt(B,A,data);


switch up
    case 1
        A=find(data>threshold);
        idx= find(diff(A)~=1);
        start_idx = [1 idx+1];
        end_idx= [idx numel(A)];
        up_idx= A([start_idx ; end_idx]');
        up_idx= filterShortStates(up_idx,fs,minStateLength,mingap);
    case 0
        A=find(data<threshold);
        idx= find(diff(A)~=1);
        start_idx = [1 idx+1];
        end_idx= [idx numel(A)];
        up_idx= A([start_idx ; end_idx]');
        up_idx= filterShortStates(up_idx,fs,minStateLength,mingap);
end

if Num>0
    up_idx=up_idx(1:Num,:);
end

end

function states= filterShortStates(states,fs,minStateLength,mingap)

D=states(:,2)-states(:,1);
r= D>round(minStateLength*fs);
states=states(r,:);

for i=1:size(states,1)-1
    if (states(i+1,1)-states(i,2))<round(mingap*fs)
        states(i+1,1)=0;
        states(i,2)=0;
    end
end
states=reshape(states',1,2*size(states,1));
r=states~=0;
states=states(r);
states=reshape(states,2,length(states)/2)';

end

