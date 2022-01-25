
path2saveplots='/Users/natalia/Documents/Research/Projects_active/AllenInstitute/Data Analysis/e-phys/e-phys_all_labs/Figures';
%features=fields(data);
features={'upstroke','downstroke'};
bin2plot=fields(data.upstroke);
yaslabel={'AP upstroke, mV/ms', 'AP downstroke, mV/ms','AP peak voltage, mV',...                   }
'AP half-width, s','AP threshold, mV','upstroke/downstroke ratio'}
% set colors for t-types
%%
% for i=1:length(t_types)
%     rgb_colors3(i,1:3)=uisetcolor;
% end
%%

          mean_data=table;
          std_data=table;
for n=1:length(features)

   fig1=figure('Position',[182,377,1000,420]);
   axes1 = axes('Parent',fig1,'Position',[0.1 0.1 .7 0.8]);
   
    for k=1:size(t_types)
          allbin_tempdata=[];

        for i=2:7
            idx=[];
            tempdata=[];
            idx=strcmp(data.(features{n}).t_type,t_types(k));
            tempdata=table2array(data.(features{n})(idx,bin2plot{i}));
            var_names{i}=sprintf('%s','mean_',bin2plot{i});
            data.(features{n}).(var_names{i})(idx)=mean(tempdata,'omitnan');
            mean_data.t_type{k}=t_types{k};
            mean_data.(var_names{i})(k)=mean(tempdata,'omitnan');
            t_types{k}
            size(~isnan(tempdata))
            std_data.t_type{k}=t_types{k};
            std_data.(var_names{i})(k)=std(tempdata,'omitnan');
            data_all_mean.(features{n})=mean_data;
            data_all_std.(features{n})=std_data;
            allbin_tempdata(:,i-1)=tempdata;
            
        end
      if size(allbin_tempdata,1)>1
      %p1=stdshade(allbin_tempdata,[0.1],rgb_colors(k,:));
     
      text( 6,nanmean(allbin_tempdata(:,end)),t_types{k});
      else
      p1=plot(allbin_tempdata,'-', 'LineWidth',2, 'Color',...
       rgb_colors(k,:), 'MarkerFaceColor', rgb_colors(k,:)); 
      text(6,allbin_tempdata(1,end-1),t_types{k});
     
      end
      if  k==1
              p=p1;
      else 
              p=[p p1];
      end
      hold on
      
    end
    
    title(features{n})
    set(axes1,'XTick',[1:6],'XTickLabel',...
    {'First AP','0-10 Hz','10-20 Hz','20-30 Hz','30-40 Hz','40-50 Hz'});
    legend(p, t_types, 'Location','WestOutside')
    axes1.FontName='Arial';
    axes1.FontSize= 14;
    ylabel(yaslabel{n});
    name=sprintf('%s',features{n},'_fig_Huib_new.pdf');
    filename=fullfile(path2saveplots,name);

    exportgraphics(fig1,filename,'Resolution',300)
    %print(filename,'-depsc');
    xlabel('Instantaneous frequency, Hz');
end
%% plot all freq bins

%color1={'#0FF4F4','#3377EF','#BA6315','#0698EA','#0ADD69','#111010',...
  %  '#2E3192','#8C288F','#EA8C68','#F40C3E','#87EF0F','#F24F0F','#DC088C',...
   % '#301EE0','#97268F','#2E3976'};
f=figure;
axes2 = axes('Parent',f);
i=1;
for  k=1:size(t_types)
   plot(data_all_mean.(features{i}){k,2:7},'o-', 'LineWidth',2, 'Color',...
       colors{k,2}, 'MarkerFaceColor', colors{k,2})
   hold on
   legend(t_types, 'Location','EastOutside')
   set(axes2,'XTick',[1:6],'XTickLabel',...
     {'First AP','0-10 Hz','10-20 Hz','20-30 Hz','30-40 Hz','40-50 Hz'});
end

%%
% make a violin scatter plot


for n=1:length(features)
   
    x=categorical(data.(features{n}).t_type, t_types);
    for nn=1:size(data.(features{n}).t_type,1)
        idx=strcmp(data.(features{n}).t_type{nn},t_types);
        if sum(idx)~=0
        c(nn,:)=rgb_colors(idx,:);
        else
        c(nn,:)=[0 0 0];
        end
    end
    % make all Mansvelder data (from the row number 259) black
    for i=259:298;
        c(i,:)=[0 0 0];
    end
    for i=2:6
    figure
    y=[];
    y=data.(features{n}).(bin2plot{i});
    y2=data.(features{n}).(var_names{i});
    swarmchart(x,y,30,c, 'LineWidth',1);
    %colormap(colors(:,2))
    hold on
    scatter(x,y2,'k_','LineWidth',3)
    titlename=sprintf('%s_%s',features{n},bin2plot{i});
    title(titlename);
    ax=gca;
    ax.FontName='Arial';
    ax.FontSize= 12;
    ylabel(yaslabel{n});
    filename=fullfile(path2saveplots,titlename);
    %print(filename,'-depsc');
    %p = kruskalwallis(data.(features{n}).(bin2plot{i}),data.(features{n}).t_type,'off')
    end
end
%% scatter of sag and upstroke
%t_b=t_types(Blue);
 %t_b=t_types(Pax);
 %t_b=t_types(Green);
 %t_b=t_types([Green; Blue]);
 %t_b=t_types([Orange; Pink]);
 t_b=t_types;
 colors2plot=rgb_colors;
figure('Position', [100 100 1000 700]);
k=1;
for i=1:size(t_b,1)
    idx3plot=[];
    idx3plot= strcmp(t_b{i},data.upstroke.t_type);
    hold on
    scatter(data.sag.Sag(idx3plot),data.upstroke.s30to40Hz(idx3plot),100,colors2plot(i,:),'filled')
end
legend(t_b,'Location','EastOutside')
ylabel('AP upstroke at 30-40 Hz, mV/ms')
xlabel('Sag ratio')
%  figure; scatter(data.sag.Sag,data.downstroke.s20to30Hz,100,c,'filled')
% figure; scatter(data_norm.sag.Sag,data_norm.upstroke.s20to30Hz,100,c,'filled')
% figure; scatter(data_norm.sag.Sag,data_norm.downstroke.s20to30Hz,100,c,'filled')
ylim([0 500]);
xlim([0 0.6]);
filename=fullfile(path2saveplots,'Upstroke_vs_sag_all_VIP');
print(filename,'-depsc');
 