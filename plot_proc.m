% plot and data process
%clear all;clc;close all
if (0)%save relaxation data
        mmxstart=mmx(end,:);mmystart=mmy(end,:);mmzstart=mmz(end,:);
        save('startm_TT160_natom2000_pc4','mmxstart','mmystart','mmzstart')
end
if (0)
    dwplotstep=3;
    figure%initial magnetization
    hold on
    for ct=1:dwplotstep:natom
        if mark_(ct)==1
            quiver3(0,loc_(ct)*1e9,0,mmx(1,ct),mmy(1,ct),mmz(1,ct),'r');
        else
            quiver3(0,loc_(ct)*1e9,0,mmx(1,ct),mmy(1,ct),mmz(1,ct),'b');
        end
    end
    xlabel('x axis');ylabel('y axis');zlabel('z axis');
    xlim([-1 1]);ylim([-5 450]);zlim([-2 2]);
    view(-27,20)
    
    figure%magnetization after relaxation
    hold on
    for ct=1:dwplotstep:natom
        if mark_(ct)==1
            quiver3(0,loc_(ct)*1e9,0,mmx(round(end/2),ct),mmy(round(end/2),ct),mmz(round(end/2),ct),'r');
        else
            quiver3(0,loc_(ct)*1e9,0,mmx(round(end/2),ct),mmy(round(end/2),ct),mmz(round(end/2),ct),'b');
        end
    end
    xlabel('x axis');ylabel('y axis');zlabel('z axis');
    xlim([-1 1]);ylim([-5 850]);zlim([-2 2]);
    view(-27,20)
    
    figure%final magnetization
    hold on
    for ct=1:dwplotstep:natom
        if mark_(ct)==1
            quiver3(0,loc_(ct)*1e9,0,mmx(end,ct),mmy(end,ct),mmz(end,ct),'r');
        else
            quiver3(0,loc_(ct)*1e9,0,mmx(end,ct),mmy(end,ct),mmz(end,ct),'b');
        end
    end
    xlabel('x axis');ylabel('y axis');zlabel('z axis');
    xlim([-1 1]);ylim([-5 850]);zlim([-1 1]);
    view(-45,53)

    if (0)
        figure%begin magnetization, mx,my,mz
        plotstep=2;
        plot(loc_(1:plotstep:end)*1e9,mmz(round(end/2),1:plotstep:end),...
            loc_(1:plotstep:end)*1e9,mmy(round(end/2),1:plotstep:end),...
            loc_(1:plotstep:end)*1e9,mmz(round(end/2),1:plotstep:end),'linewidth',1);
        xlabel('location');ylabel('m');
        legend('mx','my','mz')
        xlim([-5 850]);ylim([-1 1]);
    end
    if (0)
        figure%final magnetization, mx,my,mz
        plotstep=2;
        plot(loc_(1:plotstep:end)*1e9,mmz(end,1:plotstep:end),...
            loc_(1:plotstep:end)*1e9,mmy(end,1:plotstep:end),...
            loc_(1:plotstep:end)*1e9,mmz(end,1:plotstep:end),'linewidth',1);
        xlabel('location');ylabel('m');
        legend('mx','my','mz')
        xlim([-5 450]);ylim([-1 1]);
    end
    if (0)
        figure%final magnetization, mx
        plotstep=2;
        hold on
        plot(loc_(1:plotstep:end)*1e9,mmx(end,1:plotstep:end),'linewidth',1);
        xlabel('location');ylabel('mx');
        xlim([-5 450]);ylim([-1 1]);
    end
    if (0)
        figure%final magnetization, my
        plotstep=2;
        hold on
        plot(loc_(1:plotstep:end)*1e9,mmy(end,1:plotstep:end),'linewidth',1);
        xlabel('location');ylabel('my');
        xlim([-5 450]);ylim([-1 1]);
    end
    if (1)
        figure%final magnetization, mz
        plotstep=2;
        hold on
        plot(loc_(1:plotstep:end)*1e9,mmz(end,1:plotstep:end),'linewidth',1);
        xlabel('location');ylabel('mz');
        xlim([-5 450]);ylim([-1 1]);
    end
    if (1)
        figure%mx, begin vs final
        plotstep=2;
        hold on
        plot(loc_(1:plotstep:end)*1e9,mmx(round(end/2),1:plotstep:end),'*r','linewidth',2);
        plot(loc_(1:plotstep:end)*1e9,mmx(end,1:plotstep:end),'ob','linewidth',1);
        xlabel('location');ylabel('mx');
        legend('begin','final')
        xlim([-5 450]);ylim([-1 1]);
    end
    if (1)
        figure%my, begin vs final
        plotstep=2;
        hold on
        plot(loc_(1:plotstep:end)*1e9,mmy(round(end/2),1:plotstep:end),'*r','linewidth',2);
        plot(loc_(1:plotstep:end)*1e9,mmy(end,1:plotstep:end),'ob','linewidth',1);
        xlabel('location');ylabel('my');
        legend('begin','final')
        xlim([-5 450]);ylim([-1 1]);
    end
    if (1)
        figure%mz, begin vs final
        plotstep=2;
        hold on
        plot(loc_(1:plotstep:end)*1e9,mmz(round(end/2),1:plotstep:end),'r','linewidth',2);
        plot(loc_(1:plotstep:end)*1e9,mmz(end,1:plotstep:end),'b','linewidth',1);
        xlabel('location');ylabel('mz');
        legend('begin','final')
        xlim([-5 450]);ylim([-1 1]);
    end
    
        plotatom=498;
    figure%dynamics
    plot(t*1e12,mmx(:,plotatom),t*1e12,mmy(:,plotatom),t*1e12,mmz(:,plotatom),'linewidth',2);
    xlabel('time(ps)');ylabel('m');legend('mx','my','mz');
    
    figure;%3d plot
    plot3(mmx(:,plotatom),mmy(:,plotatom),mmz(:,plotatom))
    xlabel('mx');ylabel('my');ylabel('mz');
    
    %FFT
%     [~,loctmp]=max(mmx(:,plotatom));
%     t(loctmp)
%     nT0=size(mmx(:,plotatom));
%     FFT_module(nT0,runTime0,y,rminit,rmlast,plotfft)
%     clear loctmp
end

if (1)%sweepTT
    TT_=[108,130,149,217,260];
    %TT_=[10:10:80];
    debug=1;
    machinesel=2;%1:optiplex7040 2:laptop
    for ctTTproc=1:size(TT_,2)
        TT=TT_(ctTTproc);
        switch machinesel
            case 1
                foldname=sprintf('C:\\Users\\Brian\\Desktop\\tmp\\2017Sehyeok\\dmi100\\TT%d',TT);
                cd(foldname)
                %datset=[0,1,2,4,8,16,32,64,128,256];
                datset=[100,200,300,400,500,600];
            case 2
                %foldname=sprintf('C:\\Users\\ACER\\Desktop\\tmp\\2017Sehyeok\\dmi\\dmi100_smallH\\TT%d',TT);
                %foldname=sprintf('C:\\Users\\ACER\\Desktop\\tmp\\2017Sehyeok\\dmi\\TT%d\\dmi0',TT);
                foldname='C:\Users\ACER\Desktop\tmp\2017Sehyeok\dmi\natom1000\dmipc4\TT130\dmi0';
                cd(foldname)
                %datset=[1,2,4,8,16,32,64,128,256,512];
                datset=[2:2:64];
                %datset=[223,282,342,402,517,573,629];
        end
        dwspeed_=zeros(size(datset,2),1);
        dw_width=zeros(size(datset,2),1);
        for ctproc=1:size(datset,2)
            close all
            datname=sprintf('finalgpu_hext%d.mat',datset(ctproc));
            %datname=sprintf('finalgpu_dmi%d.mat',datset(ctproc));
            load(datname)
            datset(ctproc)
            plotstep=2;
            if debug
                figure%my, begin vs final
                hold on
                plot(loc_(1:plotstep:end)*1e9,mmy(1,1:plotstep:end),'*r','linewidth',2);
                plot(loc_(1:plotstep:end)*1e9,mmy(end,1:plotstep:end),'ob','linewidth',1);
                xlabel('site(nm)');ylabel('my');
                legend('begin','final')
                xlim([-5 450]);ylim([-1 1]);
            end
            [~,locstart]=max(abs(mmy(1,1:plotstep:end)));
            [~,locend]=max(abs(mmy(end,1:plotstep:end)));
            
            dwspeed_(ctproc)=plotstep*(loc_(locend)-loc_(locstart))/runtime2;
            if (0)%dw width
                if debug
                    figure%my, begin vs final
                    hold on
                    %plot(loc_(1:plotstep:end)*1e9,mmx(end,1:plotstep:end),'r','linewidth',1);
                    %plot(loc_(1:plotstep:end)*1e9,mmy(end,1:plotstep:end),'b','linewidth',1);
                    plot(loc_(1:plotstep:end)*1e9,mmz(end,1:plotstep:end),'*m','linewidth',1);
                    xlabel('site(nm)');ylabel('m');
                    %legend('mx','my','mz')
                    xlim([-5 850]);ylim([-1 1]);
                end
                tmp=mmz(end,1:plotstep:end);
                for ctdwwidth=2:(size(tmp,2)-1)
                    if abs(tmp(ctdwwidth))>0.9 && abs(tmp(ctdwwidth+1))<0.9
                        dwstart=ctdwwidth;
                    end
                    if abs(tmp(ctdwwidth))<0.9 && abs(tmp(ctdwwidth+1))>0.9
                        dwend=ctdwwidth+1;
                    end
               
                end

            end
%                 plotatom=plotstep*(locend-50);
%     figure%dynamics
%     plot(t*1e12,mmx(:,plotatom),t*1e12,mmy(:,plotatom),t*1e12,mmz(:,plotatom),'linewidth',2);
%     xlabel('time(ps)');ylabel('m');legend('mx','my','mz');
        end
    end
end
if (0)%tstep
    clear all;clc
    foldnam='C:\Users\ACER\Desktop\tmp\2017Sehyeok\dmi\natom2000\dmipc4\TT160\tstep\dmi0_hext400';
    cd(foldnam)
    
    plotatom=2;
    figure%dynamics
    hold on
    load('finalgpu_tt8.mat')
    plot(t*1e12,mmx(:,plotatom),'r','linewidth',3);
    load('finalgpu_tt4.mat')
    plot(t*1e12,mmx(:,plotatom),'b','linewidth',2);
    load('finalgpu_tt2.mat')
    plot(t*1e12,mmx(:,plotatom),'k','linewidth',1);
    xlabel('time(ps)');ylabel('mz');legend('0.8fs','0.4fs','0.2fs');

    plotatom=996;
    figure%dynamics
    hold on
    load('finalgpu_tt8.mat')
    plot(t*1e12,mmx(:,plotatom),'r','linewidth',3);
    load('finalgpu_tt4.mat')
    plot(t*1e12,mmx(:,plotatom),'b','linewidth',2);
    load('finalgpu_tt2.mat')
    plot(t*1e12,mmx(:,plotatom),'k','linewidth',1);
    xlabel('time(ps)');ylabel('mz');legend('0.8fs','0.4fs','0.2fs');
    
%     plotatom=498;
%     figure%dynamics
%     hold on
%     load('finalgpu_tt4.mat')
%     plot(t*1e12,mmx(:,plotatom),'r','linewidth',3);
%     load('finalgpu_tt2.mat')
%     plot(t*1e12,mmx(:,plotatom),'b','linewidth',2);
%     load('finalgpu_tt1.mat')
%     plot(t*1e12,mmx(:,plotatom),'k','linewidth',1);
%     xlabel('time(ps)');ylabel('mz');legend('0.4fs','0.2fs','0.1fs');
end
if (0)%compare with previous
    figure%dynamics
    plot(t*1e12,mmx(:,2),'r','linewidth',3);
    hold on
    load('t1.mat')
    plot(t*1e12,mmx(:,2),'b','linewidth',1);
    xlabel('time(ps)');ylabel('m');
    legend('new','old');
end
if (0)%compare cpu and gpu result
    figure%dynamics
    plot(t*1e12,mmx(:,2),'r','linewidth',3);
    hold on
    addpath('../cpu')
    load('finalcpu.mat')
    plot(t*1e12,mmx(:,2),'b','linewidth',1);
    xlabel('time(ps)');ylabel('m');
    legend('gpu','cpu');
end
if (0)%compare diff tstep
    figure%dynamics
    plot(t*1e12,mmx(:,2),'r','linewidth',3);
    hold on
    plot(t*1e12,mmx(:,2),'k','linewidth',2);
    plot(t*1e12,mmx(:,2),'b','linewidth',1);
    plot(t*1e12,mmx(:,2),'m','linewidth',1);
    xlabel('time(ps)');ylabel('m');
    legend('t1','t2','t3','t4');
end
if (0)%generate movies
    close all
    plotstep=2;
    ctplotend=469;
    for ctplot = 1:ctplotend
        figure%my, begin vs final
        titlename=sprintf('t=%fps',(ctplot-1)/ctplotend*runtime2*1e12);
        title(titlename);
        hold on
        plot(loc_(1:plotstep:end)*1e9,mmx(1+(ctplot-1)*size(mmy,1)/ctplotend,1:plotstep:end),'k','linewidth',1);
        plot(loc_(1:plotstep:end)*1e9,mmy(1+(ctplot-1)*size(mmy,1)/ctplotend,1:plotstep:end),'r','linewidth',1);
        plot(loc_(1:plotstep:end)*1e9,mmz(1+(ctplot-1)*size(mmy,1)/ctplotend,1:plotstep:end),'b','linewidth',1);
        xlabel('site(nm)');ylabel('m');
        legend('mx','my','mz')
        xlim([-5 420]);ylim([-1 1]);
        M(ctplot) = getframe(gcf);
    end
    clear ctplot
    figure
    movie(M)
    v=VideoWriter('myavi.avi')
    open(v)
    writeVideo(v,M)
end

if (0)
    figure%final magnetization
    hold on
    for ct=480:520
        if mark_(ct)==1
            q=quiver3(0,loc_(ct)*1e9,0,mmx(end,ct),mmy(end,ct),mmz(end,ct),'r');
        else
            q=quiver3(0,loc_(ct)*1e9,0,mmx(end,ct),mmy(end,ct),mmz(end,ct),'b');
        end
        q.MaxHeadSize = 0.5;
    end
    xlabel('x axis');ylabel('y axis');zlabel('z axis');
    xlim([-1 1]);ylim([192 208]);zlim([-1 1]);
    view(-103,52)
end
if (0)
    figure% magnetization
    hold on
    for ct=480:540
        if mark_(ct)==1
            q=quiver3(0,loc_(ct)*1e9,0,mmx(end,ct),mmy(end,ct),mmz(end,ct),'r');
            q=quiver3(0.5,loc_(ct)*1e9,0,mmx(round(end/2),ct),mmy(round(end/2),ct),mmz(round(end/2),ct),'r');
        else
            q=quiver3(0,loc_(ct)*1e9,0,mmx(end,ct),mmy(end,ct),mmz(end,ct),'b');
            q=quiver3(0.5,loc_(ct)*1e9,0,mmx(round(end/2),ct),mmy(round(end/2),ct),mmz(round(end/2),ct),'b');
        end
        q.MaxHeadSize = 0.5;
    end
    xlabel('x axis');ylabel('y axis');zlabel('z axis');
    xlim([-1 1]);ylim([180 220]);zlim([-1 1]);
    view(-103,52)
end