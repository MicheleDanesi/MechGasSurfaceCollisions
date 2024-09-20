function harmolen_plot_hits(nn,pausa,pars,in,out,hits) 
    % display hit diagram
    
    hit_type={'crest-1','2-crests','crest-2','face-1','2-faces','face-2'};
    rp = pars.r;
    
    % plot zig-zag
    gg=3;
    xg=[-gg:gg]*pars.gCos;
    zg=[ones(1,gg);zeros(1,gg)];zg=[zg(:); 1]'*pars.gSin;
    
    figure(1)
    f1pos = get(gcf,'Position');
    close(1)
    figure(1)
    set(gcf,'Position',f1pos)
    %hold off
    
    patch([xg(1) xg xg(end)],[1 zg 1],[0.2 0.2 0.2],'edgecolor','none','FaceAlpha',0.2)
    hold on
    
    plot(xg,zg,'k')
    axis equal
    axis([-2 2 -3*rp 1.2*pars.gSin])
    %hold on
    
    plot([-1/sqrt(2) sqrt(2)],-rp*[1 1],'k:')
    
    % plot particles
    [xxp0,zzp0,~]=cylinder(rp,50);
    
    if numel(nn)==1
        if nn==0
            nn=(1:size(in.v,1)); % this plots all hits 
        else 
            nn=(1:nn);
        end
    end 
    
    hh=[];
    for iii=1:numel(nn) %min(20,numel(nn))
        
        % get hit index and type
        ii=nn(iii);
        % delete previous if exists
        if sum(ishandle(hh))
            for jj=1:numel(hh)
                delete(hh(jj));
            end
        end
        
        parsp = [];
        
        % start position
        parsp.acolor = [0 0 1];
        hh=plot_particle(rp,in.p.x(ii),in.p.z(ii),in.v(ii,:),parsp);
        
        % hits position
        parsp.color = 0.4*[1 1 1];
        k=0;
        for jj=1:numel(hits)
            fii = find(hits(jj).i==ii);
            if ~isempty(fii)
                parsp.color = 1-(1-parsp.color)*0.8;
                parsp.acolor = parsp.color;
                k=k+1;
                hh1=plot_particle(rp,hits(jj).p.x(fii),hits(jj).p.z(fii),hits(jj).v(fii,:),parsp);
                hh = [hh; hh1];
            end
        end
        
        % exit position
        parsp.color = 'm';
        parsp.mcolor = 'm';
        parsp.acolor = 'm';
        parsp.marker = '.';
        hh1=plot_particle(rp,out.p.x(ii),out.p.z(ii),out.v(ii,:),parsp);
    
        hh = [hh; hh1];
        
        ftitl = ['particle #' num2str(ii) ' - ' num2str(k) ' hits'];

        if out.n(ii)
            if out.n(ii)==1
                ftitl = ftitl(1:end-1);
            end
            ftitl = [ftitl ' ['];
            for kk=1:out.n(ii)
                fii = find(hits(kk).i==ii);
                ftitl = [ftitl hit_type{hits(kk).w(fii)} ':'];
            end
            ftitl = [ftitl(1:end-1) ']'];
        end
        if sign(in.p.x(ii)<0)
            ftitl = [ftitl '*'];
        end

        title(ftitl)
        drawnow
        pause(pausa)
    end
    hold off
end
