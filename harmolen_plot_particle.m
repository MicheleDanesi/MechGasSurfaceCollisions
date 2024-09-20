function hh=harmolen_plot_particle(rp,xp,zp,vp,parsp)
    % draw a particle with velocity arrow

    % set default parameters
    if nargin<5
        parsp = [];
    end
    if ~isfield(parsp,'color')
        parsp.color = 'b';
    end
    if ~isfield(parsp,'marker')
        parsp.marker = '.';
    end
    if ~isfield(parsp,'mcolor')
        parsp.mcolor = parsp.color;
    end
    if ~isfield(parsp,'alpha')
        parsp.alpha = 0.2;
    end
    if ~isfield(parsp,'acolor')
        parsp.acolor = parsp.color;
    end
    
    % define points
    [xxp0,zzp0,~]=cylinder(rp,500);
    xxp=xxp0(1,:)+xp;zzp=zzp0(1,:)+zp; % particle
    ap=repmat([xp+vp(1);zp+vp(3)],1,3)+0.1*[[sin(pi+pi/3) cos(pi+pi/3);-cos(pi+pi/3) sin(pi+pi/3)]*[vp(1);vp(3)] [0;0] [sin(-pi/3) cos(-pi/3);-cos(-pi/3) sin(-pi/3)]*[vp(1);vp(3)]]; % arrow
    
    % do plot
    d=[[xp;xp+10*max(rp,1)*vp(1)/norm(vp([1 3]))],[zp;zp+10*max(rp,1)*vp(3)/norm(vp([1 3]))]];
    %norm(d(2,:)-d(1,:))
    hh=plot(d(:,1),d(:,2),'k:');
    if ischar(parsp.color)
        hh(end+1:end+2,1)=plot(xxp,zzp,parsp.color,xp,zp,[parsp.mcolor parsp.marker]);
    else
        hh(end+1,1)=plot(xxp,zzp,'color',parsp.color);
        hh(end+1,1)=plot(xp,zp,parsp.marker,'color',parsp.mcolor);
    end
    hh(end+1,1)=patch(xxp,zzp,parsp.color,'edgecolor','none','FaceAlpha',parsp.alpha);
    
    hh(end+1,1)=line([xp,xp+vp(1),ap(1,:)],[zp,zp+vp(3),ap(2,:)],'color',parsp.acolor);

end