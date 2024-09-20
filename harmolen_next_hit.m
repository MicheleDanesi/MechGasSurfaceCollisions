function [position_out, velocity_out, whts] = harmolen_next_hit(position_in, velocity_in, whts, pars)
    
    vTan = velocity_in(:,3)./velocity_in(:,1);
    vAng = atan(vTan);
    vAng(velocity_in(:,1)==0)=pi/2;
    
    % initialize some arrays
    velocity_out = velocity_in;
    position_out.x = position_in.x;
    position_out.z = position_in.z;
    xc_l = nan(size(velocity_in,1),1);
    xc_r = nan(size(velocity_in,1),1);
    zc_l = nan(size(velocity_in,1),1);
    zc_r = nan(size(velocity_in,1),1);
    
    % find hit positions on wall faces (check both sides)
    a = vTan;
    b = -ones(size(xc_l));
    c_r = -(position_in.x+pars.r*pars.gSin).*vTan + (position_in.z+pars.r*pars.gCos); % trajectory of particle's right tangent point (a*x+b*z+c_r=0)
    c_l = -(position_in.x-pars.r*pars.gSin).*vTan + (position_in.z+pars.r*pars.gCos); % trajectory of particle's left tangent point  (a*x+b*z+c_l=0)
    
    xf_l = (b*pars.walls.left.c-c_l*pars.walls.left.b)./(a*pars.walls.left.b-b*pars.walls.left.a)+pars.r*pars.gSin; % center at intersection...
    xf_r = (b*pars.walls.right.c-c_r*pars.walls.right.b)./(a*pars.walls.right.b-b*pars.walls.right.a)-pars.r*pars.gSin; 
    
    zf_l = (c_l*pars.walls.left.a-a*pars.walls.left.c)./(a*pars.walls.left.b-b*pars.walls.left.a)-pars.r*pars.gCos; 
    zf_r = (c_r*pars.walls.right.a-a*pars.walls.right.c)./(a*pars.walls.right.b-b*pars.walls.right.a)-pars.r*pars.gCos;


    % discard solutions out of grooves z range
    zf_l(zf_l+pars.r*pars.gCos>pars.gSin-pars.r*(1./pars.gCos-pars.gCos))=nan;
    zf_r(zf_r+pars.r*pars.gCos>pars.gSin-pars.r*(1./pars.gCos-pars.gCos))=nan;
    zf_l(zf_l+pars.r*pars.gCos<0)=nan;
    zf_r(zf_r+pars.r*pars.gCos<0)=nan;
    
    % find left-crest collision position
    A = (1+vTan.^2);
    B = (position_in.x.*vTan-position_in.z).*vTan;
    C = position_in.z.^2+(position_in.x.*vTan).^2-pars.r.^2-2*position_in.x.*position_in.z.*vTan; 
    D = B.^2-A.*C;
    xc_l(D>0) = (B(D>0)-sign(velocity_in(D>0,1)).*sqrt(D(D>0)))./A(D>0);
    zc_l(D>0) = position_in.z(D>0)+(xc_l(D>0)-position_in.x(D>0)).*vTan(D>0);
    
    % find right-crest hit collision position
    %A=(1+vTan.^2); no need to update A
    B = (position_in.x.*vTan-position_in.z).*vTan+2*pars.gCos;
    C = position_in.z.^2+(position_in.x.*vTan).^2-pars.r.^2-2*position_in.x.*position_in.z.*vTan+4*pars.gCos^2;
    D = B.^2-A.*C;
    xc_r(D>0) = (B(D>0)-sqrt(D(D>0)))./A(D>0);
    zc_r(D>0) = position_in.z(D>0)+(xc_r(D>0)-position_in.x(D>0)).*vTan(D>0);
    
    % discard unacceptable solutions for crest hits
    zc_l(zc_l>=0)=nan; % negative z condition should be enough
    zc_r(zc_r>=0)=nan; 
    

    % get first hit type
    htsz=[zc_l zc_r zf_l zf_r]; 
        
    % discard inconsistent solutions based on velocity direction
    for ii = 1:4
        htsz(sign(htsz(:,ii)-position_in.z)~=sign(velocity_in(:,3)),ii)=nan; 
    end

    % discard solutions based on last hit side
    whts = real([whts==1|whts==4 whts==3|whts==6]);
    whts(whts==1)=nan;
    htsz=htsz+[whts whts];

    whts_up = htsz==minnonan(htsz,[],2); %,'omitnan');
    whts_down = htsz==maxnonan(htsz,[],2); %,'omitnan');
    
    whts = whts_up;
    whts(sign(velocity_in(:,3))<0,:) = whts_down(sign(velocity_in(:,3))<0,:);
    
    cc_hits = (zc_r<0) & (zc_l==zc_r); % double crest hits
    lc_hits=whts(:,1);
    rc_hits=whts(:,2);
    lf_hits=whts(:,3);
    rf_hits=whts(:,4);
    ff_hits = (zf_l==zf_r) & ~cc_hits; % groove hits (double face)
    
    whts = [lc_hits cc_hits rc_hits lf_hits ff_hits rf_hits]; % hit types
    whts = sum(whts.*repmat([1:6],size(whts,1),1),2);

    
    % handle face bounce
    velocity_out(lf_hits,:) = (pars.Rl^-1*(((pars.Px+pars.Py-pars.Pz)*pars.Rl*velocity_in(lf_hits,:)').*harmolen_CoR(pars.Rl*velocity_in(lf_hits,:)',pars.CoR)))'; 
    position_out.x(lf_hits)=xf_l(lf_hits);
    position_out.z(lf_hits)=zf_l(lf_hits);
    
    velocity_out(ff_hits|rf_hits,:) = (pars.Rr^-1*(((pars.Px+pars.Py-pars.Pz)*pars.Rr*velocity_in(ff_hits|rf_hits,:)').*harmolen_CoR(pars.Rr*velocity_in(ff_hits|rf_hits,:)',pars.CoR)))'; 
    position_out.x(ff_hits|rf_hits)=xf_r(ff_hits|rf_hits);
    position_out.z(ff_hits|rf_hits)=zf_r(ff_hits|rf_hits);
    
    
    % find left-crest hit position (first collision)
    position_out.x(lc_hits) = xc_l(lc_hits);
    position_out.z(lc_hits) = zc_l(lc_hits);
    
    % scatter from left-crest
    ang_p1 = -atan(position_out.x(lc_hits)./position_out.z(lc_hits));
    R1l=[permute([cos(ang_p1) zeros(numel(ang_p1),1) sin(ang_p1)],[3 2 1]);repmat([0 1 0],[1 1 numel(ang_p1)]);permute([-sin(ang_p1) zeros(numel(ang_p1),1) cos(ang_p1)],[3 2 1])]; % handle this via managing the single hits in 2D
    ii1cl=find(lc_hits);
    vp1sl=nan(sum(lc_hits),3);
    for ii=1:size(vp1sl,1)
        vp1sl(ii,:) = (R1l(:,:,ii)^-1*((pars.Px+pars.Py-pars.Pz).*harmolen_CoR(R1l(:,:,ii)*velocity_in(ii1cl(ii),:)',pars.CoR,ang_p1))*R1l(:,:,ii)*velocity_in(ii1cl(ii),:)')'; 
    end
    velocity_out(lc_hits,:)=vp1sl;
    
    % find right-crest hit position (first collision)
    position_out.x(rc_hits) = xc_r(rc_hits);
    position_out.z(rc_hits) = zc_r(rc_hits);
    
    % scatter from right-crest
    ang_p1 = -atan((2*pars.gCos-position_out.x(rc_hits))./position_out.z(rc_hits));
    R1r=[permute([cos(ang_p1) zeros(numel(ang_p1),1) sin(ang_p1)],[3 2 1]);repmat([0 1 0],[1 1 numel(ang_p1)]);permute([-sin(ang_p1) zeros(numel(ang_p1),1) cos(ang_p1)],[3 2 1])]; % handle this via managing the single hits in 2D
    ii1cr=find(rc_hits);
    vp1sr=nan(sum(rc_hits),3);
    for ii=1:size(vp1sr,1)
        vp1sr(ii,:) = (R1r(:,:,ii)*((pars.Px+pars.Py-pars.Pz).*harmolen_CoR(R1r(:,:,ii)^-1*velocity_in(ii1cr(ii),:)',pars.CoR,ang_p1))*R1r(:,:,ii)^-1*velocity_in(ii1cr(ii),:)')'; 
    end
    velocity_out(rc_hits,:)=vp1sr;
end
