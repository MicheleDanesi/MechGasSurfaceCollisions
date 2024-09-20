%% HARvesting MOLecular ENergy via corrugation-induced asymmetry of momentum transfer 
%  A simple geometrical model

clear all

% PARAMETERS: ------------------------------------

% geometry
    groove_angle_tangent = 1;     % groove angle's tangent 

% coefficient of restitution 

    s = 10;                      % scaling constant; 0 -> no loss
    k_nf = 0.02;                 % normal on face
    k_tf = 0.01;                 % tangential on face
    k_nc = 0.03;                 % normal on crest
    k_txc = 0.04;                % transverse tangential on crest
    k_tyc = 0.01;                % longitudinal tangential on crest

%   kinematics
    u= -0.01;                     % plate velocity 

% simulation numbers    
    nparts = 1e6;                 % total number of particles
    rp = 10;                      % particle radius (starting value); units of the groove side
    rp_min = 0.0001;              % particle radius (stop value)

    nsteps=1000;                    % number of simulations (varying radius between rp and rp_min)

% display    
    nplots = NaN;                   % number of hits to display (NaN for no plots; 0 for all, but beware this is nparts...)
    plot_pause = 0.01;             % display time

% -------------------------------------------------

% pack parameters
step_ratio=(rp/rp_min)^(1/nsteps);
groove_angle = atan(groove_angle_tangent);

pars.r=rp;
pars.u=u;
pars.gTan = groove_angle_tangent;
pars.gAng = atan(groove_angle_tangent);
pars.gCos = cos(pars.gAng);
pars.gSin = sin(pars.gAng);

pars.CoR.k_nf = s*k_nf; % normal on face
pars.CoR.k_tf = s*k_tf; % tangential on face
pars.CoR.k_nc = s*k_nc; % normal on crest
pars.CoR.k_txc = s*k_txc; % transverse tangential on crest
pars.CoR.k_tyc = s*k_tyc; % longitudinal tangential on crest

pars.walls.left.a = groove_angle_tangent; % left wall   (a*x+b*z+c=0)
pars.walls.right.a = -groove_angle_tangent;
pars.walls.left.b = -1;
pars.walls.right.b = -1;
pars.walls.left.c = 0;
pars.walls.right.c = 2*pars.gSin; %right wall  (a*x+b*z+c=0)

% common projections and rotations

pars.Px=[1 0 0;0 0 0;0 0 0];
pars.Py=[0 0 0;0 1 0;0 0 0];
pars.Pz=[0 0 0;0 0 0;0 0 1];
pars.Rl=[cos(groove_angle) 0 sin(groove_angle); 0 1 0; -sin(groove_angle) 0 cos(groove_angle)];
pars.Rr=[cos(-groove_angle) 0 sin(-groove_angle); 0 1 0; -sin(-groove_angle) 0 cos(-groove_angle)];

% main

for jj=1:nsteps  % number of simulations (scaling radius at each step)
    
    tic
    
    in.v = randn(nparts,3); % generate random velocity distribution
    in.v(:,3)=(in.v(:,3)-u); % adopt the moving plate frame of reference
    
    out.v = nan(size(in.v)); % prepare empty output 
    out.v(in.v(:,3)<=0,:)=in.v(in.v(:,3)<=0,:);
    out.v(in.v(:,3)<=0,:)=out.v(in.v(:,3)<=0,:).*harmolen_CoR(out.v(in.v(:,3)<=0,:)',pars.CoR)'; % downgoing particles hit the flat face
    
    out.n = nan(size(in.v,1),1); % prepare empty output for tracking the number of hits 
    out.n(in.v(:,3)<=0) = 0; 

    in.p.x = (rand(size(in.v,1),1)-0.5)*2*pars.gCos; % generate random hit position
    in.p.z = -ones(size(in.v,1),1)*pars.r;
    out.p.x = nan(size(in.v,1),1);
    out.p.x(in.v(:,3)<=0)=in.p.x(in.v(:,3)<=0);
    out.p.z = nan(size(in.v,1),1);
    out.p.z(in.v(:,3)<=0) = -pars.r;
    
    % use symmetry to limit handled events within a single groove, i.e. hande only positive positions
    ip = find(in.v(:,3)>0);
    sp = sign(in.p.x(ip));  % save sign of particle start position
    p.x = in.p.x(ip).*sp;   % mirror velocity
    p.z=zeros(size(p.x,1),1)-pars.r;
    v = [in.v(ip,1).*sp in.v(ip,2:3)];

    c=0;
    
    to_go = numel(ip);
    wht=zeros(numel(ip),1);
    hits = [];
    while to_go
        [p,v_out,wht]=harmolen_next_hit(p,v,wht,pars);
        no_hits = find(wht==0);
        if ~isempty(no_hits)
            vo = v_out(no_hits,:);
            xo = p.x(no_hits)-(pars.r+p.z(no_hits)).*v(no_hits,1)./v(no_hits,3);
            xo = xo.*sp(no_hits);
            vo = [vo(:,1).*sp(no_hits) vo(:,2:3)];

            out.p.x(ip(no_hits)) = xo;
            out.p.z(ip(no_hits)) = -pars.r;
            
            out.v(ip(no_hits),:) = vo;
            out.n(ip(no_hits),:) = c;

            v_out(no_hits,:) = [];
            p.x(no_hits) = [];
            p.z(no_hits) = [];
            ip(no_hits) = [];
            sp(no_hits) = [];
            wht(no_hits) = [];
            %numel(no_hits)
        end
        v = v_out;
        c=c+1;
        fhits(c).p.x = p.x.*sp;
        hits(c).p.z = p.z;
        hits(c).v = [v(:,1).*sp v(:,2:3)];
        hits(c).i = ip;
        hits(c).w = wht;
        to_go = numel(ip);
    end
    
    toc
    particle_radius = pars.r
    disp(['step ' num2str(jj) '/' num2str(nsteps)])
    
    rrp(jj) = pars.r;
    dvmg(jj) = mean((out.v(out.n>0,3)-in.v(out.n>0,3)).*in.v(out.n>0,3)); % pressure on the grooves side
    dvmf(jj) = mean((out.v(out.n==0,3)+in.v(out.n==0,3)).*in.v(out.n==0,3)); % pressure on the flat side  
    numhts(jj) = numel(hits)-1; % maximum number of hits
    nnht = []; 
    for nht=1:numel(hits)
        nnht(nht) = numel(hits(nht).i); % number of hits at hit index
        nnhtc(nht) = sum(hits(nht).w<4); % number of crest hits at hit index
    end
    meanhts(jj) = sum([1:numel(hits)].*nnht)/sum(nnht);  % mean number of hits
    snglhts(jj) = (nnht(1)-nnht(2))/nnht(1); % fraction of single hits 
    csthts(jj) = sum(nnhtc)/sum(nnht);

    % display a few hits
    if ~isnan(nplots)
        selected=find(in.v(:,3)>0);
        harmolen_plot_hits(selected(1:nplots),plot_pause,pars,in,out,hits) %groove_angle,rp,xp,zp,vp,xp1,zp1,vp1,xp2,zp2,vp2,xo,zo,vo,hts,hts2)
    end
    
    % display velocity distributions
    
    harmolen_plot_distributions(-in.v(out.n==0,:),out.v(out.n>0,:))
    
    % scale particle radius before next step
    pars.r=pars.r/step_ratio;

end

% display summary
figure(3)
semilogx(rrp',[-dvmg' dvmf']/2)
ylabel('pressure')
xlabel('particle size')

figure(4)
subplot(2,1,1)
semilogx(rrp',[numhts' meanhts'])
ylabel('collisions #')
xlabel('particle size')
legend({'maximum','mean'})
subplot(2,1,2)
semilogx(rrp',[snglhts' csthts'])
ylabel('collisions fraction')
xlabel('particle size')
legend({'single hits','crest hits'})


