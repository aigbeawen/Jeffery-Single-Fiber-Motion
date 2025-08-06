[Optfunc, Bfunc, DU, Rot]=Main_Func;
%
stl=6               ; % flowtype flag; See function DU for flags
GE=1                ; % shear-to-ext.;
G =1.               ; % G-shear rate ;
E =G/GE             ; % Ext. rate    ;
dU=DU(stl,[E,G])    ; % Vel Grad.    ;
%
ar=6.               ; % ar - aspect ratio;
tolc=1.e-6; c=1-tolc; re=[ar*(1+tolc),1+tolc]; 
Xf=[re(1);re(2);1]*c; % Xf=[a,b,c] - ellipsoid axis lengths 
Rinf=50*max(Xf)     ; % flow domain size
Xtip=[Xf(1); 0; 0]  ; % ellipsoid major axis tip location
p0=0                ; % Ref. pressure;
mu=1.               ; % Viscosity    ;
%
invar={[],Xf,Rinf,mu,p0,dU}; % Bfunc Input variables
%
% dUm=1/2*(dU+dU'); [Q,eigv]=eig(dUm); 
% for ak=1:3; Qk(:,ak)=acos(dot(Q(:,ak*[1,1,1]),eye(3),2)/norm(Q(:,ak))); end
%  --------------------- Optimization -------------------------------------
X0=[5.85; 0.15; 0.15]; Qt0=[0.;-pi/2; 0.]; QX0=[Qt0;X0]; 
QXa=[0;-pi/2;0;0;0;0]; QXb=[2*pi;pi/2;2*pi;6.;1.;1.];
[QXopt,Popt]=Optfunc(QX0,1, invar{:},QXa,QXb)
% ----------------------------ODE Solver-----------------------------------
t0=0.; dt=.5; tn=2*pi/G*(re(1)+1/re(1))+2.5; tp=t0:dt:1; 
Qs=pi/2; Q0=[0;-Qs; 0.]; % Initial Orientation
[tp, Qt]=ode23s(@(t,Q) Bfunc(t,[Q;Xtip],invar{:}),tp,Q0);
nt=length(tp);
for jt=1:nt
    [p(jt,:), V(jt,:), dQt(jt,:)]=Bfunc([],[Qt(jt,:)';Xtip],invar{:},1);
end
% -------------------------------------------------------------------------
dirname='.\Plots_3D'; if ~exist(dirname,'dir'), mkdir(dirname); end
npts=200;for j=1:3; pts{j}=linspace(-Xf(j)-c,Xf(j)+c,npts); end
[X,Y,Z]=meshgrid(pts{:}); Vt=X.^2/Xf(1)^2+Y.^2/Xf(2)^2+Z.^2/Xf(3)^2-1;
[els, nds]=isosurface(X,Y,Z,Vt,1e-4); nx=size(nds,1);
fid=fopen(sprintf('Results_0_%g.txt',stl),'w+');
% vid=fopen(sprintf('Vdata_3D_pi_%g.txt',Qs),'w+');
fprintf(fid,'t\tdphi\tdtta\tdgma\tphi\ttta\tgma\tPtip\tPmin\tPmax\n');
titles={'Vel Mag','Pres'};
v1 = VideoWriter('VelM_3D.avi'); v2 = VideoWriter('Pres_3D.avi');
v={v1,v2}; open(v1) ; open(v2) ;
for jt=1:nt
    [Rtj,Rtwj]=Rot(Qt(jt,:)); Xj=nds*Rtj';
    x=Xj(:,1); x=x(els); y=Xj(:,2); y=y(els); z=Xj(:,3); z=z(els); 
    for jx=1:nx
        [Pj(jx,1),Uj1(jx,:)]=Bfunc([],[Qt(jt,:),nds(jx,:)]',invar{:},1);
        % Uj2(jx,:)=cross(Rtj'*Rtwj*dQt(jt,:)',nds(jx,:)');
        % Vj(jx,:)=cross(Rtwj*dQt(jt,:)',Xj(jx,:)');
    end
    Vj=vecnorm(Rtj*Uj1')'; % err(jt)=norm(Uj2-Uj1);
    % N:B- Velocity Magnitude is independent on coordinate system
    [mpj,npj]=max([-1,1].*Pj,[],1); % xpj=nds(npj,:)';disp(xpj(:)')
    hstr=[tp(jt),dQt(jt,:),Qt(jt,:),p(jt),[-1,1].*mpj];
    fprintf(fid,['%.3f', repmat(',%10.4f',1,9),'\n'],hstr');
    %
    % Vtip=norm(V(jt,:)); [mvj,nvj]=max(Vj,[],1); pvj=Vj(npj(1),:);
    % vstr=[tp(jt),Vtip,mvj,nds(nvj,:),pvj,nds(npj(1),:)]';
    % fprintf(vid,['%.2f', repmat(',%10.4f',1,9),'\n'],vstr');
    for jf=1:2
        figj=figure(jf);  clf; figj.Color='w';
        if jf==1, fOut=Vj; else, fOut=Pj; end; Vt=fOut(els);
        img=patch(x',y',z',Vt'); set(img,'EdgeAlpha',0); 
        set(gca,'Colormap',jet); set(img,'FaceLighting',"gouraud");
        axis("off"); view(60,30); axis("equal");
        title(['\rm\it',titles{jf}]); cb=colorbar;
        saveas(figj,sprintf('%s\\%s_%g',dirname,titles{jf},jt),'png');
        pause(0.2); F=getframe(figj);  writeVideo(v{jf},F);
    end
end
fclose(fid);close(v1); close(v2); %fclose(vid); 
%

