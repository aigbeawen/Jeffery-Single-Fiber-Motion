function varargout=Jfunc
varargout={@Optfunc, @Bfunc,@DU};
end
%
function [QXopt,Popt]=Optfunc(QX0,flg, varargin)
Xf=varargin{2}; I=eye(3); O=zeros(3);
kq=[0;1;0]; vq=[0;pi/2;0]; 
switch flg
    case 1 % fmincon
        [QXa,QXb]=varargin{7:8};
        opts=optimoptions('fmincon','Algorithm','interior-point',...
            'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
            'HessianFcn', @hessfcn);
        [QXopt,Popt]=fmincon(@Obj,QX0,[],[],[],[],QXa,QXb,@Con,opts);
    case 2 % my_quadprog-
        % minimize f=c'*d+1/2*d'*H*d; S.t. Aeq'*d+beq=0
        x=QX0; err=1; dk=zeros(7,1); r=1.e+08;
        while err>1e-6
            d=dk(1:end-1); x=x+d; v=dk(end);
            [~,c,H] = Obj(x); [~,beq,~,Aeq] = Con(x);
            f=c'*d+1/2*d'*H*d; df=c+H*d; ceq=Aeq'*d+beq;
            R=[df+Aeq*v+1/2*r*Aeq*ceq;ceq]; J=[H+1/2*r*(Aeq*Aeq'),Aeq;Aeq',0];
            dk=dk-2.*J\R; err=norm(R)
        end
        QXopt=x;Popt=Obj(x);
    case 3 % MLB_quadprog
        opts=optimoptions('quadprog' ,'Display','off');
        x=QX0; err=1;
        while err > 1.e-6
            [~, beq, ~, Aeq]=Con(x); [f,c,H] = Obj(x);
            d=quadprog(H,c,[],[],Aeq',-beq,[],[],[],opts);
            err=norm(d), x = x + d;
        end
        QXopt=x;Popt=Obj(x);
end
% ------------------------ Objective & Constraint Funcs.-------------------
    function [p,Dp,Hp]=Obj(X)
        [p,Dp,Hp]=Bfunc([],X,varargin{1:6},2);
    end
%
    function [c,ceq,dc,dceq]=Con(X)
        Qc=X(1:3); Xc=X(4:6); c=[]; dc=[];
        ceqX=(1./Xf.^2)'*Xc.^2-1; dceqX=[0;0;0;2*Xc./Xf.^2];
        ceqQ=kq.*Qc+vq; dceqQ=[kq.*I; O];
        ceq=[ceqQ;ceqX]; dceq=[dceqQ,dceqX];
    end
%
    function H=hessfcn(X,lambda)
        [~,~,Hf]=Bfunc([],X,varargin{1:6},2);
        Heqx=[O,O;O,2*(1./Xf.^2).*I];
        H = Hf + lambda.eqnonlin(4)*Heqx;
    end
%
    function [f,df,Hf]=fd_Obj(X)
        del=1.e-4; flg=3; % del - fin. dif. pert. size and type
        f=Bfunc([],X,varargin{1:6});
        [Hf,df]=dMdA(@(X) dMdA(@(X) Bfunc([],X,varargin{1:6}), X, del, flg), X, del, flg);
        Hf=reshape(Hf,6,6); df=reshape(df,6,1);
    end
end
%
% --------------------------- Main Function -------------------------------
function [p,Dp,Hp]=Bfunc(tx,QX,dQt,Xf,Rinf,mu,p0,dU,flg)
Qt=QX((1:3)'); X=QX((4:6)');
n1=[1,5,9]; m1=[8,7,4]; n23=[2,3;3,1;1,2]; Xf32=Xf(n23(:,[2,1]));
I=eye(3); nI=1-I; [Rt, Rtw]=Rot(Qt); 
%
eijk=zeros(3,3,3);
for k1=1:3
    for k2=1:3
        for k3=1:3
            eijk(k1,k2,k3)=(k1-k2)*(k2-k3)*(k3-k1)/2;
        end
    end
end
%
DU=Rt'*dU*Rt; DU=DU-1/3*trace(DU)*I;
%
dUm=1/2*(DU+DU'); dQm=1/2*(DU-DU');
abc=dUm(n1')    ; fgh=dUm(m1')    ;
dQv=dQm([6;7;2]) ; 
% -------------------------------------------------------------------------
if isempty(dQt)
    Qflg=1;
    nj=1:3; pj=nj-nj'; mI=pj.*(1.-.5*(pj.^2-1));
    dQ=dQv+((mI*Xf.^2)./(nI*Xf.^2)).*fgh;  dQt=Rtw\Rt*dQ;
    
else
    Qflg=2;
    dQ=Rt'*Rtw*dQt; 
end

% -------------------------------------------------------------------------
if ~isempty(tx)
    p=dQt;
else
    R=@(v)  (X'.^2)*(1./(Xf.^2+v))-1 ;
    J=@(v) -(X'.^2)*(1./(Xf.^2+v).^2);
    v=0; err=1;
    while err>1e-7
        v=v-R(v)/J(v);  err=abs(R(v));
    end
    %
    K=1./(Xf.^2+v);
    [lbg0_0, lbg1_0, lbg2_0]=lbg(0, Rinf);
    [lbg0_v, lbg1_v, lbg2_v]=lbg(v, Rinf);
    Vx=integral(@(v) R(v)./D(v),v,inf);
    GVx = 2*lbg0_v.*X; HVx=2*lbg0_v.*I-1/(P(2)*D(v))*(Gv()*Gv()');
    %
    ABC_2=sum(prod(lbg2_0(n23),2)); ABC=1/6*(3*I-1)*(lbg2_0.*abc)/ABC_2;
    FGH_2=(2*lbg1_0.*nI*(lbg0_0.*Xf.^2));
    FGH=(lbg0_0(n23).*fgh+[-1,1].*Xf32.^2.*lbg1_0.*(dQv-dQ))./FGH_2;
    sFGH=sum(FGH,2);
    %
    p =p0+2*mu*(HVx(n1)*ABC+ HVx(m1)*sFGH);
    %
    if nargout>1
        switch flg
            case 1
                [Xi,dXi]=fXi();
                RST=-fgh./lbg1_0; UVW=2*sum(eijk,3)*(ABC.*Xf.^2);
                vXi=dXi'*RST;
                for k1=1:3
                    for k2=1:3
                        for k3=1:3
                            vXi(k1)=vXi(k1)-eijk(k1,k2,k3)*UVW(k2)*dXi(k2,k3);
                        end
                    end
                end

                Cof=[ABC(1),FGH(3,1),FGH(2,2) ;
                    FGH(3,2),ABC(2),FGH(1,1) ;
                    FGH(2,1),FGH(1,2),ABC(3)];

                Vel=DU*X+vXi+HVx*Cof'*X-Cof*GVx;
                Velf=cross(dQ,X); % if (X'.^2)*(Xf.^-2)-1=0
                %
                Dp=Vel; if nargout>2, Hp=dQt; end
            case 2
                [~, ~, GRt,GRtw]=Rot(Qt);
                for j1=1:3
                    GRtj=GRt(:,:,j1);
                    GDUj=GRtj'*dU*Rt+Rt'*dU*GRtj;
                    GDUj=GDUj-1/3*trace(GDUj)*I;
                    %
                    GdUmj=1/2*(GDUj+GDUj'); GdQmj=1/2*(GDUj-GDUj');
                    Gabc(:,j1)=GdUmj(n1)  ; Gfgh(:,j1)=GdUmj(m1)  ;
                    GdQv(:,j1)=GdQmj([6;7;2]);
                    if nargout>2
                        for j2=1:3
                            [~, ~, ~, ~, HRt,HRtw]=Rot(Qt);
                            GRtk=GRt(:,:,j2); HRtjk=HRt(:,:,j1,j2);
                            HDUjk=HRtjk'*dU*Rt+GRtj'*dU*GRtk+GRtk'*dU*GRtj+Rt'*dU*HRtjk;
                            HDUjk=HDUjk-1/3*trace(HDUjk)*I;
                            %
                            HdUmjk=1/2*(HDUjk+HDUjk'); HdQmjk=1/2*(HDUjk-HDUjk');
                            Habc(:,j1,j2)=HdUmjk(n1)  ; Hfgh(:,j1,j2)=HdUmjk(m1)  ;
                            HdQv(:,j1,j2)=HdQmjk([6;7;2]);
                        end
                    end
                end
                %
                switch Qflg
                    case 1
                        GdQ=GdQv+((mI*Xf.^2)./(nI*Xf.^2)).*Gfgh;
                        if nargout>2
                            HdQ=HdQv+((mI*Xf.^2)./(nI*Xf.^2)).*Hfgh;
                        end
                    case 2
                        GdQ=Rt'*tensorprod(GRtw,dQt,2,1)+...
                                tensorprod(GRt,Rtw*dQt,1,1);
                        if nargout>2
                            HdQ=tensorprod(Rt,tensorprod(HRtw,dQt,2,1),1,1)+...
                                2*permute(tensorprod(GRt,tensorprod(GRtw,dQt,2,1),1,1),[1,3,2])+...
                                tensorprod(HRt,Rtw*dQt,1,1);
                        end
                end
                %
                GHVx=GHV(); if nargout>2, HHVx=HHV(); end
                for j1=1:3
                    GABC(:,j1)=1/6*(3*I-1)*(lbg2_0.*Gabc(:,j1))/ABC_2;
                    GFGH(:,j1)=sum((lbg0_0(n23).*Gfgh(:,j1)+...
                        [-1,1].*Xf32.^2.*lbg1_0.*(GdQv(:,j1)-GdQ(:,j1)))./FGH_2,2);
                    if nargout>2
                        for j2=1:3
                            HABC(:,j1,j2)=1/6*(3*I-1)*(lbg2_0.*Habc(:,j1,j2))/ABC_2;
                            HFGH(:,j1,j2)=sum((lbg0_0(n23).*Hfgh(:,j1,j2)+...
                                [-1,1].*Xf32.^2.*lbg1_0.*(HdQv(:,j1,j2)-HdQ(:,j1,j2)))./FGH_2,2);
                        end
                    end
                end
                %
                n2=n1+3^2*(0:2)'; m2=m1+3^2*(0:2)';
                dpdX=mu*(GHVx(n2)*ABC+GHVx(m2)*sFGH);
                dpdQ =mu*(HVx(n1)*GABC+ HVx(m1)*GFGH)';
                Dp =2*[dpdQ;dpdX];
                if nargout>2
                    n3=n2+3^3*permute(0:2,[1,3,2]); m3=m2+3^3*permute(0:2,[1,3,2]);
                    HpHX=mu*(sum(HHVx(n3).*ABC',2)+sum(HHVx(m3).*sFGH',2));
                    HpHX=permute(HpHX,[1,3,2]);
                    HpHQ=mu*(tensorprod(HVx(n1),HABC,2,1)+ tensorprod(HVx(m1),HFGH,2,1));
                    HpHQ=reshape(HpHQ,3,3);
                    HpHQX=mu*(GHVx(n2)*GABC+GHVx(m2)*GFGH);
                    Hp=2*[HpHQ,HpHQX';HpHQX,HpHX];
                end
        end
    end
end
% -------------------------------------------------------------------------
    function f=dP(n)
        f=-P(2)*(K.^(n-1)-n*P(2)/P(n+1)).*Gv;
    end
% Fourth Derivative - Omega
    function f=HHV
        dim=1:3; G=Gv(); H=Hv(); GH=GHV();
        for i=dim
            for j=dim
                for k=dim
                    for l=dim
                        ijkl=[i,j,k,l]; [frq,ord]=sort(histc(ijkl,dim),'descend');
                        n=sum(frq.^2); m=n/2-2; m=m-2*(m-1)*(m-2)*(m-3)/(5*4*3);
                        rst=dim(ord); r=rst(1);s=rst(2);t=rst(3);
                        G4=prod(G(rst).^frq(:));
                        switch m 
                            case 1
                                fn=(1/X(r)-(3/2*Y(1)+2*(K(r)-2*P(2)/P(3)))*G(r))*GH(r,s,t)+1/(P(2)*D(v))*...
                                    (-3/2*Y(2)-2*K(r)^2+6*P(2)/P(4)+2*P(2)/P(3)*(K(r)-2*P(2)/P(3)))*G4;
                            case 2
                                fn=(-Y(1)/2*G(s)+2/G(r)*H(r,s)+1/G(s)*H(s,s))*GH(r,r,s)+...
                                    2/D(v)*K(r)*G(s)*(K(r)*G(s)+2/G(r)*H(r,s))+...
                                    1/(P(2)*D(v))*(-(Y(2)/2+2*K(r)^2+K(s)^2)+(Y(1)/2+2*K(r)+K(s))*(K(s)-2*P(2)/P(3))-...
                                    2*(K(s)^2-3*P(2)/P(4)))*G4;
                            case {3,4}
                                switch m
                                    case 4
                                        s=r;
                                end
                                fn=(-Y(1)/2*G(s)+3/G(r)*H(r,s))*GH(r,r,r)+6/D(v)*K(r)*(K(r)*G(r)*G(s)+2*H(r,s))+...
                                    1/(P(2)*D(v))*(-(Y(2)/2+3*K(r)^2)+(Y(1)/2+3*K(r))*(K(s)-2*P(2)/P(3))-...
                                    2*(K(s)^2-3*P(2)/P(4)))*G4;
                        end
                        f(i,j,k,l)=fn;
                    end
                end
            end
        end
    end
% Third Derivative - Omega
    function f=GHV
        dijk=2*permute(I,[1,3,2]).*I+I;
        %
        G=Gv();
        for i =1:3
            f=G(i)*K.*dijk(:,:,i);
            for j=1:3
                if i~=j
                    f=f+K(i)*G(6-i-j)*abs(eijk(:,:,j));
                end
            end
            GHv1(:,:,i)=f;
        end
        GHv2=-1/(2*P(2))*(Y(1)/2+((K+K')+...
            permute(K,[3,2,1]))-2*P(2)/P(3)).*...
            ((G*G').*permute(G,[3,2,1]));
        f=-2/D(v)*(GHv1+GHv2);
    end
% f1=[gradlbg0_v,gradbg1_v,gradbg2_v]; f2=[hessbg0_v::hesslbg1_v::hesslbg2_v];
    function [f1,f2]=Dlbg
        bnI=logical(nI); G1=Glbg(); H1=Hlbg(); 
        G2=zeros(3); G3=zeros(3); H2=zeros(3,3,3); H3=zeros(3,3,3);
        for j=1:3
            bol=bnI(j,:); Xfj=Xf(bol).^2; Gj=G1(:,bol); Hj=H1(:,:,bol);
            dXfj=diff(Xfj);
            G2(:,j)=diff(Gj,1,2)/dXfj; G3(:,j)=diff(Xfj'.*Gj,1,2)/dXfj; 
            H2(:,:,j)=diff(Hj,1,3)/dXfj;
            H3(:,:,j)=diff(permute(Xfj,[3,2,1]).*Hj,1,3)/dXfj; 
        end
        f1=cat(3,G1,-G2,G3); f2=cat(4,H1,-H2,H3); 
    end
%
    function f=Hlbg
        K3=permute(K,[3,2,1]);
        f=K3/D(v).*((Y(1)/2+((K+K')+K3)-...
            2*P(2)/P(3)).*(Gv*Gv')-Gv./X.*I);
    end
%   
    function f=Glbg % [Glfa,Gbta,Ggma]
        f=-K'/D(v).*Gv;
    end
% Hessian - Lambda
    function f=Hv
        I=eye(3);
        f=2*P(2)*K.*I-((K+K')-2*P(2)/P(3)).*(Gv.*Gv');
    end
% Gradient - Lambda
    function f=Gv
        f=2*X.*K*P(2);
    end
%
    function [f,df]=fXi
        XY=prod(I+X'.*nI,2); ab=prod(I+(Xf'.^2+v).*nI,2);
        df1=zeros(3);
        for r1=1:3
            for r2=1:3
                if r2~=r1
                    df1(r1,r2)=X(6-r1-r2);
                end
            end
        end
        df1=lbg1_v.*df1;
        df2=-2*P(2)/D(v)*(XY*X')./(ab*(Xf'.^2+v));
        f=lbg1_v.*XY; df=df1+df2;
    end
%
    function f=P(n)
        f=(X'.^2)*(K.^n); f=1/f;
    end
%
    function f=Y(n)
        f=[1,1,1]*(K.^n);
    end
%
    function varargout=lbg(varargin)
        [v0,v1]=varargin{:};
        for j=1:3
            f0(j,1)=integral(@(v) 1./((Xf(j)^2+v).*D(v)),v0,v1);
        end
        M=[Xf.^2,f0,f0.*Xf.^2]; M=diff(M([1:3,1],:),1,1);
        f1=-M(:,2)./M(:,1); f2=M(:,3)./M(:,1);
        f1=f1([2,3,1]); f2=f2([2,3,1]);
        varargout={f0, f1, f2};
    end
%
    function f=D(v)
        f=sqrt(prod(Xf.^2+v));
    end
end
% -------------------------------------------------------------------------
function [Rt,Rtw,GRt,GRtw,HRt,HRtw]=Rot(Q)
cQ=cos(Q); sQ=sin(Q);
% -------------------------------------------------------------------------
Rx1=[1     0     0;  0     cQ(1)  sQ(1);0 -sQ(1) cQ(1)];
Rz =[cQ(2) sQ(2) 0; -sQ(2) cQ(2)  0    ;0  0     1    ];
Rx2=[1     0     0;  0     cQ(3)  sQ(3);0 -sQ(3) cQ(3)];
Rt=Rx1'*Rz'*Rx2';
%
Rtw=[1   0      cQ(2)        ;
     0  -sQ(1)  sQ(2)*cQ(1)	 ;
     0   cQ(1)  sQ(2)*sQ(1)	];
%
% Rtw2=[cQ(2),-sQ(2)*cQ(3),sQ(2)*sQ(3);
%           0,       sQ(3),      cQ(3);
%           1,           0,         0];
% -------------------------------------------------------------------------
if nargout>2
GRx1=[0      0     0;  0     -sQ(1)  cQ(1);0 -cQ(1) -sQ(1)];
GRz =[-sQ(2) cQ(2) 0; -cQ(2) -sQ(2)  0    ;0  0      0    ];
GRx2=[0      0     0;  0     -sQ(3)  cQ(3);0 -cQ(3) -sQ(3)];
GRt=cat(3,GRx1'*Rz'*Rx2',Rx1'*GRz'*Rx2',Rx1'*Rz'*GRx2');
%
GRtw=zeros(3,3,3);
GRtw(:,:,1)=[0, 0,     0 ; 0,-cQ(1),-sQ(2)*sQ(1);0,-sQ(1),sQ(2)*cQ(1)];
GRtw(:,:,2)=[0, 0, -sQ(2); 0,    0 , cQ(2)*cQ(1);0,     0,cQ(2)*sQ(1)];
% ------------------------------------------------------------------------
if nargout>4
HRx1=[0       0     0;  0     -cQ(1) -sQ(1);0  sQ(1) -cQ(1)];
HRz =[-cQ(2) -sQ(2) 0;  sQ(2) -cQ(2)  0    ;0  0      0    ];
HRx2=[0       0     0;  0     -cQ(3) -sQ(3);0  sQ(3) -cQ(3)];
HRt(:,:,:,1)=cat(3,HRx1'*Rz'*Rx2' ,GRx1'*GRz'*Rx2',GRx1'*Rz'*GRx2');
HRt(:,:,:,2)=cat(3,GRx1'*GRz'*Rx2',Rx1'*HRz'*Rx2' ,Rx1'*GRz'*GRx2');
HRt(:,:,:,3)=cat(3,GRx1'*Rz'*GRx2',Rx1'*GRz'*GRx2',Rx1'*Rz'*HRx2');
%
HRtw=zeros(3,3,3,3);
HRtw(:,:,1,1)=[0, 0,     0 ; 0, sQ(1), -sQ(2)*cQ(1);0,-cQ(1),-sQ(2)*sQ(1)];
HRtw(:,:,2,1)=[0, 0,     0 ; 0,     0, -cQ(2)*sQ(1);0,     0, cQ(2)*cQ(1)];
HRtw(:,:,1,2)=[0, 0,     0 ; 0,     0, -cQ(2)*sQ(1);0,     0, cQ(2)*cQ(1)];
HRtw(:,:,2,2)=[0, 0, -cQ(2); 0,     0, -sQ(2)*cQ(1);0,     0,-sQ(2)*sQ(1)];
end
end
end
% -------------------------------------------------------------------------
function [J,R]=dMdA(func, A, del, flg, varargin)
A=A(:); n2=numel(A); nA=size(A); nA=num2cell(nA);
D=zeros(n2,1); tol=eps;
f=@(A) func(A,varargin{:}); R=f(A); n1=numel(R);
%
J=zeros(n1*n2,1);
Fac={ 
    1,         [-1  1 ; -1  0] ;
    1,         [-1  1 ;  0  1] ;
    2,         [-1  1 ; -1  1] ;
    2,     [1 -4  3 ; -2 -1 0] ;
    2,     [-3 4 -1 ;  0  1 2] ;
   12, [1 -8  8 -1 ; -2 -1 1 2]
                                    };
%
M1=Fac{flg,1}; M2=Fac{flg,2}; m=length(M2); 
for i=1:n2
    na=(i-1)*n1+1; nb=na+n1-1;
    D(i,1)=1; J(na:nb,1)=zeros(n1,1); dA=del*D; h=del; 
    if abs(A(i))>tol,  dA=dA.*A(:); h=h*A(i);  end
    for j=1:m
        Jj=M2(1,j)*f(A+M2(2,j)*reshape(dA,nA{:})); 
        J(na:nb,1)=J(na:nb,1)+Jj(:);
    end
    J(na:nb,1)=J(na:nb,1)/(M1*h);
    D(i,1)=0;
end
nR=size(R); nR=num2cell(nR);J=reshape(J,nR{:},n2);
end
% -------------------------------------------------------------------------
function f=DU(flg,val)
f=zeros(3);
switch flg
    case 1 % Simple shear
        G=val(2); f(6)=G;
    case 2 % Shearing/stretching
        E=val(1); G=val(2); f([1,5,6,9])=[-E,-E, G, 2*E];
    case 3 % Uniaxial
        E=val(1); f([1,5,9])=[-E,-E,2*E];
    case 4 % Biaxial
        E=val(1); f([1,5,9])=[-2*E, E, E];
    case 5 % Shearing/planar stretching
        E=val(1);G=val(2); f([1, 6, 9])=[-E,G, E];
    case 6 % Balanced shear/biaxial elongation flow
        E=val(1);G=val(2); f([1,5,6,9])=[-2*E, E, G, E];
    case 7 % Triaxial
        E=val(1); f([1,5,9])=[E, E, E];
    case 8 % Balanced shear/triaxial elongation flow
        E=val(1);G=val(2); f([1,5, 6, 9])=[E, E, G, E];
end
end