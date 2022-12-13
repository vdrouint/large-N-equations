for asd4=1
    freshstart=1;
    if freshstart
        clear all 
        name='yk_mv_hund';
        name1=[name '.m'];
        name2=[name '_data'];
        
        options = optimset('Display','off');
        %---------------------------------------------------------------
        valence=0;
        %0 means Kondo
        %1 means infinite U Anderson model
        
        logridflag=1;
        becon=0.0; %necessary for 1D chain
        geom=0.5; 
        %geom = 0 is single impurity
        %geom = 0.1 is 2-site Hund impurity
        %geom = 0.5 is a Hund's impurity (3 sites)

        %------------------------------------------------------------
        if logridflag
            N=1000;
            D=500; ran = 2*D;
            minFreq = -4;
            maxFreq = 5;
        else
            %N=12000;
            N=4000;
            D=50; ran=2*D;
        end
        rho=1/D/2;
        howoften0=20;
        dispflag=1;
        disp('initializing the grid...');
        %---------------------------------------------------------
        %-------------------------------
        f=@(xoT) 1./(exp(xoT)+1);
        nB=@(xoT) 1./(exp(xoT)-1);
        dfomfun=@(xoT,T) xoT/(4*T)./(cosh(xoT/2).^2);
        d2fomfun=@(xoT,T) -1/(4*T)./(cosh(xoT/2).^2);
        dnBomfun=@(xoT,T) xoT/(4*T)./(sinh(xoT/2).^2);
        omcut=100.3243;
        cut=@(xoT) xoT.*(abs(xoT)<(omcut) & abs(xoT)>eps)+omcut*sign(xoT).*(abs(xoT)>=omcut)+eps*sign(xoT).*(abs(xoT)<eps);
        %-----------------------------
        %values of Tk/Delta
        %x0 = [0.1];
        %x0 = 10.^(0:-0.07:-1.4);
        x0 = 1;
        %x0 = [0.0035];

        
        %initial Q/N
        q00=0.3;
        %gamma = gammaFactor*q00;
        gammaFactor = 1.0;
        
        %use B as gammafactor (? sometimes)
        B0= [0.99];
        
        %initial hybridization
        %V20=1/(pi*rho);
        V20 = 25/rho;
        %epsf0= -0.3:0.05:-0.1;
        %epsf0 = [-3];
        epsf0 = [-0.100/rho]; 
        %epsmin = -0.2;
        epsmin = 10;
        %epsf0 = 0.2;
        
        maxphaseshift = 0.9995;

        %-----------------------------
        k0=(-1:.5:1)*pi; %k0=(-1:.05:1)*pi;
        TK = 1;
        T0i=50;
        Tmax = 1000;
        Tmin=1/20000000;
        lowTsweep=1;
        dec01=1.2;
        decEta = 1.2; %make deceta = dec01 in mv
        %dec02=1.15;
        %etaFactor = T0;
        tdown0 = 1;
        
        maxIterations1 = 501;
        maxIterations2 = 500;
        
        if lowTsweep
            Tmin2=Tmin;
        else
            Tmin2=Tmin;
        end
        %-------------------------------
        qq0=0;
        for B=B0
            for q0=q00
                for V2=V20
                    for epsf=epsf0
                        for x=x0
                            qq0=qq0+1;
                            Bq(qq0)=B;
                            q0q(qq0)=q0;
                            xq(qq0)=x;
                            V2q(qq0)=V2;
                            epsfq(qq0)=epsf;
                        end
                    end
                end
            end
        end
    end
end % initialization
%-----------------------------------
disp('computing gc...');
if logridflag
    %omp=logrid2(.1,TK*5,N,0.4,0.5);
    omp=10.^(minFreq:((maxFreq - minFreq)/(N-1)):maxFreq);    
else
    del=ran/N;
    omp=del:del:ran;
end
om=[-fliplr(omp) omp]; omp0=omp;
if logridflag
    dom=gradient(om);
    BmX=(om')*ones(1,length(om))-ones(length(om),1)*om;
else
    dom=om(2)-om(1);
end
izom=find((om)==min(abs(om))); lzom=length(izom);
%------------------------------
len=length(om);
IRind=[izom(1)+1-(100:-1:1) izom(end)-1+(1:100)];
UVind=[1:100 len+1-(100:-1:1)];
%------------------------------
eta=1e-2*D;
if 0                    % This is a 2D DoS for the c-electrons
    gc=0; gcBmX=0; gc2=0; gc2BmX=0;
    e0=-.995:.01:.995; de=e0(2)-e0(1);
    for ex=e0
        H=ex*D; dos=2/(pi^2)/D*ellipke(1-ex.^2);
        gc=gc+(D*de*dos)*1./(om+1i*eta-H);
        gc2=gc2-(D*de*dos)*1./((om+1i*eta-H).^2);
        if logridflag
            gcBmX=gcBmX+(D*de*dos)*1./(BmX+1i*eta-H);
            gc2BmX=gc2BmX-(D*de*dos)*1./((BmX+1i*eta-H).^2);
        end
        fprintf('.');
    end
    fprintf('\n');
else
    gcf=@(z) -rho*(log(D-z)-log(-D-z));
    gc2f=@(z) -rho*(1./(z-D)-1./(z+D));
    gc=gcf(om+1i*1e-3*D);
    if logridflag; gcBmX=gcf(BmX+1i*1e-3*D); end
    if 1 % charge
        gc2f=@(z) -rho*(1./(z-D)-1./(z+D));
        gc2=gc2f(om+1i*1e-3*D);
        if logridflag; gc2BmX=gc2f(BmX+1i*1e-3*D); end
    end
end
igc=imag(gc);
if logridflag
    gcmBmX=-conj(gcBmX);
    gcmXmB=transpose(gcBmX);
end
%iGp=-om.*exp(-abs(om)/10);
%Gp=yhilbert(iGp);
disp('done here!')
%-------------------------------
qq=0;
while qq<qq0
    qq=qq+1;
    B=Bq(qq);
    q0=q0q(qq);
    gamma=gammaFactor*q0;  %can use up to 1.03 in there
    x=xq(qq);
    Del=1/x;     % This is FM coupling
    
    if x > 5000
        Del = 0.0;
    end
        
    %-------------------
    %dec0 = dec01 + (0.2 - epsf)/4;
    dec0 = dec01;
    
    V2=V2q(qq);
    
    %modify T0
    if epsmin==10
        T0 = T0i;
    else     
        T0 = T0i*10^(1 - 5*(-epsmin - epsfq(qq)));
    end
    
    
    %epsf=epsfq(qq)/rho;
    epsf = epsfq(qq);
    J=1/rho/log(D/TK);
    etaFactor = T0;
    
    if ~valence
        epsf=-1/J;
        V2=1;
        etaFactor = 1e-0*T0;
    end
    
    
 
    err = etaFactor/V2;    

    %lambda initial for Kondo model
    lam=T0*log(1/q0+1);
    %more accurate for anderson model
    %funSolve1 = @(u) 1./(exp(u/T0) -1) + q0./(exp((u+epsf)/T0) +1) - q0;
    %lam = fsolve(funSolve1,T0*log(1/q0+1), options);
    

    disp(['< ' num2str(qq) ' =================='])
    disp(['New T0=' num2str(T0)]);
    disp(['epsf * rho=' num2str(epsf)]);
    disp(['q=' num2str(q0)]);
    disp(['k=' num2str(gamma)]);
    disp(['TK/Del=' num2str(x)]);
    disp(['Del=' num2str(Del)]);
    disp(['pi*rho*V^2=' num2str(rho*pi*V2)]);
    disp(['rho*J=' num2str(rho*J)]);
    if dispflag
        disp('---------------------')
    end
    tic
    %------------------------------
    %setting up the important parameters
    eta=err;
    howoften=howoften0;
    mixing=0.5;
    bec=0;
    T=T0; nBom=nB(cut(om/T)); fom=f(cut(om/T)); 
    etaB=yhilbert(1-2*fom); etaX=yhilbert(ones(size(om)));
    etaConst = (V2*5);
    if ~valence
        etaB = yhilbert(1-2*fom); etaX = zeros(size(om));
        etaConst = V2*100;
    end
    tgc=gc;

    if logridflag
        fBmX=f(cut(BmX/T)); fmXmB=transpose(fBmX);
        dfomBmX=dfomfun(cut(BmX/T),T);
    end
    %-------------

    %bosonic spinons 
    if geom == 0
        SigB=-eta.*etaB;
        GB=1./(om-lam-V2*SigB);
    end

    if geom == 0.1
        SigB=-eta.*etaB;
        GB=1./(om-lam-V2*SigB + 2*Del)+ 1./(om-lam-V2*SigB - 2*Del); 
        GB = GB./2;
    end

    if geom == 0.5
        SigB=-eta.*etaB;
        GB=1./(om-lam-V2*SigB + 2*Del)+ 2./(om-lam-V2*SigB - Del); 
        GB = GB./3;
    end

    %imaginary parts
    iGB=imag(GB); 
    tt=0;
    tdown=tdown0;
    
    %holon part  
    SigX=zeros(size(om));  
    igX=@(lam,eta,epsf) (om-lam+eta*etaX)*valence + epsf;
    GX=1./(igX(lam,eta,epsf)-V2*SigX);  iGX=imag(GX);
    
    %important stability stuff
    if qq==1
        istable=howoften;
    else
        istable=howoften0;
    end
    ee=0; cont=1; readytogo=0;
    jlast=0; jj=0; clear errarc2 specer1 qarc2 lamarc2 pnt minsp maxsp becarc1

    
    while cont
        jj=jj+1;        
        %dec=dec0+(T/T0)*0.25;
        dec = dec0;
        decEtaF = decEta;
        if jj-jlast>0
            mixing=0.5;
            eta=max(eta/decEtaF,T/etaConst);
            %eta=min(1e-3,T/10000);            
        end
        stable=(jj>=istable);
        if tdown
            tt=tt+1;
            if tt>1
                if tdown==1
                    T=T/dec;
                    err=err/dec;
                end
                if tdown==-1
                    T=T*dec;
                    err=dec*err;
                end
                eta=max(err,eta);
            end
            tdown=0;
            ee=0; howoften=howoften0;
            if dispflag
                disp('------------------');
            end
            disp([num2str(tt) ') T/TK=' num2str(T/TK) ' & T=' num2str(T)]);
            if dispflag
                disp('------------------');
            end
            %------------------------------------------------------------
            fom=f(cut(om/T)); nBom=nB(cut(om/T));
            dfom=dfomfun(cut(om/T),T);
            dnBom=dnBomfun(cut(om/T),T);
            etaB=yhilbert(1-2*fom);
            %etaB=1i*(0.5-fom);

            tgc=gc;
            tgc=tgc./sum(-imag(tgc).*dom/pi);
            %plot(om,-imag(gc),om,-imag(tgc),'.'); pause
            %---------
            if logridflag
                fBmX=f(cut(BmX/T)); fmXmB=transpose(fBmX);
                dfomBmX=dfomfun(cut(BmX/T),T);
            end
        end

        %-----------
        %-----------
        GBold=GB; GXold=GX; becold=bec;
        %-------- ConvX1 ------------
        SigX3=0; SigX4=0;
        if logridflag
            SigX1=-(GB.*dom)*(imag(gcBmX).*fBmX)/pi;
            SigX2=-(iGB.*nBom.*dom)*gcmBmX/pi;
        else
            SigX1=-conv(GB.*dom,fliplr((fom).*imag(tgc)),'same')/pi;
            SigX2=-conv(iGB.*(nBom).*dom,tgc,'same')/pi; %This was gc2=gc
        end
        
        SigX=mixing*SigX+(1-mixing)*(SigX1+SigX2);
        
        if logridflag
            SigB1=gamma*(GX.*dom)*(imag(gcmXmB).*fmXmB)/pi;
            SigB2=-gamma*(iGX.*(1-fom).*dom)*gcmXmB/pi;
        else
            SigB1=gamma*conv(GX.*dom,imag(tgc).*(fom),'same')/pi;
            SigB2=-gamma*conv(iGX.*fliplr(fom).*dom,tgc,'same')/pi;
        end
        
        SigB=mixing*SigB+(1-mixing)*(SigB1+SigB2);
        
        SigBbr=SigB-eta*etaB; % add broadening...
        
        %implement broadening
        SigB = SigBbr; 
        
        if J<1e-4
            SigB=zeros(size(GB));
            GX=zeros(size(GX));
        end
        
        %---------- Green's functions definition ---------------
        %----------------------
        %---------------- The lambda loop ----------------------
        %----- what is mimlam?
        
        if 1
            if geom==0
                minlam=-V2*real(sum(SigB(izom))/lzom);
            end
            if geom==0.5
                minlam=-V2*real(sum(SigB(izom))/lzom) + 2*Del;
            end

            if geom==0.1
                minlam=-V2*real(sum(SigB(izom))/lzom) + 2*Del;
            end
        

            
            qerr=1;
            maxk=100;
            lam1=minlam;  lam2=max(D,T); kk=0; q=1; qX=0;
            while (abs(lam-minlam)>1e-8) || abs(q-q0)>1e-8
                kk=kk+1;
                if kk==2*maxk
                    %disp('error converging...');
                    %disp([log10(eta1) log10(err))])
                    %keyboard
                    break;
                end
                
                lam=(lam1+lam2)/2;
                if geom==0
                    GB=1./(om-lam-V2*SigB);
                elseif geom==0.1
                    GB=1./(om-lam-V2*SigB - 2*Del) + 1./(om-lam-V2*SigB + 2*Del);
                    GB = GB./2;
                elseif geom==0.5
                    GB=2./(om-lam-V2*SigB - Del) + 1./(om-lam-V2*SigB + 2*Del);
                    GB = GB./3;
                end
                iGB=imag(GB);
                sumruleBnow=sum(-iGB.*dom)/pi;
                %GB=GB/sumruleBnow; 
                iGB=imag(GB);                
                qB=-sum(iGB.*nBom.*dom)/pi;
                %--------
                GX=1./(igX(lam,eta,epsf)-V2*SigX); iGX=imag(GX);
                sumruleXnow=sum(-iGX.*dom)/pi;
                %GX=GX/sumruleXnow; 
                iGX=imag(GX);
                qX=-gamma*sum(iGX.*fom.*dom)/pi;
                %---------
                q=qB+valence*qX;
                qerr=q-q0;
                %----------------------
                if qerr<0 && q>0
                    lam2=lam;
                else
                    lam1=lam;
                end
                qarc1(kk)=q;
                lamarc1(kk)=lam;
            end
        end
        
        indnan=[find(isnan(GX)) find(isnan(GB))];
        if find(indnan)
            stable=1;
            break;
        end
        flagB=iGB.*om>0; flagX=iGX>0;
        wrongB=sum(abs(flagB)); wrongX=sum(abs(flagX));
        if wrongB+wrongX && 0
            disp('spectral problem'); keyboard
            GB(flagB)=0;
            GX(flagX)=0;
        end
        %-------------------------------------------------
        specerr1(jj,:)=(abs(GB-GBold))./sum(abs(GBold).*dom)+(abs(GX-GXold))./sum(abs(GXold).*dom)+abs(bec-becold)*0.1;
        error=sum(specerr1(jj,:));
        errarc2(jj)=error;
        disp([num2str(qq) ',' num2str(tt) ',' num2str(jj) '(' num2str(jj-jlast) '/' num2str(istable-jlast) '). log10_error=' num2str(log10(error))  ', log10(eta)= ' num2str(log10(eta)) ', log10(V2*eta)= ' num2str(log10(V2*eta))]);
        if (log10(error))<-3 && jj-jlast>(howoften0/5)  % even if it is not stable, assume so...
            stable=~~1; readytogo=1;
        end
        regrid=0;
        converge=0;
        if stable
            %converge=mean(errarc2(ind2))<1e-3;
            converge=1;
            %readytogo=1;
        end
        if jj==maxIterations1
            eta=eta*50;
        end
        if jj>maxIterations2
            if qq==1
                tt=tt-1;
                T=TT(qq,tt);
                %------------------------------------------------------------
                fom=f(cut(om/T)); nBom=nB(cut(om/T));
                dfom=dfomfun(cut(om/T),T);
                dnBom=dnBomfun(cut(om/T),T);
                %---------
                %                 fBmX=f(cut(BmX/T)); fmXmB=transpose(fBmX);
                %                 dfomBmX=dfomfun(cut(BmX/T),T);
                %plotting the phase shift at the end
                %semilogx(TT,deltaX/pi); pause
            end
            cont=0;
            tdown=0;
            %T=0;
            %qq=length(x0);
            %SigB=SigBarc(qq,:);
            %SigX=SigXarc(qq,:);
            %GB=GBarc(qq,:);
            %GX=GXarc(qq,:);
        end
        if (converge && ~~readytogo && T>0) % ready to record and move on
            tdown=tdown0;

            
            if T/dec<Tmin  || T*dec>Tmax %|| qq~=1
                cont=0; %don't continue, you've reached the limit
                tdown=0;
            end

            phase_shift_test = -imag(log(-lzom./sum(GX(izom))))/pi;
            %disp(num2str(phase_shift_test))
            if (phase_shift_test >= maxphaseshift)
                cont=0;
                tdown=0;                
            end

            readytogo=0;
            %aminsp(qq,tt)=mean(minsp(ind));
            %amaxsp(qq,tt)=mean(maxsp(ind));
            jlast=0; jj=0; clear errarc2 specerr1 qarc2 lamarc2 pnt minsp maxsp becarc1
            howoften=howoften0; istable=howoften;
            dGXr=gradient(real(GX))./dom;  dGBr=gradient(real(GB))./dom;
            
            TT(qq,tt)=T;
            errtarc(tt)=error;
            qarc(qq,tt)=q;
            
            
            qL(qq,tt)=-sum((iGB.*nBom.*(om<0).*dom))/pi;
            qR(qq,tt)=-sum((iGB.*nBom.*(om>0).*dom))/pi;
            mv(qq,tt)=-gamma*sum((iGX.*fom.*dom))/pi/q0;
            
            IRleak(qq,tt)=sum(abs(iGB(IRind).*dom(IRind))/sum(abs(iGB.*dom))+abs(iGX(IRind).*dom(IRind))/sum(abs(iGX.*dom)))/2;
            UVleak(qq,tt)=sum(abs(iGB(UVind).*dom(UVind))/sum(abs(iGB.*dom))+abs(iGX(UVind).*dom(UVind))/sum(abs(iGX.*dom)))/2;
            
            lameff(qq,tt)=lam+V2*real(SigB(izom));
            iJeff(qq,tt)=1/J+V2*real(SigX(izom));
  
            
            %sumruleB(qq,tt)=sumruleBnow;
            sumruleB(qq,tt) = sum(-imag(GB).*dom)/pi;
            %sumruleX(qq,tt)=sumruleXnow;
            sumruleX(qq,tt) = sum(-imag(GX).*dom)/pi;
            deltaX(qq,tt)=-imag(log(-1./GX(izom)));
            
            if dispflag
                disp(['area (sum rule) of X,B=' num2str(sumruleX(qq,tt)) ' , ' num2str(sumruleB(qq,tt))]);
                disp(['q=' num2str(q)]);
                disp(['qB/q=' num2str(qB/q)]);
                disp(['qX/q=' num2str(qX/q)]);
                disp(['deltaX/pi=' num2str(deltaX(qq,tt)/pi)])
                disp('observables...')
            end
            
            
            errarc(tt)=error;
            lamarc(tt)=lam;
            
            if find(isnan(iGB+iGX))
                disp('NaN!!')
                keyboard
            end
                        
            expX(qq,tt)=sum(dom.*fom.*imag(GX.^2))/pi;
            expB(qq,tt)=sum(dom.*nBom.*imag(GB.^2))/pi;
            
            %ichic{tt}(qq,:)=conv(fom.*iGX,fliplr(iGX),'same')*dom/pi-conv(fliplr(fom.*iGX),iGX,'same')*dom/pi;
            
            deltaomX=gradient(SigX)./gradient(om).*GX;
            chicX(qq,tt)=sum(fom.*imag(deltaomX.^2).*dom)/pi;
            deltaomB=gradient(SigB)./gradient(om).*GB;
            chicB(qq,tt)=sum(nBom.*imag(deltaomB.^2).*dom)/pi;
            
            %ichic{qq}(tt,:)=(ichic{qq}(tt,:)-fliplr(ichic{qq}(tt,:)))/2;
            
            
            
            if geom==0
                suscBloc=GB.^2;
                suscBzq = suscBloc;
                suscBfq = suscBloc; 
                %contribution to the entropy
                calSB=log(-1./GB);
                suscpB=suscBloc;
                suscpBpar=0; suscpBperp=0;
                JH(qq,tt) = 0.0;
            end
            
            if geom==0.5
                %need GB with q: to do that need Del, SigB and lam
                GBks=@(y) 1./(om + 2*Del*cos(2*pi*y/3) - lam - V2*SigB);
                epks=@(y) 2*cos(2*pi*y/3);
                suscBzq=GBks(0).^2 + 2*GBks(1).*GBks(1);
                suscBfq=2*GBks(0).*GBks(1) + GBks(1).*GBks(1);
                suscBloc=GB.^2;
                %contribution to the entropy
                calSB = 0.0;
                preJh = 0.0;
                for sumn=-1:1
                    %suscBzq = suscBzq + (GBks(sumn).^2)./3;
                    calSB = calSB + log(-1./GBks(sumn))./3;
                    preJh = preJh + epks(sumn)*imag(GBks(sumn))/3;
                end
                suscpB=suscBloc;
                suscpBpar=0; suscpBperp=0;
                %JH calculation
                sumOcGB = sum(dom.*nBom.*preJh./(2*pi));
                JH_0=-Del./sumOcGB;
                %JH with an xi to lift the 1st order
                %xi_1st = 2;
                %polynom = [xi_1st*Del^2 0 1 -sumOcGB];
                %polynom2 = @(x) (x.^3)*xi_1st*Del^2 + x - sumOcGB;
                %oneOvJH = roots(polynom);
                %compare_old = abs(oneOvJH - sumOcGB);
                %ind_min = find(compare_old == min(compare_old));
                %oneOvJH = fzero(polynom2, sumOcGB);
                %JH(qq,tt) = 1./real(oneOvJH(ind_min));
                JH(qq,tt) = JH_0;
            end      
            
            if geom==0.1
                %need GB with q: to do that need Del, SigB and lam
                GBks=@(y) 1./(om + 2*Del*cos(pi*y) - lam - V2*SigB);
                epks=@(y) 2*cos(pi*y);
                suscBzq=0;
                suscBloc=GB.^2;
                %contribution to the entropy
                calSB = 0.0;
                preJh = 0.0;
                for sumn=0:1
                    suscBzq = suscBzq + (GBks(sumn).^2)./2;
                    calSB = calSB + log(-1./GBks(sumn))./2;
                    preJh = preJh + epks(sumn)*imag(GBks(sumn))/2;
                end
                suscpB=suscBloc;
                suscpBpar=0; suscpBperp=0;
                %JH calculation
                sumOcGB = sum(dom.*nBom.*preJh./(2*pi));
                JH_0=-Del./sumOcGB;
                %JH with an xi to lift the 1st order
                %xi_1st = 2;
                %polynom = [xi_1st*Del^2 0 1 -sumOcGB];
                %polynom2 = @(x) (x.^3)*xi_1st*Del^2 + x - sumOcGB;
                %oneOvJH = roots(polynom);
                %compare_old = abs(oneOvJH - sumOcGB);
                %ind_min = find(compare_old == min(compare_old));
                %oneOvJH = fzero(polynom2, sumOcGB);
                %JH(qq,tt) = 1./real(oneOvJH(ind_min));
                JH(qq,tt) = JH_0;
            end      
            
            %susceptibilities calculations
            %suscpar(qq,tt)=sum((nBom+1/2).*imag(suscpBpar).*dom)/pi;
            %suscperp(qq,tt)=sum((nBom+1/2).*imag(suscpBperp).*dom)/pi;
            
            %
            suscploc(qq,tt)=sum((nBom).*imag(suscBloc).*dom)/pi;
            % zero momentum susc relates to total spin susceptibility
            suscpzq(qq,tt)=sum((nBom).*imag(suscBzq).*dom)/pi;
            suscpfq(qq,tt)=sum((nBom).*imag(suscBfq).*dom)/pi;
            %the moment
            tnBom = T*nBom;
            momentMuLoc(qq,tt)= sum((tnBom).*imag(suscBloc).*dom)/pi;
            momentMuBzq(qq,tt)= sum((tnBom).*imag(suscBzq).*dom)/pi;
            
            %%---------------
            %enropy calculations
            %%---------------
            
            
            if logridflag
                SC =-gamma*real(GB.*dom)*(dfomBmX.*imag(gcBmX))*transpose(iGX.*fom.*dom)/pi/pi+gamma*(nBom.*iGB.*dom)*(dfomBmX.*imag(gcBmX))*transpose(real(GX).*dom)/pi/pi;
                d2fomBmX=d2fomfun(cut(BmX/T),T);
                Res(qq,tt)=2*pi*(iGB.*dom)*(d2fomBmX)*transpose(iGX.*fom.*dom)/pi/pi+2*pi*(nBom.*iGB.*dom)*(d2fomBmX)*transpose(iGX.*dom)/pi/pi;

                Sig2B=gamma*((GX.*dom)*transpose(d2fomBmX.*imag(gcBmX)+fBmX.*imag(gc2BmX))-(iGX.*fliplr(fom).*dom)*transpose(d2fomBmX))/pi;
                Sig2X=((GB.*dom)*(d2fomBmX.*imag(gcBmX)+fBmX.*imag(gc2BmX))-(iGB.*nBom.*dom)*conj(gc2BmX))/pi;
                calGB=GB.*Sig2B; calGX=GX.*Sig2X;
            else
                calGX=(gradient(SigX)./gradient(om)).*GX;
                F1=imag(calGX).*fom; F2=calGX;
                chic=conv(fliplr(F1),F2,'same')*dom/pi; chic=chic+conj(fliplr(chic));
                pchicX(qq,:)=chic;
                calGB=(gradient(SigB)./gradient(om)).*GB;
                F1=imag(calGB).*nBom; F2=calGB;
                chic=conv(fliplr(F1),F2,'same')*dom/pi; chic=chic+conj(fliplr(chic));
                pchicB(qq,:)=chic;
                
                SigC1=conv((nBom).*iGB,fliplr(conj(GX)),'same')*doml/pi;
                SigC2=-conv(fliplr(fom.*iGX),GB,'same')*doml/pi;
                SigC=SigC1+SigC2;
                SC = gamma*doml/pi*sum(dfom.*imag(gc).*real(SigC));

            end
            NdeltaC(qq,tt)=-pi*rho*(sum(dom.*nBom.*iGB.*real(GX)-dom.*fom.*iGX.*real(GB))/pi); % This is SigC(0+i\eta)
            %attempt to make a T matrix
            %SigX1=-(GB.*dom)*(imag(gcBmX).*fBmX)/pi;
            %SigX2=-(iGB.*nBom.*dom)*gcmBmX/pi;
            %tmat1 = -(GB.*dom)*(imag(GX(BmX)).*fBmX)/pi;
            %tmat = tmat1 + tmat2;
            NdZ(qq,tt)=-sum(dom.*nBom.*imag(GB).*dGXr+dom.*fom.*imag(GX).*dGBr)/pi; % This is \partial_\omega \Sig_C(0+i\eta)
            
            
            %---------------- Calculating SB ----------------------
            
            SB=-sum(dnBom.*imag(calSB).*dom)/pi;
            %disp(SB)
            %-----------------------------------------------------
            SX=-gamma*sum(dfom.*imag(log(-1./GX-1i*8)).*dom)/pi;
            
            SB1=-sum(dnBom.*real(GB).*imag(SigB).*dom)/pi;
            SX1=-gamma*sum(dfom.*real(GX).*imag(SigX).*dom)/pi;
            
            SB2=-sum(dnBom.*iGB.*real(SigB).*dom)/pi;
            SX2=-gamma*sum(dfom.*iGX.*real(SigX).*dom)/pi;
            
            S = SB+SX+V2*(SB1+SX1)+SC;
            % or 
            %S = SB+SX+SB1+SX1+SC;
            
            %Res(qq,tt)=2*pi*(iGB*dom)*(d2fomBmX)*transpose(iGX.*fomX.*domX)/pi/pi+2*pi*(nBomB.*iGB*dom)*(d2fomBmX)*transpose(iGX.*domX)/pi/pi;
            Sarc(qq,tt)=S;
            SBarc(qq,tt) = SB;
            SXarc(qq,tt) = SX;
            SCarc(qq,tt) = SC;
            SB1arc(qq,tt) = SB1;
            SX1arc(qq,tt) = SX1;
            
                                 
            
            GBtarc(tt,:)=GB;
            GXtarc(tt,:)=GX;
            SigBtarc(tt,:)=SigB;
            SigXtarc(tt,:)=SigX;
            %tscatteringarc(tt,:) = tmat;
            
            GBarc(qq,:)=GB;
            GXarc(qq,:)=GX;
            SigBarc(qq,:)=SigB;
            SigXarc(qq,:)=SigX;
            
            %-----------------------------
            iGB1=iGB.*dom; iGX1=iGX.*dom;
            IRleak(qq,tt)=sum(abs(iGB1(IRind))/sum(abs(iGB1))+abs(iGX1(IRind))/sum(abs(iGX1)))/2;
            UVleak(qq,tt)=sum(abs(iGB1(UVind))/sum(abs(iGB1))+abs(iGX1(UVind))/sum(abs(iGX1)))/2;
            if dispflag
                disp(['IR,UV leaks:' num2str(IRleak(qq)) ', ' num2str(UVleak(qq))]);
            end
            %-----------------------------
            
            
            qBarc(qq,tt)=qB;
            qXarc(qq,tt)=qX;
            xarc(qq,tt)=x;
            yarc(qq,tt)=T/TK;
            qLarc(qq,tt)=-sum((iGB.*nBom.*(om<0).*dom))/pi;
            qRarc(qq,tt)=-sum((iGB.*nBom.*(om>0).*dom))/pi;
            lamarc(qq,tt)=lam;
            
            TKarc(qq)=TK;
            Jarc(qq)=J;
            Barc(qq)=B;
            q0arc(qq)=q0;
            T0arc(qq) = T0;
            
            

        end
    end
    
    saveDelta = deltaX(qq,:);
    saveEpsf = epsf;
    saveq = q;
    saveqX = qXarc(qq,:);
    saveqB = qBarc(qq,:);
    saveT = TT(qq,:);
    saveSusc = suscpzq(qq,:);
    saveSuscAlt = suscpfq(qq,:);
    saveSuscLoc = suscploc(qq,:);
    saveS = Sarc(qq,:);
    saveSB = SBarc(qq,:);
    saveSX = SXarc(qq,:);
    saveJ = JH(qq,:);
    saveLameff = lameff(qq,:);
    saveLam = lamarc(qq,:);
    saveOm = om;
    saveQplus = qRarc(qq,:);
    saveQminus = qLarc(qq,:);
    momentLoc = momentMuLoc(qq,:);
    momentzq = momentMuBzq(qq,:);

    name3 =[name '_data_x0_' num2str(fix(x)) '-' num2str(fix(rem(x,1)*1000)) '_' num2str(fix(gammaFactor)) '-' num2str(fix(rem(gammaFactor,1)*1000)) '_' num2str(fix(epsf)) '-' num2str(fix(rem(epsf,1)*1000)) ];
    save(name3,'saveOm','saveT','saveq','saveEpsf','saveDelta','saveqX','saveqB','saveSusc','saveSuscAlt','saveSuscLoc', 'saveS', 'saveSB','saveSX','saveJ','saveLameff','saveLam','GBtarc', 'GXtarc','SigBtarc', 'SigXtarc', 'saveQplus','saveQminus','momentLoc','momentzq')
    
    toc
    disp(['================== ' num2str(qq) ' >'])
   
    
end
%-----------------------
save(name2,'-regexp','^(?!(BmX|XmB|gcBmX|gcmBmX|gcmXmB|gc2BmX|fBmX|fmXmB|dfomBmX|f2fomBmX|gc2BmX)$).')      