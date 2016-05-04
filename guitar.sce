
//init(44100,0.001,110,4,0.1,1,[0.3 0.7],[100,10;1000,8])

function init(SR,B,f,TF,x0,c0,rp,loss)
    mu=1;
    N=50;
    
    k=1/SR;
    Gamma2=(2*f)^2;
    kappa2=((2*f*sqrt(B))/%pi)^2;
    T=mu*(f^2);
    omega_1=loss(1,1)*2*%pi;
    T60_1=loss(1,2);
    omega_2=loss(2,1)*2*%pi;
    T60_2=loss(2,2);
    xi_1=(-Gamma2+sqrt((Gamma2^2)+4*kappa2*omega_1*omega_1))/(2*kappa2);
    xi_2=(-Gamma2+sqrt((Gamma2^2)+4*kappa2*omega_2*omega_2))/(2*kappa2);
    

    sigma_0=((6*log(10))/(xi_2-xi_1))*((xi_2/T60_1)-(xi_1/T60_2));
    sigma_1=((6*log(10))/(xi_2-xi_1))*(-(1/T60_1)+(1/T60_2));

    TMP1=eye(N-1,N-1)-3*eye(N-1,N-1);
    TMP2=TMP1(1,:);
    TMP2(2)=1;
    TMP1=TMP1(:,1);
    TMP1(2)=1;
    Dxx=(N^2)*toeplitz(TMP1,TMP2);
    mata=(1+(sigma_0*k))*eye(N-1,N-1)-sigma_1*k*Dxx;
    matb=-2*eye(N-1,N-1)-(Gamma2*k*k*Dxx)+kappa2*k*k*Dxx*Dxx;
    matc=(1-(sigma_0*k))*eye(N-1,N-1)+sigma_1*k*Dxx;
    //SolveG is very slow because of the plots, use Solve for sound
    //SolveG(mata,matb,matc,x0,c0,N,rp,SR);
    Solve(mata,matb,matc,x0,c0,N,rp,SR);
endfunction

function xi(omega,kappa2,Gamma2)
    return (-Gamma2+sqrt((Gamma2^2)+4*kappa2*omega*omega))/(2*kappa2);
endfunction

function SolveG(mata,matb,matc,x0,c0,N,rp,SR)
    NF=100000;
    u_nm1=zeros(N-1,1);
    h=1/N;
    x=h;i=1;
    while x<x0
        u_nm1(i)=(c0/x0)*x;
        x=x+h;
        i=i+1;
    end
    
    x=x0;
    while x<1
        u_nm1(i)=c0*(x-1)/(x0-1);
        x=x+h;
        i=i+1;
    end
    u_n=u_nm1;
    
    x=0;i=0;
    while(x<rp(1))
        x=x+h;
        i=i+1;
    end
    a1=i-1;
    a2=i;   
    while(x<rp(2))
        x=x+h;
        i=i+1;
    end
    a3=i-1;
    a4=i;
    
    u_np1=zeros(N-1,1);
    //T=chol(mata);
    out = zeros(2,NF)
    smata=sparse(mata)
    for i=1:NF
        drawlater;
//        y=inv(T)*((-matb*u_n)-(matc*u_nm1));
//        u_np1=inv(T')*y;
        //y=Gauss(T,((-matb*u_n)-(matc*u_nm1)));
        //u_np1=Gauss(T',y);
        u_np1=lusolve(smata,(-matb*u_n)-(matc*u_nm1));
        u_nm1=u_n;
        u_n=u_np1;
        clf();
        plot([0,u_nm1',0]');
        a=gca();
        a.data_bounds=[0,-1;N-1,1];
        drawnow;
        p1=u_n(a1);
        p2=u_n(a2);
        out(1,i)=p1+(p2-p1)*(rp(1)-h*a1)/h;
        p1=u_n(a3);
        p2=u_n(a4);
        out(2,i)=p1+(p2-p1)*(rp(2)-h*a3)/h;
    end
    plot(1:NF,out);
    playsnd(out,SR);
    disp(fft(out(:,1)));
endfunction

function Solve(mata,matb,matc,x0,c0,N,rp,SR)
    NF=100000;
    u_nm1=zeros(N-1,1);
    h=1/N;
    x=h;i=1;
    while x<x0
        u_nm1(i)=(c0/x0)*x;
        x=x+h;
        i=i+1;
    end
    
    x=x0;
    while x<1
        u_nm1(i)=c0*(x-1)/(x0-1);
        x=x+h;
        i=i+1;
    end
    u_n=u_nm1;
    
    x=0;i=0;
    while(x<rp(1))
        x=x+h;
        i=i+1;
    end
    a1=i-1;
    a2=i;   
    while(x<rp(2))
        x=x+h;
        i=i+1;
    end
    a3=i-1;
    a4=i;
    
    u_np1=zeros(N-1,1);
    out = zeros(2,NF)
    smata=sparse(mata)
    for i=1:NF
        u_np1=lusolve(smata,(-matb*u_n)-(matc*u_nm1));
        u_nm1=u_n;
        u_n=u_np1;

        p1=u_n(a1);
        p2=u_n(a2);
        out(1,i)=p1+(p2-p1)*(rp(1)-h*a1)/h;
        p1=u_n(a3);
        p2=u_n(a4);
        out(2,i)=p1+(p2-p1)*(rp(2)-h*a3)/h;
    end
    plot(1:NF,out);
    playsnd(out,SR);
    plot(fftshift(fft(out(1,:))));
endfunction

function y=Gauss(A,b)
    sA=size(A);
    sA=sA(1);
    y=zeros(sA,1);
    s=0;
    for i=sA:-1:1
        for j=sA:-1:i
            s=s+A(i,j)*y(j);
        end;
        y(i)=(b(i)-s)/A(i,i);
        s=0;
    end;
endfunction
