c 6-7-2020
       parameter(Nb=1000,Mdim=1000,Nff=100,Nc=Nb/Nff,Nt=2000)
       implicit real*8(a-h,o-z)
      complex*16 ci,cs,eta,czero,cx,CPOT,comPOT,coeff
       real*8 X(Mdim,Mdim),T(Mdim,Mdim),X_grid(Mdim)
      complex*16 H(nb,nb),E(nb),C(Nb,Nb),D(Nb,Nb)
      real*8 Hr(Nb,Nb),Hi(Nb,Nb),Er(Nb),Ei(Nb)
      real*8 VecR(Nb,Nb),VecI(Nb,Nb)
      real*8 V1(Nb),v2(Nb),V3(Nb)
      integer*4 space,Indx_i(Nb),Indx_ic(Nb)
      complex*16 W(Nb,Nb),W2(Nb,Nb),EX(Nb,Nb)
      complex*16 U0(Nff,Nff),U1(Nff,Nff)
       real*8 Ur(Nff,Nff),Ui(Nff,Nff),QEr(Nff),QEi(Nff)
      real*8 QEVR(Nff,Nff),QEVI(Nff,Nff),QEv1(Nff),QEv2(Nff),QEv3(Nff)
      complex*16 CQE(Nff,Nff),QE(Nff)

       write(60,*) 'Nb=1000,Mdim=Nb,Nff=100,Nc=Nb/Nff,Nt=2000)',
     &Nb,Mdim,Nff,Nc,Nt
      ci=(0.d0,1.d0)
      czero=(0.d0,0.d0)
      pi=4.d0*datan(1.d0)
      theta=0.1d0
      eta=cdexp(ci*theta)
      omega=0.057d0*10
      Time=2.d0*pi/omega
      eps0=0.0015d0
c     eps0=0.d0
      Dt=Time/dfloat(Nt)
        ispecFLOWUET=1
c =1 calculate spectrum of Floquet NOT by ttp and stop
c =0 do bot calculate spectrum of Floquet NOT by ttp
c
c DVRham
c     DVR using sin as a basis set
       hbar=1.d0
       space=1
       am=1.d0
       facKIN=hbar**2/(2.d0*am)
      write(6,*) 'facKIN=',fackin
       Abox=50.d0
      call GMAT(Mdim,Mdim,Pi,Abox,X)
      write(6,*) ' GMAT was DONE'
            ivec=1
      call RS(Mdim,Mdim,X,X_grid,ivec,T,v1,v2,ierr)
        if(ierr.ne.0) write(6,*) ' end of RS ierr = ',ierr
       write(6,*) ' RS was done'
             if(space.eq.1) then
         do i=1,Mdim
      X_grid(i)=X_grid(i)-Abox/2.d0
         end do
             end if
           do k=1,Mdim
       cx=x_grid(k)*(1.d0,0.d0)
       cs=CPOT(cx)
       write(100,*) x_grid(k),dreal(cs)
         end do
c
        do i=1,Nb
        do j=1,i
        cs=czero
           do k=1,Mdim
       cx=x_grid(k)*eta
       TT=T(i,k)*T(j,k)
       comPOT=CPOT(cx)
       cs=cs+TT*comPOT
           end do
      H(i,j)=cs
c     if(i.eq.j)H(i,j)=H(i,j)+facKIN*(dfloat(i)*pi/Abox/ETA)**2
      if(i.eq.j)H(i,j)=H(i,j)+0.5d0*(dfloat(i)*pi/Abox/ETA)**2
       H(j,i)=H(i,j)
       end do
       end do
c
       ivec=1
      call  CEIGEN
     &(H,C,E,Hr,Hi,Er,Ei,ivec,VecR,VecI,v1,v2,v3,Nb,Nb)
       if(ivec.eq.0) then
       call ORDERC(E,Nb,Nb)
         do i=1,Nb
           write(10,*) dreal(E(i)),dimag(E(i))
         cs=cdexp(-ci*E(i)*Time/hbar)
           write(11,*) dreal(cs),dimag(cs)
          end do
       else
         call CORDER2(E,C,nb,nb)
         do i=1,Nb
           write(10,*) dreal(E(i)),dimag(E(i))
         cs=cdexp(-ci*E(i)*Time/hbar)
           write(11,*) dreal(cs),dimag(cs)
          end do
          call ProdXrealMat(X,C,W,Nb)
          call ProdTranMat(C,W,D,Nb)
c D dipole transiton matrix
          end if
c END SPECTRUM OF FIELD FREE HAMILTONIAN
         if(ivec.eq.0) stop
         do i=1,Nff
         write(200,*) dreal(E(i)),dimag(E(i))
         end do
c
c   below calculation of Floquet matrix
c
           do i=1,Nb
           do j=1,Nb
          H(i,j)=czero
           end do
           end do
          ii=0
          do ic=-Nc/2+1,+Nc/2
          do i=1,Nff
          ii=ii+1
          H(ii,ii)=E(i)+ic*hbar*omega
           Indx_i(ii)=i
           Indx_ic(ii)=ic
           write(300,*) dreal(H(ii,ii)),dimag(H(ii,ii))
           jj=0
           do icp=-Nc/2+1,+Nc/2
           do ip=1,Nff
           jj=jj+1
           iii=icp-ic
           if(iii.lt.0) iii=-iii
           if(iii.eq.1) then
           H(ii,jj)=eps0*eta*D(i,ip)/2.d0
           end if
          end do
          end do
          end do
          end do
c
          if(ii.ne.jj.and.ii.ne.Nb) then
          write(12,*) 'GEVALD AND STOP'
          stop
          end if
          if(ivec.eq.0) stop
c
c BELOW CALCULATE SPECTRUM of FLOQUET not by ttp
        if(ispecFLOWUET.eq.1) then
       ivec=0
      call  CEIGEN
     &(H,C,E,Hr,Hi,Er,Ei,ivec,VecR,VecI,v1,v2,v3,Nb,Nb)
       if(ivec.eq.0) then
       call ORDERC(E,Nb,Nb)
         do i=1,Nb
           write(12,*) dreal(E(i)),dimag(E(i))
           cs=cdexp(-ci*E(i)*Time/hbar)
           write(13,*) dreal(cs),dimag(cs)
          end do
       else
         call CORDER2(E,C,nb,nb)
         do i=1,Nb
           write(12,*) dreal(E(i)),dimag(E(i))
           cs=cdexp(-ci*E(i)*Time/hbar)
           write(13,*) dreal(cs),dimag(cs)
          end do
          end if
          stop
           end if
c
c    BELOW CALCULATE Exp(-iH_floquet*dt/hbar)
           cfac=-ci*dt/hbar
           do i=1,Nb
           do j=1,Nb
          EX(i,j)=czero
          W(i,j)=H(i,j)
          if(i.eq.j) EX(i,i)=1.d0
           end do
           end do
          kk=1
                    do k=1,20
          kk=kk*k
           coeff=cfac**k/dfloat(kk)
          if(cdabs(coeff).lt.1.d-10) then
          write(500,*) k,coeff
          if(k.eq.20) stop
          go to 9999
          end if
          call AddMat(EX,W,coeff,Nb)
          call ProdMat(W,H,W2,Nb)
          do i=1,Nb
          do j=1,Nb
          W(i,j)=W2(i,j)
          end do
          end do
 9999     continue
          end do
c  above end EXPONENT of FLOQUET - Floquet evolution operator
          do i=1,Nff
          U0(i,j)=0.d0
          if(i.eq.j) U0(i,i)=1.d0
          end do
          tt=0.d0
          do it=1,Nt
          tt=tt+dt
                         do j=1,Nff
          do  ii=1,Nb
          i=Indx_i(ii)
          ic=Indx_ic(ii)
          cs=czero
          do jj=1,Nb
          ip=Indx_i(jj)
          icp=Indx_ic(jj)
          if(icp.eq.0) then
          cs=cs+cdexp(ci*omega*ic*tt)*EX(ii,jj)*U0(ip,j)
          end if
          end do
          U1(i,j)=cs
          end do
                  do i=1,Nff
           U0(i,j)=U1(i,j)
                  end do
                         end do
c
           end do
          if(tt.ne.TIME)then
          write(6,*) 'GEVALD tt NE Time'
          stop
          end if
c  BELOW CALCULATE FLOQUET SPECTRUM BY ttp
       ivec=0
      call  CEIGEN
     &(U1,CQE,QE,Ur,Ui,QEr,QEi,ivec,QEVR,QEVI,QEv1,QEv2,QEv3,Nff,Nff)
       if(ivec.eq.0) then
       call ORDERC(QE,Nff,Nff)
         do i=1,Nff
           write(130,*) dreal(QE(i)),dimag(QE(i))
          end do
       else
         call CORDER2(QE,CQE,nb,nb)
          CMAX=-1.d+10
         do i=1,Nff
           write(130,*) dreal(QE(i)),dimag(QE(i))
          if(cdabs(CQE(1,i)).gt.CMAX) then
          CMAX=cdabs(CQE(1,i))
          imax=i
          end if
           end do
           write(131,*) imax,QE(imax),CMAX
          end if
c
           stop
         end


      SUBROUTINE GMAT(MDIM,M,PI,A,GR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION GR(MDIM,MDIM)
      write(6,*) ' Abox pi ',a,pi
      do i=1,m
         do j=1,m
            gr(I,J)=0.D0
         END DO
      END DO
      C=2.D0*A/PI**2
      DO 1 I=1,M
         GR(I,I) =A/2.D0
         IF(I.NE.1) THEN
            DO 2 J=1,I-1
               GR(I,J)= 0.D0
               IJ=I+J
               IJM=I-J
              K=MOD(IJ,2)
               IF(K.NE.0) THEN
                  BB = 1.D0/IJ**2 - 1.D0/IJM**2
                  GR(I,J)=C*BB
               END IF
               GR(J,I)=GR(I,J)
 2          CONTINUE
         END IF
 1    CONTINUE
      RETURN
      END


        function CPOT(cx)
        implicit real*8 (a-h,o-z)
        complex*16 CPOT,cx
        cpot=-0.63d0*cdexp(-0.142d0*cx**2)
c       cpot=(cx**2/2.d0-0.8d0)*cdexp(-0.1d0*cx**2)
        return
        end

          subroutine ProdXrealMat(A,B,C,N)
          implicit real*8(a-h,o-z)
          real*8 A(N,N)
          complex*16 B(N,N),C(N,N),cs
          do i=1,N
          do j=1,N
          cs=(0.d0,0.d0)
          do k=1,N
          cs=cs+A(i,k)*B(k,j)
          end do
          C(i,j)=cs
          end do
          end do
                    return
          end

          subroutine ProdMat(A,B,C,N)
          implicit real*8(a-h,o-z)
          complex*16 A(N,N)
          complex*16 B(N,N),C(N,N),cs
          do i=1,N
          do j=1,N
          cs=(0.d0,0.d0)
          do k=1,N
          cs=cs+A(i,k)*B(k,j)
          end do
          C(i,j)=cs
          end do
          end do
          return
          end

          subroutine ProdTranMat(A,B,C,N)
          implicit real*8(a-h,o-z)
          complex*16 A(N,N),B(N,N),C(N,N),cs
          do i=1,N
          do j=1,N
          cs=(0.d0,0.d0)
          do k=1,N
          cs=cs+A(k,i)*B(k,j)
          end do
          C(i,j)=cs
          end do
          end do
          return
          end

            subroutine AddMat(A,B,c,N)
            implicit real*8(a-h,o-z)
            complex*16 A(N,N),B(N,N),c
            do i=1,N
            do j=1,N
            A(i,j)=A(i,j)+c*B(i,j)
            end do
            end do
                        return
            end

      SUBROUTINE ORDERC(E,n,Ndim)
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 E(Ndim),ww,e1,e2
C     ********SMALL TO BIG **********************
   1  IX=0
      I=1
   2  E1= E(i)
      E2= E(i+1)
c     IF(dreal(E2).GE.dreal(E1)) GO TO 5
       de1=dabs(dimag(E1))
       de2=dabs(dimag(E2))
           de1=cdabs(E1)
           de2=cdabs(E2)
       de1=dreal(E1)
       de2=dreal(E2)
      IF(dE2.GE.dE1) GO TO 5
      IX=1
      WW=E(I)
      E(I)=E(I+1)
      E(I+1)=WW
  5   I=I+1
      IF((I+1).LE.N)GO TO 2
      IF(IX.GT.0)GO TO 1
      RETURN
      END


      SUBROUTINE CORDER2(E,V,n,ndim)
      IMPLICIT REAL*8(C-H,O-Z)
      complex*16 E(Ndim),V(Ndim,Ndim),ww
      real*8 E1,E2
C     ********SMALL TO BIG **********************
   1  IX=0
      I=1
   2  E1= dreal(E(i))
      E2= dreal(E(i+1))
      IF(E2.GE.E1) GO TO 5
      IX=1
      WW=E(I)
      E(I)=E(I+1)
      E(I+1)=WW
      do j=1,n
      ww=V(j,i)
      V(j,i)=V(j,i+1)
      V(j,i+1)=ww
      end do
  5   I=I+1
      IF((I+1).LE.N)GO TO 2
      IF(IX.GT.0)GO TO 1
      RETURN
      END

      SUBROUTINE ORDER1(E,n,ndim)
      IMPLICIT REAL*8(C-H,O-Z)
      real*8 E(Ndim),ww
      real*8 E1,E2
C     ******** BIG to small **********************
   1  IX=0
      I=1
   2  E1= E(i)
      E2= E(i+1)
      IF(E1.GE.E2) GO TO 5
      IX=1
      WW=E(I)
      E(I)=E(I+1)
      E(I+1)=WW
  5   I=I+1
      IF((I+1).LE.N)GO TO 2
      IF(IX.GT.0)GO TO 1
      RETURN
      END

      SUBROUTINE ORDER2(E,V,n,ndim)
      IMPLICIT REAL*8(C-H,O-Z)
      real*8 E(Ndim),ww
       complex*16 V(Ndim),wc
      real*8 E1,E2
C     ******** BIG to SMALL **********************
   1  IX=0
      I=1
   2  E1= E(i)
      E2= E(i+1)
      IF(E1.GE.E2) GO TO 5
      IX=1
      WW=E(I)
      E(I)=E(I+1)
      E(I+1)=WW
      wc=V(i)
      V(i)=V(i+1)
      V(i+1)=wc
  5   I=I+1
      IF((I+1).LE.N)GO TO 2
      IF(IX.GT.0)GO TO 1
      RETURN
      END




      subroutine CEIGEN
     &(H,T,E,Hr,Hi,Er,Ei,ivec,VecR,VecI,w1,w2,w3,N,Ndim)
      implicit real*8(a-h,o-z)
      complex*16 H(Ndim,Ndim),T(Ndim,Ndim)
      complex*16 E(Ndim),S,V
      real*8 Hr(Ndim,Ndim),Hi(Ndim,Ndim),VecR(Ndim,Ndim)
      real*8 VecI(Ndim,Ndim),w1(Ndim),W2(Ndim),w3(Ndim)
      real*8 Er(Ndim),Ei(Ndim)
      integer C_prod
      do i=1,n
      do j=1,n
      a=dreal(H(i,j))
      b=dimag(H(i,j))
      Hr(i,j)=a
      Hi(i,j)=b
      end do
      end do
      call CG
     & (Ndim,N,Hr,Hi,Er,Ei,Ivec,VecR,VecI,w1,w2,w3,ierr)
        write(6,*) ' end of CG ierr = ',ierr

              if(ivec.eq.0) then
       do i=1,n
        a=er(i)
                b=ei(i)
        E(i)=dcmplx(a,b)
        end do
         else
c       write(6,*) ' number of calculated eigevecotrs',N
c       write(6,*) ((i,er(i),ei(i)),i=1,n)
        do i=1,n
        a=er(i)
        b=ei(i)
        E(i)=dcmplx(a,b)
        end do
        DO 1 K=1,N
        S=(0.D0,0.D0)
        DO 2 I=1,N
        a=VecR(i,K)
        b=VecI(i,K)
        V=dcmplx(a,b)
        T(i,K)=V
        S=S+V*V
   2  CONTINUE
        S=CDSQRT(S)
        DO 3 I=1,N
        V=T(i,k)/S
        T(i,k)=V
   3  CONTINUE
   1  CONTINUE
        end if
        RETURN
        END


      function dconjug(Z)
      implicit real*8(a-h,o-z)
      complex*16 z,dconjug
      a=dreal(z)
      b=dimag(z)
      dconjug=dcmplx(a,-b)
      return
      end
