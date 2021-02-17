     	implicit none
      	integer, parameter :: k12=selected_int_kind(12)
        integer nb1p,chnln1p,nb,nb2,chnln,nc,nop,nop2

      parameter(nb1p=160)       ! number of beads 124
      parameter(chnln1p=40)     ! chain length  31
      parameter(nb=160)          ! number of beads
      parameter(nb2=154)         ! number of beads without glycines
      parameter(chnln=40)       ! chain length
      parameter(nc=2)           ! chain number	
      parameter(nop=nb*nc)      !total number of beads
      parameter(nop2=nb2*nc)    !total number of beads minus glycines
    	
      integer(kind=k12) coll
      real*8 drandm,xr
      integer iflag
      external drandm,srand

      integer aa(nb),identity(nop),bptnr(nop2)
      integer i,j,k,l,iii,iiii,numbeads,jjjj

      real*8 old_rx(nb1p), old_ry(nb1p), old_rz(nb1p)     
      real*8 dave_rx(124), dave_ry(124), dave_rz(124)     	  
      real*8 new_rx(nb2), new_ry(nb2), new_rz(nb2)  
      real*8 xtmp(nop2),ytmp(nop2),ztmp(nop2),xtmp2(nop2),ytmp2(nop2),ztmp2(nop2)              
      real*8 sv(6,nop),boxl_orig,t
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      real*8 dislocx(nb2),dislocy(nb2),dislocz(nb2)
      real sig_ij,sig_max,bds(400,400),wel(400,400),bdtemp,wltemp
 
      real*8 rmin,r3,rca,rnh,rco,move(chnln),DS1, DS2, DS3                              
      real*8 rcax,rcay,rcaz,rnhx,rnhy,rnhz,rcox,rcoy,rcoz
      real*8 del,pi,setemp,reducedtemp
      real*8 drca(20),drnh(20),drco(20)
      real*8 del_rca(20),del_rnh(20),del_rco(20)
      real*8 bdln(chnln),bl_rn(chnln),bl_rc(chnln)
      real*8 del_bdln(chnln),del_blrn(chnln),del_blrc(chnln)
      real*8 scalerca,rx2,ry2,rz2
      real*8 d1,d2,d3,d1min,d2min,d3min
      real*8 rxij,ryij,rzij,rijsq,boxl
      REAL*8 SQZ610(5,20)
      REAL*8 SZ6,SZ7,SZ8,SZ9,SZ10
      REAL*8 SZ6X,SZ7X,SZ8X,SZ9X,SZ10X
      REAL*8 SZ6Y,SZ7Y,SZ8Y,SZ9Y,SZ10Y
      REAL*8 SZ6Z,SZ7Z,SZ8Z,SZ9Z,SZ10Z
      INTEGER ISZ,ID
	  
      real*8 bm(nop),bmass_temp,bmass(28)                      
      real*8 sumvel,tred,const,sumx,sumy,sumz,x1,x2,y1,y2       

      iflag = 1058472402
      xr=drandm(iflag)
      call srand(iflag)
      boxl=158.54D0
    
      reducedtemp=0.50
      setemp = reducedtemp*12

!      open(7,file='inputs/peptidex.inp',status='old')  
!      read(7,*) numbeads,boxl_orig,dave_rx
!      close(7)
!      open(7,file='inputs/peptidey.inp',status='old')  
!      read(7,*) numbeads,boxl_orig,dave_ry
!      close(7)
!      open(7,file='inputs/peptidez.inp',status='old')  
!      read(7,*) numbeads,boxl_orig,dave_rz
!      close(7)
7       format(A4,3X,I4,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3)
	  
!      open(7,file='ala31.pdb',status='unknown',form='formatted')
!      do i=1,31
!      write(7,7) 'ATOM',(i-1)*4+1,' N  ','ALA','A',i,dave_rx(i+31),dave_ry(i+31),dave_rz(i+31)
!      write(7,7) 'ATOM',(i-1)*4+2,' CA ','ALA','A',i,dave_rx(i),dave_ry(i),dave_rz(i)
!      write(7,7) 'ATOM',(i-1)*4+3,' C  ','ALA','A',i,dave_rx(i+2*31),dave_ry(i+2*31),dave_rz(i+2*31)
!      write(7,7) 'ATOM',(i-1)*4+4,' CB ','ALA','A',i,dave_rx(i+3*31),dave_ry(i+3*31),dave_rz(i+3*31)
!      end do
!      close(7)	  
	  
      open(7,file='abeta1-40.pdb',status='old')
      do i=1,chnln
      read(7,'(30x,3(f8.3))') old_rx(i+chnln),old_ry(i+chnln),old_rz(i+chnln)
      read(7,'(30x,3(f8.3))') old_rx(i),old_ry(i),old_rz(i)
      read(7,'(30x,3(f8.3))') old_rx(i+2*chnln),old_ry(i+2*chnln),old_rz(i+2*chnln)
      read(7,'(30x,3(f8.3))') old_rx(i+3*chnln),old_ry(i+3*chnln),old_rz(i+3*chnln)
      end do
      close(7)

	!read in identities
      open(7,file='inputs/identity.inp',status='old')   
      read (7,*) aa
      close(7)
      
      do l=1,nop,nb
       do k=1,nb 
		identity(l+k-1)=aa(k)
      enddo
      enddo

	!write file to verify identities
      open(7,file='checks/identity.out')
      do i=1,nop
		write(7,*) i,identity(i)
      enddo
      close(7)
	
	!this writes the hydrophobicity inputs
      open(7,file='parameters/hp1.inp')
      do i=1,nb
      if (identity(i) .le. 4) then
			write(7,*)'0'
      elseif (identity(i) .gt. 9) then
			write(7,*)'1'			
      endif
      enddo
      close(7)

	!this writes the sidechain inputs
      open(7,file='parameters/firstside1.data')
      do i=chnln*3+1,nb
      if (identity(i) .eq. 9) then
			write(7,*)'0'
      else
			write(7,*)'1'
      endif
      enddo
      close(7)

	!this writes the identity input for the simulation
      open(7,file='parameters/identity.inp')   
      do l=1,nb
      	if (identity(l) .ne. 9) write(7,*)identity(l)
      enddo
      close(7)

	del = 0.02375d0

	!read in l/d-isomer pseudobond constraints
      open(7,file='parameters/rcarnrco.data',status='old')
      do i=1,20
		read(7,"(6(f6.3,2x))") drca(i),drnh(i),drco(i),del_rca(i),del_rnh(i),del_rco(i)
	      if(del_rca(i).lt.del) del_rca(i)=del
      		if(del_rnh(i).lt.del) del_rnh(i)=del
      		if(del_rco(i).lt.del) del_rco(i)=del
      end do
      close(7)

	!reassign pseudobond identities
      do i=1,chnln
		iii=identity(3*chnln+i)-8
        	bdln(i)=drca(iii)
        	bl_rn(i)=drnh(iii)
        	bl_rc(i)=drco(iii)
        	del_bdln(i)=del_rca(iii)
        	del_blrn(i)=del_rnh(iii)
        	del_blrc(i)=del_rco(iii)
!!!!!!!!!!!!!
!         if(identity(i).eq.9) then
!	   move(i) = 0.7
!	      else
!        move(i) = 2.5
!		    endif
!!!!!!!!!!!!!!
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
	  move(1) = 1.5
	  move(2) = 0.5
	  move(3) = 0.6
	  move(4) = 1.5
	  move(5) = 2.0  
	  move(6) = 3.0  
  	  move(7) = 1.7  
  	  move(8) = 0.5	
  	  move(9) = 0.5		
  	  move(10) = 2.8
  	  move(11) = 0.7
  	  move(12) = 0.5
  	  move(13) = 3.0
  	  move(14) = 3.0
  	  move(15) = 1.0	
  	  move(15) = 1.0
  	  move(16) = 0.6
  	  move(17) = 0.5
  	  move(18) = 0.5
  	  move(19) = 1.5
  	  move(20) = 1.5
  	  move(21) = 0.5
  	  move(22) = 0.6
  	  move(23) = 2.0
  	  move(24) = 0.5
  	  move(25) = 0.6
  	  move(26) = 0.5
  	  move(27) = 0.5
  	  move(28) = 0.7
  	  move(29) = 0.7
  	  move(30) = 0.5
  	  move(31) = 0.5
  	  move(32) = 0.5
  	  move(33) = 0.6
  	  move(34) = 0.5
  	  move(35) = 1.5	
  	  move(36) = 0.5
  	  move(37) = 0.6
  	  move(38) = 0.6
  	  move(39) = 0.5
  	  move(40) = 0.5	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       OPEN(UNIT=7,FILE='parameters/sqz6to10.data',STATUS='OLD')
      DO I=1,20
       READ(7,'(5(f6.3,2x))') (SQZ610(K,I),K=1,5)
      END DO
       CLOSE(UNIT=7)
	   
	!this entire section takes each sidechain and moves the location of the bead in a box until it satisfies the pseudobond contraints
      do k=3*chnln1p+1,3*chnln1p+chnln
		rx2=(old_rx(k)-old_rx(k-3*chnln1p))
		ry2=(old_ry(k)-old_ry(k-3*chnln1p))
		rz2=(old_rz(k)-old_rz(k-3*chnln1p))
		scalerca=bdln(k-3*chnln1p)/dsqrt(rx2**2+ry2**2+rz2**2)
      		old_rx(k)=old_rx(k-3*chnln1p)+scalerca*rx2
      		old_ry(k)=old_ry(k-3*chnln1p)+scalerca*ry2
      		old_rz(k)=old_rz(k-3*chnln1p)+scalerca*rz2
      end do
                DS1 = 0.02
                DS2 = 0.02
                DS3 = 0.02
      write(6,*)
      write(6,'(a50)') "First value should be less than the 0.01 tolerance"
      write(6,'(a62)') "If not, adjust d1 d2 d3 values in sidechain adjustment section"

      do k=1,chnln
		rmin=100.d0		
         	do d1=-move(k),move(k),0.005d0
			rcax=old_rx(k+3*chnln1p)+d1-old_rx(k)
          		rnhx=old_rx(k+3*chnln1p)+d1-old_rx(k+chnln1p)
          		rcox=old_rx(k+3*chnln1p)+d1-old_rx(k+2*chnln1p)
             IF(K.GE.2) SZ6X=OLD_RX(K+3*CHNLN1P)+D1-OLD_RX(K-1+2*CHNLN1P)
             IF(K.GE.2) SZ8X=OLD_RX(K+3*CHNLN1P)+D1-OLD_RX(K-1)
             IF(K.LT.CHNLN) SZ7X=OLD_RX(K+3*CHNLN1P)+D1-OLD_RX(K+1+CHNLN1P)
             IF(K.LT.CHNLN) SZ9X=OLD_RX(K+3*CHNLN1P)+D1-OLD_RX(K+1)
         		do d2=-move(k),move(k),.005d0
           			rcay=old_ry(k+3*chnln1p)+d2-old_ry(k)
           			rnhy=old_ry(k+3*chnln1p)+d2-old_ry(k+chnln1p)
           			rcoy=old_ry(k+3*chnln1p)+d2-old_ry(k+2*chnln1p)
             IF(K.GE.2) SZ6Y=OLD_RY(K+3*CHNLN1P)+D2-OLD_RY(K-1+2*CHNLN1P)
             IF(K.GE.2) SZ8Y=OLD_RY(K+3*CHNLN1P)+D2-OLD_RY(K-1)
             IF(K.LT.CHNLN) SZ7Y=OLD_RY(K+3*CHNLN1P)+D2-OLD_RY(K+1+CHNLN1P)
             IF(K.LT.CHNLN) SZ9Y=OLD_RY(K+3*CHNLN1P)+D2-OLD_RY(K+1)
         			do d3=-move(k),move(k),.005d0
           				rcaz=old_rz(k+3*chnln1p)+d3-old_rz(k)
           				rnhz=old_rz(k+3*chnln1p)+d3-old_rz(k+chnln1p)
           				rcoz=old_rz(k+3*chnln1p)+d3-old_rz(k+2*chnln1p)
             IF(K.GE.2) SZ6Z=OLD_RZ(K+3*CHNLN1P)+D3-OLD_RZ(K-1+2*CHNLN1P)
             IF(K.GE.2) SZ8Z=OLD_RZ(K+3*CHNLN1P)+D3-OLD_RZ(K-1)
             IF(K.LT.CHNLN) SZ7Z=OLD_RZ(K+3*CHNLN1P)+D3-OLD_RZ(K+1+CHNLN1P)
             IF(K.LT.CHNLN) SZ9Z=OLD_RZ(K+3*CHNLN1P)+D3-OLD_RZ(K+1)
           				rca=dsqrt(rcax**2+rcay**2+rcaz**2)
           				rnh=dsqrt(rnhx**2+rnhy**2+rnhz**2)
           				rco=dsqrt(rcox**2+rcoy**2+rcoz**2)
             SZ6=DSQRT(SZ6X**2+SZ6Y**2+SZ6Z**2)
             SZ7=DSQRT(SZ7X**2+SZ7Y**2+SZ7Z**2)
             SZ8=DSQRT(SZ8X**2+SZ8Y**2+SZ8Z**2)
             SZ9=DSQRT(SZ9X**2+SZ9Y**2+SZ9Z**2)
!            SZ10=DSQRT(SZ10X**2+SZ10Y**2+SZ10Z**2)
                DS1 = DABS(BDLN(K)-RCA)
                DS2 = DABS(BL_RN(K)-RNH)
                DS3 = DABS(BL_RC(K)-RCO)
           				r3=dabs(bdln(k)-rca)+dabs(bl_rn(k)-rnh)+dabs(bl_rc(k)-rco)
						 ISZ=0
             IF(K.GE.2) THEN
             IF(SZ6.LT.SQZ610(2,IDENTITY(3*CHNLN+K)-8)) ISZ=1
             IF(SZ8.LT.SQZ610(1,IDENTITY(3*CHNLN+K)-8)) ISZ=1
             ENDIF
             IF(K.LT.CHNLN) THEN
             IF(SZ7.LT.SQZ610(3,IDENTITY(3*CHNLN+K)-8)) ISZ=1
             IF(SZ9.LT.SQZ610(4,IDENTITY(3*CHNLN+K)-8)) ISZ=1
             ENDIF
!            IF(K.GE.3) THEN
!            IF(SZ10.LT.SQZ610(5,IDENTITY(3*CHNLN+K)-8)) ISZ=1
!            ENDIF
             IF(ISZ.eq.0) ID=ID+1
             IF((R3.LT.RMIN).AND.(ISZ.EQ.0)) THEN
            					d1min=d1
          					    d2min=d2
               					d3min=d3
               					rmin=r3
           				endif
         			end do        
         		end do        
      end do        
		
       		old_rx(k+3*chnln1p)=old_rx(k+3*chnln1p)+d1min
       		old_ry(k+3*chnln1p)=old_ry(k+3*chnln1p)+d2min
       		old_rz(k+3*chnln1p)=old_rz(k+3*chnln1p)+d3min
       		rcax=old_rx(k+3*chnln1p)-old_rx(k)
       		rnhx=old_rx(k+3*chnln1p)-old_rx(k+chnln1p)
       		rcox=old_rx(k+3*chnln1p)-old_rx(k+2*chnln1p)
       		rcay=old_ry(k+3*chnln1p)-old_ry(k)
       		rnhy=old_ry(k+3*chnln1p)-old_ry(k+chnln1p)
       		rcoy=old_ry(k+3*chnln1p)-old_ry(k+2*chnln1p)
       		rcaz=old_rz(k+3*chnln1p)-old_rz(k)
       		rnhz=old_rz(k+3*chnln1p)-old_rz(k+chnln1p)
       		rcoz=old_rz(k+3*chnln1p)-old_rz(k+2*chnln1p)
       		rca=dsqrt(rcax**2+rcay**2+rcaz**2)
       		rnh=dsqrt(rnhx**2+rnhy**2+rnhz**2)
       		rco=dsqrt(rcox**2+rcoy**2+rcoz**2)
             IF(K.GE.2) SZ6X=old_RX(K+3*CHNLN)-old_RX(K-1+2*CHNLN)
             IF(K.GE.2) SZ8X=old_RX(K+3*CHNLN)-old_RX(K-1)
             IF(K.LT.CHNLN) SZ7X=old_RX(K+3*CHNLN)-old_RX(K+1+CHNLN)
             IF(K.LT.CHNLN) SZ9X=old_RX(K+3*CHNLN)-old_RX(K+1)
             IF(K.GE.3) SZ10X=old_RX(K+3*CHNLN)-old_RX(K-2+2*CHNLN)
             IF(K.GE.2) SZ6Y=old_RY(K+3*CHNLN)-old_RY(K-1+2*CHNLN)
             IF(K.GE.2) SZ8Y=old_RY(K+3*CHNLN)-old_RY(K-1)
             IF(K.LT.CHNLN) SZ7Y=old_RY(K+3*CHNLN)-old_RY(K+1+CHNLN)
             IF(K.LT.CHNLN) SZ9Y=old_RY(K+3*CHNLN)-old_RY(K+1)
             IF(K.GE.3) SZ10Y=old_RY(K+3*CHNLN)-old_RY(K-2+2*CHNLN)
             IF(K.GE.2) SZ6Z=old_RZ(K+3*CHNLN)-old_RZ(K-1+2*CHNLN)
             IF(K.GE.2) SZ8Z=old_RZ(K+3*CHNLN)-old_RZ(K-1)
             IF(K.LT.CHNLN) SZ7Z=old_RZ(K+3*CHNLN)-old_RZ(K+1+CHNLN)
             IF(K.LT.CHNLN) SZ9Z=old_RZ(K+3*CHNLN)-old_RZ(K+1)
             IF(K.GE.3) SZ10Z=old_RZ(K+3*CHNLN)-old_RZ(K-2+2*CHNLN)
             SZ6=DSQRT(SZ6X**2+SZ6Y**2+SZ6Z**2)
             SZ7=DSQRT(SZ7X**2+SZ7Y**2+SZ7Z**2)
             SZ8=DSQRT(SZ8X**2+SZ8Y**2+SZ8Z**2)
             SZ9=DSQRT(SZ9X**2+SZ9Y**2+SZ9Z**2)
             SZ10=DSQRT(SZ10X**2+SZ10Y**2+SZ10Z**2)			
            if (rmin .ge. 0.01) write(6,*)'ERROR in residue',k
        	write(6,'(i3,7(f8.4,1x))') k,rmin,d1min,d2min,d3min,rca,rnh,rco
      end do

      pi = 4.d0*datan(1.d0)
      coll = 0
      t=0.d0
  
	!this gets rid of glycine residues
        	l=0
    	 do i=0,3
    		do k=1,chnln
        		if(identity(i*chnln+k).ne.9) then
         			l=l+1
         			new_rx(l)=old_rx(i*chnln1p+k)
         			new_ry(l)=old_ry(i*chnln1p+k)
         			new_rz(l)=old_rz(i*chnln1p+k)
              endif
       		enddo
  	    enddo
		
      OPEN(7,file='chninfo-ab40.data',status='unknown',form='formatted')
      do i=1,nb2
	      write(7,*) new_rx(i),new_ry(i),new_rz(i)
      enddo
      CLOSE(UNIT=7)
	  
	!this makes a pdb to show the structure of one of each species
      open(7,file='checks/configone.pdb')
      write(7,*) nb2
      do i=1,nb2
       write(7,'(a6,i5,a3,1x,15x,3f8.3)') 'ATOM  ',i,'N', new_rx(i), new_ry(i), new_rz(i)
      enddo
      close(7)

	!this gets the cooridantes for the peptide relative to the first bead
      do i=1, nb2-1
      do j=i+1, nb2
		dislocx(j) = new_rx(j) - new_rx(i)
		dislocy(j) = new_ry(j) - new_ry(i)
		dislocz(j) = new_rz(j) - new_rz(i)
      enddo
      enddo
      do i=1,nb2
		xtmp(i)=new_rx(i)
		ytmp(i)=new_ry(i)
		ztmp(i)=new_rz(i)
      enddo

      	open(7,file='inputs/beadwell_ha55a.data',status='unknown')
      	do i=1,400
      		read(7,711) iiii,jjjj,bdtemp,wltemp
711  		format(2(2x,i2,2x),2(f6.3,2x))
       		bds(iiii,jjjj)=bdtemp
       		wel(iiii,jjjj)=wltemp
      	end do
      close(7)

      if(bds(29,29)/2 .gt. 2.000) then
		sig_max=2.000+bds(29,29)/2
      else
		sig_max=4.0D0
      endif
      write(6,*)'Maximum Sigma',sig_max
      do i=1,29
        do j=1,29
			sig_ij=bds(i,j)
	       	if (sig_ij .gt. sig_max) then
                  		sig_max = sig_ij
      write(6,*)'Maximum Sigma',sig_max
              	endif
      enddo
      enddo
		
	!this selects a random position, build the peptide in that location, and checks for overlaps
	!if there is an overlap it repeats
      do i=1,nc
		!write(6,*)'Peptide 1',i
993		xtmp((i-1)*nb2+1)=(drandm(0))*boxl
		xtmp((i-1)*nb2+1)=xtmp((i-1)*nb2+1)-boxl*ANINT(xtmp((i-1)*nb2+1)/boxl)
		ytmp((i-1)*nb2+1)=(drandm(0))*boxl
		ytmp((i-1)*nb2+1)=ytmp((i-1)*nb2+1)-boxl*ANINT(ytmp((i-1)*nb2+1)/boxl)
		ztmp((i-1)*nb2+1)=(drandm(0))*boxl
		ztmp((i-1)*nb2+1)=ztmp((i-1)*nb2+1)-boxl*ANINT(ztmp((i-1)*nb2+1)/boxl)
       do j=2,nb2
			xtmp((i-1)*nb2+j)=xtmp((i-1)*nb2+j-1)+dislocx(j)
			xtmp((i-1)*nb2+j)=xtmp((i-1)*nb2+j)-boxl*ANINT(xtmp((i-1)*nb2+j)/boxl)
			ytmp((i-1)*nb2+j)=ytmp((i-1)*nb2+j-1)+dislocy(j)
			ytmp((i-1)*nb2+j)=ytmp((i-1)*nb2+j)-boxl*ANINT(ytmp((i-1)*nb2+j)/boxl)
			ztmp((i-1)*nb2+j)=ztmp((i-1)*nb2+j-1)+dislocz(j)
			ztmp((i-1)*nb2+j)=ztmp((i-1)*nb2+j)-boxl*ANINT(ztmp((i-1)*nb2+j)/boxl)
       enddo
        do k=(i-1)*nb2+1,i*nb2
	      		do l=1,(i-1)*nb2
                       		rxij=xtmp(k)-xtmp(l)
                       		ryij=ytmp(k)-ytmp(l) 
                       		rzij=ztmp(k)-ztmp(l)
                      	 	rxij=rxij-boxl*aint((rxij/boxl)+.5)
                       		ryij=ryij-boxl*aint((ryij/boxl)+.5)
                       		rzij=rzij-boxl*aint((rzij/boxl)+.5)
	               		rijsq=rxij*rxij+ryij*ryij+rzij*rzij
	               		rijsq=rijsq*1.0000000001d0
                    if (rijsq .le. (5)**2) goto 993
           enddo
       enddo
      enddo

	!write full pdb file to verify configuration
      open(7,file='checks/config.pdb')
      do i=1,nop2
       write(7,'(a6,i5,a3,1x,15x,3f8.3)') 'ATOM  ',i,'N', xtmp(i), ytmp(i), ztmp(i)
      enddo
      close(7)

      open(7,file='checks/emblem.input')
	j=0
      do i=1,chnln
	j=j+1
		write(7,'(2i5,15x,3f8.3,i5)') j,j, xtmp(chnln+i), ytmp(chnln+i), ztmp(chnln+i),identity(chnln+i)
	j=j+1
		write(7,'(2i5,15x,3f8.3,i5)') j,j, xtmp(i), ytmp(i), ztmp(i), identity(i)
	j=j+1
		write(7,'(2i5,15x,3f8.3,i5)') j,j, xtmp(chnln*2+i), ytmp(chnln*2+i), ztmp(chnln*2+i), identity(chnln*2+i)
	j=j+1
		write(7,'(2i5,15x,3f8.3,i5)') j,j, xtmp(chnln*3+i), ytmp(chnln*3+i), ztmp(chnln*3+i), identity(chnln*3+i)
      enddo
         close(7)
   	
      do i=1,nop2
   		if(xtmp(i).lt.xmin) xmin=xtmp(i)
   		if(xtmp(i).gt.xmax) xmax=xtmp(i)
   		if(ytmp(i).lt.ymin) ymin=ytmp(i)
   		if(ytmp(i).gt.ymax) ymax=ytmp(i)
   		if(ztmp(i).lt.zmin) zmin=ztmp(i)
   		if(ztmp(i).gt.zmax) zmax=ztmp(i)
      end do
	
	!check is total bounds are within box length	
      write(6,*) 'These should be less than',boxl
      write(6,*) 'X-direction',abs(xmin)+abs(xmax)
      write(6,*) 'Y-direction',abs(ymin)+abs(ymax)
      write(6,*) 'Z-direction',abs(zmin)+abs(zmax)

	!write a config check
      open(7,file='checks/configcheck.out')
      write(7,*) coll,t
      write(7,*)xtmp/boxl,ytmp/boxl,ztmp/boxl
      close(7)	
	
	!this writes the actual config input
      open(7,file='results/run0000.config',status='unknown',form='unformatted')
        write(7) coll,t,xtmp/boxl,ytmp/boxl,ztmp/boxl
      close(7)
	
	!read in masses and adjust for bead identity
      open(7,file='parameters/mass.data',status='old')
      do i=1,23
       read(7,'(4x,i2,2x,f8.3)') iiii,bmass_temp
        bmass(iiii)=bmass_temp
      end do
	
	bmass(3)=bmass(20)
	
      do i=1,4
		bmass(i+4)=bmass(i)
      end do

      close(7)

      do i=1,nop
		bm(i)=bmass(identity(i))
      end do

	!write to make sure masses are correct
      open(7,file='checks/masses.out')
       do i=1,nop
        write(7,*) i,identity(i),bm(i)
      enddo
       close(7)

      write(6,*)'Desired reduced temperature',setemp/12

	!assign random velocities
      do k=1,nop
       if(identity(k).ne.9) then
         	x1=drandm(0)*0.75d0+0.2d0
         	x2=drandm(0)*0.75d0+0.2d0
         	y1=drandm(0)*0.75d0+0.2d0
         	y2=drandm(0)*0.75d0+0.2d0
         	sv(4,k)=dsqrt(-2.d0*log(x1))*cos(2.d0*pi*y1)
         	sv(5,k)=dsqrt(-2.d0*log(x1))*sin(2.d0*pi*y1)
         	sv(6,k)=dsqrt(-2.d0*log(x2))*cos(2.d0*pi*y2)
       endif
      end do

	!scale velocities so that linear momentum = 0
      do k=1,nop
       if(identity(k).ne.9) then
         	sumx=sumx+sv(4,k)
         	sumy=sumy+sv(5,k)
         	sumz=sumz+sv(6,k)
      		endif
      end do
	
      do k=1,nop
       if(identity(k).ne.9) then
         	sv(4,k)=sv(4,k)-sumx/real(nop2)
         	sv(5,k)=sv(5,k)-sumy/real(nop2)
         	sv(6,k)=sv(6,k)-sumz/real(nop2)
     	endif
      end do
	
	sumvel=0.d0
	
	!check summed velocities
      open(7,file='checks/sumvelcheck.out')
      do k=1,nop
       if(identity(k).ne.9) then
			sumvel=sumvel+(sv(4,k)*sv(4,k)+sv(5,k)*sv(5,k)+sv(6,k)*sv(6,k))/bm(k)
       endif
      enddo  
	!here sv is actually momentum not velocity
	 
	tred=sumvel/3.d0/dble(nop2)
	const=dsqrt(setemp/tred)

       do k=1,nop
       if(identity(k).ne.9) then
     		sv(4,k)=const*sv(4,k)/bm(k)
     		sv(5,k)=const*sv(5,k)/bm(k)
     		sv(6,k)=const*sv(6,k)/bm(k)
        endif
      enddo
  	
	sumvel=0.d0
      
      do k=1,nop
        if(identity(k).ne.9) then
			sumvel=sumvel+bm(k)*(sv(4,k)*sv(4,k)+sv(5,k)*sv(5,k)+sv(6,k)*sv(6,k))
			write(7,*)k,sumvel
      	endif
      enddo
      close(7)
	
	tred=sumvel/3.d0/dble(nop2)
      write(6,*)'New reduced temperature',tred/12.0

	!this is the velocity input to the simulations
      open(7,file='results/run0000.lastvel',status='unknown',form='unformatted') 
	l=0
      do k=1,nop   
         if(identity(k).ne.9) then
         	l=l+1
         	xtmp(l)=sv(4,k)
         	ytmp(l)=sv(5,k)
         	ztmp(l)=sv(6,k)
        endif
      enddo
      write(7) coll,xtmp,ytmp,ztmp
      close(7)

	!these are the two blank files needed to start
      open(7,file='results/run0000.bptnr',status='unknown',form='unformatted')
      close(7)
      open(7,file='results/run0000.energy',status='unknown',form='unformatted')
      close(7)

	!check velocities
      open(7,file='checks/velocitycheck.out')
	l=0
       do k=1,nop2   
       if(identity(k).ne.9) then
         	l=l+1
         	xtmp(l)=sv(4,k)
         	ytmp(l)=sv(5,k)
         	ztmp(l)=sv(6,k)
        write(7,*) k,identity(k),xtmp(k),ytmp(k),ztmp(k)
    		endif
      enddo
      close(7)

	!check HB parnters
      open(7,file='checks/hbcheck.out')
      do i=1,nop2
		write(7,*)i,bptnr(i)
      enddo
      close(7)

      stop
      end
