#include "def.h"
#include "header.f"
        program fibrilgrow	
use global
	  implicit none  
        REAL    GAUSS,RANDF,dummy,c(1)
		    external GAUSS
		
        integer a,b,d,e,a1,b1,i,j,ii,jj,iii,iiii,jjj,k,l,m,n,p,q,o,y,temp_count	
	
	integer izero,iflag,ichar, len, target_flag
	character*64 filename
	logical exist_flag
	character*1 zero,char
	data zero/'0'/
	izero = ichar('0')
	iflag = 9
	target_flag = 113
    target_fname_digits = char(target_flag/1000+izero)//char(mod(target_flag,1000)/100+izero)//char(mod(target_flag,100)/10+izero)//char(mod(target_flag,10)+izero)

	filename = '../run'//target_fname_digits//'.config'
	!fname_digits = '0155'
	inquire(file = filename, exist = exist_flag)
		
8       format(A4,3X,I4,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3)		
2 	    format(i15,3f12.4,3i8,4f12.4)	
			
		!boxl=329.770d0
			
		    open(11,file='../../parameters/identity.inp',action='read',form='formatted')	      
             read(11,*) aa
            close(11) 
        do i=1,nc
         do j=1,nb2,1
             k=(i-1)*nb2+j-1
             identity(k)=aa(j)
         enddo
        enddo
		
      	open(unit=7,file='../../parameters/beadwell_ha55a.data',status='unknown')
      	do i=1,(20)*(20)
      		read(7,711) ii,jj,bdtemp,wltemp
711  		format(2(2x,i2,2x),2(f6.3,2x))
       		bds(ii,jj)=bdtemp
       		wel(ii,jj)=wltemp
      	enddo
	      close(unit=7)			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				 do a =1,nop 
				 do b =1,nop
				 dis(a,b) = 0.0
				 enddo
				 enddo
			N_PHI_SHL(1) = 0.0
            N_PHI_SHL(2) = 0.0
			N_PSI_SHL(1) = 0.0
            N_PSI_SHL(2) = 0.0			
				  do a = 1,chnln*nc
				  beta(a) = 0
				  enddo					
			  	do a=1,nc
				move(a) = 0
				sheet_identity(a) = 0
				pep_num(a) = 0
				pep_count(a) = 0
				pep_nbr_num(a) = 0
				num_sheet(a) = 0
				fibril_identity(a) = 0
				sheet_count(a) = 0
				sheet_nbr_num(a) = 0
				  do b=1,3
				pep_nbr(a,b) = 0
				sheet_nbr(a,b) = 0
				  enddo
			    do b=1,nc
	          hb_beta(a,b) = 0
			   temp_pep_in_sheet(a,b) = 0
			   temp_sheet_in_fibril(a,b) = 0
		        hb_contact(a,b) = 0
				hp_contact(a,b) = 0
				dimer(a,b)  = 0
				hp_dimer(a,b) = 0 	
				  contact_count(a,b) = 0				
			   enddo 
			   enddo
			    
        count_file = 0 		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! recored final fibril information  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	     open(7,file='../run'//target_fname_digits//'.config',action='read',form='unformatted')				
         open(8,file='../run'//target_fname_digits//'.bptnr',form='unformatted',status='unknown')

	     do while (.true.)	
		    read(7,end=1120) coll,t,xtmp,ytmp,ztmp
		  !  count_file = count_file +1
			t=t*dsqrt(dble(setemp))/dble(sigma)
		     read(8,end=1121) coll1, bptnr
      	enddo
1120	   continue	
1121      continue   
           close(7)
		   close(8)	   

	      open(9,file='info'//target_fname_digits//'.dat',action='write',form='formatted')		  
		! find all dimer
		    do i =1,nc-1
		    do j = i+1,nc
                 do ii =(i-1)*nb+chnln+2,(i-1)*nb+chnln*2-1
	   	           if ((bptnr(ii).ne.0).and.((bptnr(ii)-mod(bptnr(ii),nb))/nb+1.eq.j)) then
			         hb_contact(i,j) = hb_contact(i,j) + 1	
			         hb_contact(j,i) = hb_contact(j,i) + 1	
			       endif
		          enddo	
				   
			   do ii =(i-1)*nb+chnln*2+2,(i-1)*nb+chnln*3-1
	   	          if ((bptnr(ii).ne.0).and.((bptnr(ii)-mod(bptnr(ii),nb))/nb+1.eq.j))  then 
			       hb_contact(i,j) = hb_contact(i,j) + 1	
			       hb_contact(j,i) = hb_contact(j,i) + 1				 
			      endif
		        enddo	

					 do ii = (i-1)*nb+3*chnln+1,(i-1)*nb+4*chnln	 
					 do jj = (j-1)*nb+3*chnln+1,(j-1)*nb+4*chnln	 
		    dx = xtmp(ii)-xtmp(jj)
		    dx = dx - dnint(dx)
            dy = ytmp(ii)-ytmp(jj)
		    dy = dy -dnint(dy)
            dz = ztmp(ii)-ztmp(jj)
	        dz = dz - dnint(dz)
            dis(ii,jj) = sqrt(dx**2+dy**2+dz**2)
			dis(jj,ii) = dis(ii,jj)
                if (dis(ii,jj).le.wel(identity(ii),identity(jj))/boxl) then 
			hp_contact(i,j) = hp_contact(i,j) + 1	
			hp_contact(j,i) = hp_contact(j,i) + 1	
			    endif
				
					 enddo
					 enddo
	
                if ((hp_contact(i,j).gt.0) .and.(hb_contact(i,j).eq.0)) then
			   hp_dimer(i,j) = 1
			   hp_dimer(j,i) = 1
			!   write(6,*) 'i am here hpdimer ',i,j
			    endif	

           if (hb_contact(i,j) .gt. (chnln)/2) then
              !call AddToAggregate_hb(i,j) 
		        dimer(i,j) = 1
				dimer(j,i) = 1	
              do m = (i-1)*nb+1, i*nb
                 do n = (j-1)*nb+1, j*nb
                    ii = m-(i-1)*nb
                    jj = n-(j-1)*nb
                    if((ii .ne. chnln+1).and.(ii .ne. 3*chnln).and.(jj .ne. chnln+1).and.(jj .ne. 3*chnln))then
                       if(bptnr(m) == n) then
                          hb_total=hb_total+1
                          call phipsi_beta(m,n,n_phi_shl,n_psi_shl)
                  
                          IF (((N_PHI_SHL(1) .le. -30).and.(N_PHI_SHL(1) .ge. -180)).AND. &
                              ((N_PHI_SHL(2) .le. -30).and.(N_PHI_SHL(2) .ge. -180)).AND. &
                              ((N_PSI_SHL(1) .le. 180).and.(N_PSI_SHL(1) .ge. 0)).AND. &
                              ((N_PSI_SHL(2) .le. 180).and.(N_PSI_SHL(2) .ge. 0))) THEN
        
                             HB_BETA(i,j)=HB_BETA(i,j) + 1
							 beta((i-1)*chnln+ii) = 1
							 beta((j-1)*chnln+jj) = 1
                           ENDIF

	               endif
                    endif  
                 enddo
              enddo   
           endif	
       		  
			     enddo	
			   enddo   	

		   do i = 1,nc2*chnln
         if	(beta(i).eq.1) beta_content = beta_content + 1.0
			enddo   
			beta_content = beta_content/nc2/(chnln-2)*100.0
			
        !identify neighbour peptides 1 or 2 of each peptide
            do a=1,nc-1
               do b=a+1,nc			
			     if (dimer(a,b).eq.1) then
			        pep_nbr_num(a) = pep_nbr_num(a) + 1
			        pep_nbr_num(b) = pep_nbr_num(b) + 1
					pep_nbr(a,pep_nbr_num(a)) = b
					pep_nbr(b,pep_nbr_num(b)) = a	
			    endif
		       enddo
		    enddo
        !assign each dimer to a sheet
			    sheet_num = 1		
		     do a=1,nc
		     if ((pep_nbr_num(a).ne.0).and.(sheet_identity(a).eq.0)) then	 
			    call sheet_assign(a)
			   do b=1,pep_num(sheet_num)
		      final_pep_in_sheet(sheet_num,b) = temp_pep_in_sheet(sheet_num,b)
                enddo	
           sheet_num = sheet_num + 1	                 		
			   endif   
		     enddo
	 
                 sheet_num = sheet_num -1  ! sheet_num change from temp value become final value
				 final_sheet_num = sheet_num	 
				 
              do i =1,nc
			  final_sheet_identity(i) = sheet_identity(i) 
       if (final_sheet_identity(i).ne.0) final_num_pep_in_sheet(final_sheet_identity(i)) = final_num_pep_in_sheet(final_sheet_identity(i))+1					
		        enddo				 
				 
				  ! assigning sheets to fibril
				  do i = 1,sheet_num-1
				  do j = i+1,sheet_num			  
				     do ii = 1,pep_num(i)
				     do jj = 1,pep_num(j)
					a = final_pep_in_sheet(i,ii)
					 b = final_pep_in_sheet(j,jj)
			      if (hp_dimer(a,b) .eq.1)  then 
				 contact_count(i,j) = contact_count(i,j) + 1
           !          write(6,*) 'i am here contact_count ',i,j,contact_count(i,j)	 
				    endif
					 enddo
					 enddo

					 if (contact_count(i,j).ge.2)  then
				sheet_nbr_num(i) = sheet_nbr_num(i) + 1
				  sheet_nbr_num(j) = sheet_nbr_num(j) + 1				  
				  sheet_nbr(i,sheet_nbr_num(i)) = j
				  sheet_nbr(j,sheet_nbr_num(j)) = i	
				  !write(6,*) i,j,contact_count(i,j),sheet_nbr_num(i),sheet_nbr_num(j)
				  endif
				  
					 enddo
                   enddo			   
				  
				   fibril_num = 1
				  
		     do a=1,sheet_num
		     if ((sheet_nbr_num(a).ne.0).and.(fibril_identity(a).eq.0)) then	 
			    call fibril_assign(a)
              
 			  do b=1,num_sheet(fibril_num)
		        final_sheet_in_fibril(fibril_num,b) = temp_sheet_in_fibril(fibril_num,b)
                enddo		 		 				
		         
                 fibril_num = fibril_num + 1		   
			   endif	
		     enddo
			 
				fibril_num = fibril_num -1  ! fibril_num change from temp value become final value
					 final_fibril_num = fibril_num
					 
				     do i = 1, final_fibril_num ! Yiming 
			     final_num_sheet_in_fibril(i) = num_sheet(i)  
			       enddo
				
              do i =1,sheet_num 
			  final_fibril_identity(i) = fibril_identity(i)
			 final_num_pep_in_fibril(final_fibril_identity(i)) = final_num_pep_in_fibril(final_fibril_identity(i))+pep_num(i)					
			   enddo
		
				   do e = 1,nc
					  if ((sheet_identity(e).ne.0).and.(fibril_identity(sheet_identity(e)).ne.0)) then
					    a = fibril_identity(sheet_identity(e))
					    move(a) = move(a) + 1
				        final_pep_in_fibril(a,move(a)) = e
				      endif 	   
                  enddo	

	     write(9,'(A4,F12.4,F8.3,10i6)') fname_digits,ttotal+t,beta_content,fibril_num,sheet_num,num_sheet(1),final_num_pep_in_fibril(1),num_sheet(2),final_num_pep_in_fibril(2),num_sheet(3),final_num_pep_in_fibril(3)		 
		 		 
         close(9)
		t= 0.0
		 coll = 0	 
           chnnum = 0
		   
		     do i =1,nc
			    if(sheet_identity(i).eq.0) then
				  do k=1,nb
        if ((bptnr((i-1)*nb+k).eq.0).or.((bptnr((i-1)*nb+k).ge.(i-1)*nb).and.(bptnr((i-1)*nb+k).le.i*nb))) y = y + 1
				   enddo
				  if (y.eq.nb) then 
			        chnnum = chnnum + 1
				!	if (chnnum .eq. 4) write(6,*) 'hha ',chnnum,i
		          if (chnnum .gt. 10) goto 777
					  
			       do j =1,nb
				   if (bptnr((i-1)*nb+j).eq.0) then
				     newbptnr((chnnum-1)*nb+j) = 0  
				   else
				   newbptnr((chnnum-1)*nb+j) = (chnnum-1)*nb+bptnr((i-1)*nb+j)-(i-1)*nb
				   endif
			   newxtmp((chnnum-1)*nb+j) = xtmp((i-1)*nb+j)
			   newytmp((chnnum-1)*nb+j) = ytmp((i-1)*nb+j)
			   newztmp((chnnum-1)*nb+j) = ztmp((i-1)*nb+j)		
          !  if (chnnum.eq.4) write(6,*) boxl*newxtmp((4-1)*nb+j),boxl*newytmp((4-1)*nb+j),boxl*newztmp((4-1)*nb+j)		   
			    enddo
				  endif
				  y = 0
				  endif				
		     enddo
777     chnnum = chnnum -1
write(6,*) 'solubility: ',dble(chnnum)*10**7/6.023/(boxl**3)

			 y = chnnum
		     do i =1,nc			 
               if (fibril_identity(sheet_identity(i)).eq.1) then
			     chnnum = chnnum + 1
			    do j =1,nb	 
				   
			        if (pep_nbr_num(i).eq.1) then
					
						do k=1,final_num_pep_in_fibril(1)
			    if (final_pep_in_fibril(1,k).eq.pep_nbr(i,1)) then 
					 if (bptnr((i-1)*nb+j).eq.0) then
				     newbptnr((chnnum-1)*nb+j) = 0  
				   else
		      newbptnr((chnnum-1)*nb+j) = (k-1)*nb+bptnr((i-1)*nb+j)-(pep_nbr(i,1)-1)*nb
				   endif
				  endif 
				       enddo
					
                   elseif (pep_nbr_num(i).eq.2) then					

						do k=1,final_num_pep_in_fibril(1)
			    if (final_pep_in_fibril(1,k).eq.pep_nbr(i,1)) then 
					 if (bptnr((i-1)*nb+j).eq.0) then
				     newbptnr((chnnum-1)*nb+j) = 0  
				   elseif (mod(bptnr((i-1)*nb+j),nb)+1.eq. pep_nbr(i,1)) then 
		             newbptnr((chnnum-1)*nb+j) = (k-1)*nb+bptnr((i-1)*nb+j)-(pep_nbr(i,1)-1)*nb
				   endif				
				
			   elseif(final_pep_in_fibril(1,k).eq.pep_nbr(i,2)) then
					 if (bptnr((i-1)*nb+j).eq.0) then
				     newbptnr((chnnum-1)*nb+j) = 0  
				   elseif (mod(bptnr((i-1)*nb+j),nb)+1.eq. pep_nbr(i,2)) then 
		             newbptnr((chnnum-1)*nb+j) = (k-1)*nb+bptnr((i-1)*nb+j)-(pep_nbr(i,2)-1)*nb
				   endif			 
		
            		endif 
				       enddo
					   
				   endif				 
			   newxtmp((chnnum-1)*nb+j) = xtmp((i-1)*nb+j)
			   newytmp((chnnum-1)*nb+j) = ytmp((i-1)*nb+j)
			   newztmp((chnnum-1)*nb+j) = ztmp((i-1)*nb+j)						  
			    enddo   
			   endif			 
	         enddo		 
		write(6,*) 'free peptide num, peptide in fibril: ',y,chnnum-y
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! recored final config  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	     open(10,file='run0000.config',action='write',form='unformatted')				
		 open(11,file='run0000.bptnr',action='write',form='unformatted')   
	  nop3 = nb*chnnum
	  nop4 = nop3	
		    allocate(newbptnr2(nop3)) 		
		    allocate(newxtmp2(nop3))
			allocate(newytmp2(nop3))
			allocate(newztmp2(nop3))
			allocate(vx2(nop3))
			allocate(vy2(nop3))
		    allocate(vz2(nop3))		
            allocate(bm(nop3))			
		 
		     do i = 1,nop3
			 newbptnr2(i) = newbptnr(i) 
			 newxtmp2(i) = newxtmp(i)
			 newytmp2(i) = newytmp(i)
			 newztmp2(i) = newztmp(i)			 
		     enddo
			 
		    write(10) coll,t,newxtmp2,newytmp2,newztmp2
		     write(11) coll, newbptnr2	
			    close(10)
		      close(11)	 

	open(18,file='config.pdb',action='write',form='formatted')
	write(18,*) nop3
	do i=1,nop3
		write(18,'(a6,i5,a3,1x,15x,3f8.3)') 'ATOM  ',i,'N', boxl*newxtmp2(i), boxl*newytmp2(i), boxl*newztmp2(i)
	enddo
	close(18)
			  
!****************************************
      !!// velocities

      open(unit=7,file='../../parameters/mass.data',status='old')
       do i=1,23
       read(7,'(4x,i2,2x,F8.3)') iiii,bmass_temp
        bmass(iiii)=bmass_temp
       end do
        bmass(3)=bmass(20)
       do i=1,4
        bmass(i+4)=bmass(i)
       end do
      close(unit=7)

       do i=1,nop3
        bm(i)=bmass(identity(i))
       end do

      write(6,*)'desired temperature = ',setemp

      sumx=0.0
      sumy=0.0
      sumz=0.0
      do k=1,nop3
         if(identity(k).ne.9) then
        ! x1=drandm(0)*0.75d0+0.2d0
        ! x2=drandm(0)*0.75d0+0.2d0
        ! y1=drandm(0)*0.75d0+0.2d0
        ! y2=drandm(0)*0.75d0+0.2d0
        ! sv(4,k)=dsqrt(-2.d0*log(x1))*cos(2.d0*Pi*y1)
        ! sv(5,k)=dsqrt(-2.d0*log(x1))*sin(2.d0*Pi*y1)
        ! sv(6,k)=dsqrt(-2.d0*log(x2))*cos(2.d0*Pi*y2)
		  vx(k)=gauss(dummy)
          vy(k)=gauss(dummy)
          vz(k)=gauss(dummy)
         endif
      end do

!     scale velocities so that linear momentum = 0
      do k=1,nop3
         if(identity(k).ne.9) then
         sumx=sumx+vx(k)
         sumy=sumy+vy(k)
         sumz=sumz+vz(k)
         endif
      end do  
	  
      do k=1,nop3
         if(identity(k).ne.9) then
         vx(k)=vx(k)-sumx/dble(nop4)
         vy(k)=vy(k)-sumy/dble(nop4)
         vz(k)=vz(k)-sumz/dble(nop4)
         endif
      end do	  
	  
      sumvel=0.d0
      do k=1,nop3
         if(identity(k).ne.9) then
         sumvel=sumvel+(vx(k)*vx(k)+vy(k)*vy(k)+vz(k)*vz(k))/bm(k)
!        write(6,*) k,sumvel,bm(k)
         endif
      enddo  ! since Here sv is actually momentum not velocity (By Mookyung)   
	  
      tred=sumvel/3.d0/dble(nop4)

      const=dsqrt(setemp/tred)

      do k=1,nop3
         if(identity(k).ne.9) then
         vx(k)=const*vz(k)/bm(k)
         vy(k)=const*vy(k)/bm(k)
         vz(k)=const*vz(k)/bm(k)
         endif
      enddo
	  
      sumvel=0.d0
      do k=1,nop3
         if(identity(k).ne.9) then
         sumvel=sumvel+bm(k)*(vx(k)*vx(k)+vy(k)*vy(k)+vz(k)*vz(k))
         endif
      enddo 
	  
      tred=sumvel/3.d0/dble(nop4)
      write(6,*)'new scaled temperature = ',tred


      open(unit=13,file='run0000.lastvel',action='write',form='unformatted')
      l=0
      do k=1,nop3   
         if(identity(k).ne.9) then
         l=l+1
         vx2(l)=vx(k)
         vy2(l)=vy(k)
         vz2(l)=vz(k)
         endif
      enddo 

      write(13) coll,vx2,vy2,vz2
	  
      stop		 
		 
		 
	     close(7)
           close(8) 
		 close(12)
		 close(13)
           
		   deallocate(newbptnr2) 		
		deallocate(newxtmp2)
		deallocate(newytmp2)
		deallocate(newztmp2)
		deallocate(vx2)
		deallocate(vy2)
		deallocate(vz2)		
        deallocate(bm)		
		    end
			
       recursive subroutine sheet_assign(i)
#include "def.h"	   
	   use global
	      integer a,b,i,j,k,ii
		  
			sheet_identity(i) = sheet_num
			pep_count(i) = pep_count(i) + 1

			   if (pep_count(i).eq.1)  then
			   pep_num(sheet_num) = pep_num(sheet_num) + 1
			temp_pep_in_sheet(sheet_num,pep_num(sheet_num)) = i
                   endif	
			  			   
			    if(pep_nbr_num(i).eq.1) then				
				  a = pep_nbr(i,1)

					 if(sheet_identity(a).eq.0) then
			sheet_identity(a) = sheet_num
			pep_count(a) = pep_count(a) + 1
			   if (pep_count(a).eq.1) then
			   pep_num(sheet_num) = pep_num(sheet_num) + 1
			temp_pep_in_sheet(sheet_num,pep_num(sheet_num)) = a 
              endif
			         endif

					  if(pep_nbr_num(a).gt.1) then
					     b = pep_nbr(a,2)
                         if (sheet_identity(b).eq.0) call sheet_assign(b)			   
					 endif
				  
	      elseif (pep_nbr_num(i).gt.1) then
				     do ii=1,pep_nbr_num(i)
					  a = pep_nbr(i,ii)
			           if (sheet_identity(a).eq.0) call sheet_assign(a)
                      enddo			   
                 endif			
         
			
688 		   return   

	     end
		 
		    recursive  subroutine fibril_assign(i) 
#include "def.h"	
 use global
	      integer a,b,i,j,k,ii
		  
			fibril_identity(i) = fibril_num
			sheet_count(i) = sheet_count(i) + 1

			   if (sheet_count(i).eq.1) then
			   num_sheet(fibril_num) = num_sheet(fibril_num) + 1
			   temp_sheet_in_fibril(fibril_num,num_sheet(fibril_num)) = i 
              endif
				   
			    if(sheet_nbr_num(i).eq.1) then				
				    a = sheet_nbr(i,1)

					 if(fibril_identity(a).eq.0) then
					      fibril_identity(a) = fibril_num
					      sheet_count(a) = sheet_count(a) + 1
					    if (sheet_count(a).eq.1) then
					     num_sheet(fibril_num) = num_sheet(fibril_num) + 1
					     temp_sheet_in_fibril(fibril_num,num_sheet(fibril_num)) = a 
				        endif
			            endif	  				 
                
				     if(sheet_nbr_num(a).gt.1) then			 
					      b = sheet_nbr(a,2)
                       if (fibril_identity(b).eq.0) call fibril_assign(b)			   
					   endif
				  
				  elseif (sheet_nbr_num(i).gt.1) then
				     do ii=1,sheet_nbr_num(i)
					  a = sheet_nbr(i,ii)
			           if (fibril_identity(a).eq.0) call fibril_assign(a)
                      enddo		
				   		   
                 endif						  
			
687				return
               
		     end

!    *******************************************************************
!    ** TRANSLATIONAL VELOCITIES FROM MAXWELL-BOLTZMANN DISTRIBUTION  **
!    **                                                               **
!    ** THE DISTRIBUTION IS DETERMINED BY TEMPERATURE AND (UNIT) MASS.**
!    ** THIS ROUTINE IS GENERAL, AND CAN BE USED FOR ATOMS, LINEAR    **
!    ** MOLECULES, AND NON-LINEAR MOLECULES.                          **
!    **                                                               **
!    ** ROUTINE REFERENCED:                                           **
!    **                                                               **
!    ** REAL FUNCTION GAUSS ( DUMMY )                                 **
!    **    RETURNS A UNIFORM RANDOM NORMAL VARIATE FROM A             **
!    **    DISTRIBUTION WITH ZERO MEAN AND UNIT VARIANCE.             **
!    *******************************************************************

        
   
     REAL FUNCTION GAUSS (dummy)                                           
                                                                                
!    *******************************************************************
!    ** RANDOM VARIATE FROM THE STANDARD NORMAL DISTRIBUTION.         **
!    **                                                               **
!    ** THE DISTRIBUTION IS GAUSSIAN WITH ZERO MEAN AND UNIT VARIANCE.**
!    **                                                               **
!    ** REFERENCE:                                                    **
!    **                                                               **
!    ** KNUTH D, THE ART OF COMPUTER PROGRAMMING, (2ND EDITION        **
!    **    ADDISON-WESLEY), 1978                                      **
!    **                                                               **
!    ** ROUTINE REFERENCED:                                           **
!    **                                                               **
!    ** REAL FUNCTION RANF ( DUMMY )                                  **
!    **    RETURNS A UNIFORM RANDOM VARIATE ON THE RANGE ZERO TO ONE  **
!    *******************************************************************
                                                             
        REAL        A1, A3, A5, A7, A9   		
        PARAMETER ( A1 = 3.949846138, A3 = 0.252408784 )                
        PARAMETER ( A5 = 0.076542912, A7 = 0.008355968 )                
        PARAMETER ( A9 = 0.029899776                   )                
                                                                        
        REAL        SUM, R, R2                                          
        REAL        RANDF, DUMMY
        EXTERNAL    RANDF		
        INTEGER     I,c(1)                                                   
                                                                        
!    *******************************************************************
			
        SUM = 0.0                                                       
         i = 1                                                                     		
            c(1)=i**2
           i = i + 1			
           SUM = SUM + randf(c(1))                                   
            c(1)=i**2
           i = i + 1                                        
           SUM = SUM + randf(c(1)) 
            c(1)=i**2
           i = i + 1                                        
           SUM = SUM + randf(c(1)) 
            c(1)=i**2
           i = i + 1                                        
           SUM = SUM + randf(c(1)) 
            c(1)=i**2
           i = i + 1                                       
           SUM = SUM + randf(c(1)) 
            c(1)=i**2
           i = i + 1                                        
           SUM = SUM + randf(c(1)) 
            c(1)=i**2
           i = i + 1                                       
           SUM = SUM + randf(c(1)) 
            c(1)=i**2
           i = i + 1                                        
           SUM = SUM + randf(c(1)) 
            c(1)=i**2
           i = i + 1                                       
           SUM = SUM + randf(c(1)) 
            c(1)=i**2
           i = i + 1                                       
           SUM = SUM + randf(c(1)) 
            c(1)=i**2
           i = i + 1                                       
           SUM = SUM + randf(c(1)) 
            c(1)=i**2
           i = 0                                      
           SUM = SUM + randf(c(1)) 		   
        R  = ( SUM - 6.0 ) / 4.0                                        
        R2 = R * R                                                      
                                                                        
        GAUSS = (((( A9 * R2 + A7 ) * R2 + A5 ) * R2 + A3 ) * R2 +A1 ) * R                                                     
                                                                        
        RETURN                                                          
        END                                                             
                                                                        
 REAL FUNCTION RANDF ( c )                                    
                                                                        
!    *******************************************************************
!    ** RETURNS A UNIFORM RANDOM VARIATE IN THE RANGE 0 TO 1.         **
!    **                                                               **
!    **                 ***************                               **
!    **                 **  WARNING  **                               **
!    **                 ***************                               **
!    **                                                               **
!    ** GOOD RANDOM NUMBER GENERATORS ARE MACHINE SPECIFIC.           **
!    ** PLEASE USE THE ONE RECOMMENDED FOR YOUR MACHINE.              **
!    *******************************************************************
        integer,intent(in) ::c(1)
		integer b(1)
        real k  				  
 		integer factor
	    integer                            :: seedSize
		integer, dimension(:), allocatable :: seed
		integer, dimension(8)              :: dtVals

  seed=0
  dtVals=0

	 factor = 17271   ! this factor amplify the diversity of random number            
  call DATE_AND_TIME(VALUES=dtVals)
  call RANDOM_SEED(SIZE=seedSize)
      if(seedSize .gt. 8) then
     write (*,*) 'ERROR: Seed size too large to init with DATE_AND_TIME return'
     stop
  end if
       allocate(seed(seedSize)) 
  call RANDOM_SEED(GET=seed)
  b(1)=dtVals(8)
  call RANDOM_SEED(PUT=factor*(c(:)+b(:)))
  call RANDOM_SEED(GET=seed)
  CALL RANDOM_NUMBER(k)	                
	 
			 randf=k
	  
        RETURN        
      deallocate(seed) 		
        END 	  			 
			 
#include "phipsi_beta.f"		 