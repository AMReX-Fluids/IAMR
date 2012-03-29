
c ::: -----------------------------------------------------------
c ::: This routine will tag high error cells based on the 
c ::: density gradient
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: tag      <=  integer tag array
c ::: DIMS(tag) => index extent of tag array
c ::: set       => integer value to tag cell for refinement
c ::: clear     => integer value to untag cell
c ::: adv       => scalar array
c ::: DIMS(adv) => index extent of adv array
c ::: nvar      => number of components in rho array (should be 1)
c ::: lo,hi     => index extent of grid
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of tag array
c ::: problo    => phys loc of lower left corner of prob domain
c ::: time      => problem evolution time
c ::: -----------------------------------------------------------
       subroutine fort_adverror (tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3,
     &                           set,clear,
     &                           adv,advl1,advl2,advl3,advh1,advh2,advh3,
     &                           lo,hi,nvar,domlo,domhi,dx,xlo,
     &                           problo,time,level)
       implicit none

       integer          set, clear, nvar, level
       integer          tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
       integer          advl1,advl2,advl3,advh1,advh2,advh3
       integer          lo(3), hi(3), domlo(3), domhi(3)
       integer          tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
       double precision adv(advl1:advh1,advl2:advh2,advl3:advh3,nvar)
       double precision dx(3), xlo(3), problo(3), time

       integer          i,j,k
 
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                 tag(i,j,k) = merge(set,tag(i,j,k),adv(i,j,k,1).gt.0.5)
             end do
          end do
       end do
 
       end
