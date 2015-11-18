subroutine assign_draws(probs,argsort,draws,mapping,output,nx,ny,nd,nc,nc_all,nm)

 implicit none
 integer,intent(in) :: nx,ny,nd,nc,nc_all,nm
 integer,intent(in) :: argsort(nc,nx,ny),mapping(nm)
 real,intent(in) :: probs(nc,nx,ny),draws(nc_all,nd)
 real,intent(out) :: output(nd,nx,ny)
 integer :: i,j,k,l,samples(nc),tmp,args(nc)
 real :: dif,probs_cell(nc),undef
 output = 0.0
 samples = 0
 undef = -9999.0

 do i=1,nx
  do j=1,ny
   probs_cell = probs(:,i,j)
   !Clean up the probabilities (If we don't have an estimate then it makes no
   !sense to use it...)
   args = argsort(:,i,j) + 1
   do k=1,nc
    if (args(k) .lt. 0)cycle
    if (draws(mapping(args(k))+1,1) .eq. -9999)then
     probs_cell(k) = 0.0
    endif
   enddo
   !Set to undef if the sum of probabilities is now 0....
   if (sum(probs_cell) .eq. 0.0)then
     output(:,i,j) = undef
     cycle
   endif
   !Normalize probabilities 
   probs_cell = probs_cell/sum(probs_cell)
   !Define the samples
   samples = int(ceiling(nd*probs_cell))
   tmp = sum(samples)
   if (tmp .gt. nd)then
    dif = tmp - nd
    do k=1,nc
     if (samples(k).gt.0)then
      samples(k) = samples(k) - 1
      dif = dif - 1
     endif
     if (dif .eq. 0)exit
    enddo
   endif
   args = argsort(:,i,j) + 1 !Convert from python to fortran
   tmp = 1
   do k=1,nc
    if (samples(k) .le. 0)cycle
    do l=tmp,tmp+samples(k)-1
     output(l,i,j) = draws(mapping(args(k))+1,l-tmp+1)
    enddo
    tmp = tmp + samples(k)
   enddo
  enddo 
 enddo

end subroutine
