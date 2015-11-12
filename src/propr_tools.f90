subroutine assign_draws(probs,argsort,draws,output,nx,ny,nd,nc)

 implicit none
 integer,intent(in) :: nx,ny,nd,nc
 integer,intent(in) :: argsort(nc,nx,ny)
 real,intent(in) :: probs(nc,nx,ny),draws(nc,nd)
 real,intent(out) :: output(nd,nx,ny)
 integer :: i,j,k,l,samples(nc),tmp,args(nc)
 real :: dif
 output = 0.0
 samples = 0

 do i=1,nx
  do j=1,ny
   samples = int(ceiling(nd*probs(:,i,j)))
   tmp = sum(samples)
   if (tmp .gt. nd)then
    dif = tmp - nd
    do k=1,nc
     samples(k) = samples(k) - 1
     dif = dif - 1
     if (dif .eq. 0)exit
    enddo
   endif
   args = argsort(:,i,j) + 1 !Convert from python to fortran
   tmp = 1
   do k=1,nc
    if (samples(k) .eq. 0)exit !THIS COULD BE A PROBLEM IF SAMPLES(K) = 1
    do l=tmp,tmp+samples(k)
     output(l,i,j) = draws(args(k),l-tmp+1)
    enddo
    tmp = tmp + samples(k)
   enddo
  enddo 
 enddo

end subroutine
