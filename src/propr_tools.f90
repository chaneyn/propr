subroutine assign_draws(probs,argsort,draws,mapping,output,nx,ny,nd,nc,nc_all,nm)

 implicit none
 integer,intent(in) :: nx,ny,nd,nc,nc_all,nm
 integer,intent(in) :: argsort(nc,nx,ny),mapping(nm)
 real,intent(in) :: probs(nc,nx,ny),draws(nc_all,nd)
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
    if (samples(k) .le. 0)cycle
    do l=tmp,tmp+samples(k)-1
     output(l,i,j) = draws(mapping(args(k))+1,l-tmp+1)
    enddo
    tmp = tmp + samples(k)
   enddo
  enddo 
 enddo

end subroutine
