subroutine assign_draws(probs,argsort,draws,mapping,minimum,maximum,mean,&
                        alpha,beta,ncmax,nx,ny,nd,nc,nc_all,nm)

 implicit none
 integer,intent(in) :: nx,ny,nd,nc,nc_all,nm,ncmax
 integer,intent(in) :: argsort(nc,nx,ny),mapping(nm)
 real,intent(in) :: probs(nc,nx,ny),draws(nc_all,nd)
 !real,intent(out) :: output(nd,nx,ny)
 real,intent(out) :: minimum(nx,ny),maximum(nx,ny),mean(nx,ny),alpha(nx,ny),beta(nx,ny)
 integer :: i,j,k,l,samples(nc),tmp,args(nc),count
 real :: dif,probs_cell(nc),undef,array(nd),narray(nd),nmean,nvar
 real :: compute_min,compute_mean,compute_max,compute_var
 !output = 0.0
 samples = 0
 undef = -9999.0
 minimum = undef
 maximum = undef
 mean = undef
 alpha = undef
 beta = undef

 do i=1,nx
  do j=1,ny
   !Initialize array
   array = undef
   probs_cell = 0.0
   !Clean up the probabilities (If we don't have an estimate then it makes no
   !sense to use it...)
   args = argsort(:,i,j) + 1
   count = 0
   do k=1,nc
    probs_cell(k) = probs(k,i,j)
    if (args(k) .lt. 0)cycle
    if (draws(mapping(args(k))+1,1) .eq. -9999)then
     probs_cell(k) = 0.0
    else 
     count = count + 1
    endif
    if (count .eq. ncmax)exit
   enddo
   !Set to undef if the sum of probabilities is now 0....
   if (sum(probs_cell) .eq. 0.0)then
     !output(:,i,j) = undef
     array(:) = undef
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
     !output(l,i,j) = draws(mapping(args(k))+1,l-tmp+1)
     array(l) = draws(mapping(args(k))+1,l-tmp+1)
    enddo
    tmp = tmp + samples(k)
   enddo
   !Compute the metrics
   !array = undef
   !array = output(:,i,j)
   minimum(i,j) = compute_min(array,undef,nd)
   maximum(i,j) = compute_max(array,undef,nd) 
   mean(i,j) = compute_mean(array,undef,nd) 
   narray = (array - minimum(i,j))/(maximum(i,j) - minimum(i,j))
   nmean = compute_mean(narray,undef,nd)
   nvar = compute_var(narray,undef,nd)
   alpha(i,j) = ((1-nmean)/nvar - (1/nmean))*nmean**2
   beta(i,j) = alpha(i,j)*(1/nmean - 1)
  enddo 
 enddo

end subroutine

function compute_min(array,undef,nd) result(val)

 real,intent(in) :: array(nd),undef
 real :: val
 val = minval(array,mask=array.ne.undef)

end function compute_min

function compute_max(array,undef,nd) result(val)

 real,intent(in) :: array(nd),undef
 real :: val
 val = maxval(array,mask=array.ne.undef)

end function compute_max

function compute_mean(array,undef,nd) result(val)

 real,intent(in) :: array(nd),undef
 real :: val
 val = sum(array,mask=array.ne.undef)/max(1,count(array.ne.undef))

end function compute_mean

function compute_var(array,undef,nd) result(val)

 real,intent(in) :: array(nd),undef
 real :: val,mean,compute_mean
 mean = compute_mean(array,undef,nd)
 val = sum((array - mean)**2,mask=array.ne.undef)/max(1,count(array.ne.undef))

end function compute_var
