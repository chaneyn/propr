subroutine assign_draws(probs,argsort,draws,mapping,max_in,min_in,minimum,maximum,mean,&
                        alpha,beta,var,ncmax,nx,ny,nd,nc,nc_all,nm)

 implicit none
 integer,intent(in) :: nx,ny,nd,nc,nc_all,nm,ncmax
 integer,intent(in) :: argsort(nc,nx,ny),mapping(nm)
 real,intent(in) :: probs(nc,nx,ny),draws(nc_all,nd),min_in(nc_all),max_in(nc_all)
 !real,intent(out) :: output(nd,nx,ny)
 real,intent(out) :: minimum(nx,ny),maximum(nx,ny),mean(nx,ny),alpha(nx,ny),beta(nx,ny),var(nx,ny)
 integer :: i,j,k,l,samples(nc),tmp,args(nc),count
 real :: dif,probs_cell(nc),undef,array(nd),narray(nd),nmean,nvar,min_cell(nc),max_cell(nc)
 real :: compute_min,compute_mean,compute_max,compute_var
 !output = 0.0
 samples = 0
 undef = -9999.0
 minimum = undef
 maximum = undef
 mean = undef
 alpha = undef
 beta = undef
 var = undef

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
     !array(:) = undef
     cycle
   endif
   !Normalize probabilities 
   probs_cell = probs_cell/sum(probs_cell)
   !Define the samples
   samples = int(ceiling(nd*probs_cell))
   !print*,probs_cell
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
   !print*,samples
   !stop
   !args = argsort(:,i,j) + 1 !Convert from python to fortran
   tmp = 1
   do k=1,nc
    if (samples(k) .le. 0)cycle
    do l=tmp,tmp+samples(k)-1
     !output(l,i,j) = draws(mapping(args(k))+1,l-tmp+1)
     array(l) = draws(mapping(args(k))+1,l-tmp+1)
    enddo
    tmp = tmp + samples(k)
   enddo
   !HACK
   min_cell = undef
   max_cell = undef
   do k=1,nc
    if (args(k) .lt. 0)cycle
    min_cell(k) = min_in(mapping(args(k))+1)
    max_cell(k) = max_in(mapping(args(k))+1)
   enddo
   !Compute the parameters
   minimum(i,j) = minval(min_cell,mask=probs_cell.ne.0)
   maximum(i,j) = maxval(max_cell,mask=probs_cell.ne.0)
   !HACK
   !Compute the metrics
   !array = undef
   !array = output(:,i,j)
   !minimum(i,j) = compute_min(array,undef,nd)
   !maximum(i,j) = compute_max(array,undef,nd) 
   !if ((i.eq.1) .and. (j .eq. 1))then
   ! print*,array
   ! print*,probs_cell
   ! print*,maximum(i,j)
   !endif
   !stop
   mean(i,j) = compute_mean(array,undef,nd) 
   var(i,j) = compute_var(array,undef,nd)
   narray = (array - minimum(i,j))/(maximum(i,j) - minimum(i,j))
   nmean = compute_mean(narray,undef,nd)
   nvar = compute_var(narray,undef,nd)
   alpha(i,j) = ((1-nmean)/nvar - (1/nmean))*nmean**2
   beta(i,j) = alpha(i,j)*(1/nmean - 1)
   if ((i.eq.20) .and. (j .eq. 10))then
    !print*,narray
    !print*,mean(i,j),alpha(i,j),beta(i,j),nmean,nvar,minimum(i,j),maximum(i,j)
    !print*,(mean(i,j) - minimum(i,j))/(maximum(i,j) - minimum(i,j))
   endif
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
 val = sum((array - mean)**2,mask=array.ne.undef)/max(1,count(array.ne.undef)-1)

end function compute_var

subroutine compute_parameters(max_in,min_in,mean_in,var_in,&
                        probs,argsort,mapping,minimum,maximum,mean,var,&
                        alpha,beta,ncmax,nx,ny,nc,nc_all,nm)

 implicit none
 integer,intent(in) :: nx,ny,nc,nc_all,nm,ncmax
 integer,intent(in) :: argsort(nc,nx,ny),mapping(nm)
 real,intent(in) :: probs(nc,nx,ny),max_in(nc_all),min_in(nc_all),mean_in(nc_all),var_in(nc_all)
 real,intent(out) :: minimum(nx,ny),maximum(nx,ny),mean(nx,ny),alpha(nx,ny),beta(nx,ny),var(nx,ny)
 integer :: i,j,k,l,samples(nc),args(nc),count
 real :: probs_cell(nc),undef
 real :: var_cell(nc),mean_cell(nc),min_cell(nc),max_cell(nc)
 real :: nmean,nvar
 undef = -9999.0
 minimum = undef
 maximum = undef
 mean = undef
 alpha = undef
 beta = undef
 var = undef

 do i=1,nx
  do j=1,ny
   !Initialize array
   probs_cell = 0.0
   !Clean up the probabilities (If we don't have an estimate then it makes no
   !sense to use it...)
   args = argsort(:,i,j) + 1
   count = 0
   do k=1,nc
    if (args(k) .lt. 0)cycle
    if (mean_in(mapping(args(k))+1) .eq. -9999)then
     probs_cell(k) = 0.0
    else 
     probs_cell(k) = probs(k,i,j)
     count = count + 1
    endif
    if (count .eq. ncmax)exit
   enddo
   !Set to undef if the sum of probabilities is now 0....
   if (sum(probs_cell) .eq. 0.0)then
     maximum(i,j) = undef
     minimum(i,j) = undef
     var(i,j) = undef
     mean(i,j) = undef
     alpha(i,j) = undef
     beta(i,j) = undef
   else
     !Assemble the parameters
     min_cell = undef
     max_cell = undef 
     mean_cell = undef
     var_cell = undef
     do k=1,nc
      if (args(k) .lt. 0)cycle
      min_cell(k) = min_in(mapping(args(k))+1)
      max_cell(k) = max_in(mapping(args(k))+1)
      mean_cell(k) = mean_in(mapping(args(k))+1)
      var_cell(k) = var_in(mapping(args(k))+1)
     enddo
     probs_cell = probs_cell/sum(probs_cell)
     !Compute the parameters
     minimum(i,j) = minval(min_cell,mask=probs_cell.ne.0)
     maximum(i,j) = maxval(max_cell,mask=probs_cell.ne.0)
     mean(i,j) = sum(probs_cell*mean_cell)
     var(i,j) = sum(probs_cell*(var_cell + mean(i,j)**2)) - mean(i,j)**2
     !Normalize the mean and variance
     nmean = (mean(i,j) - minimum(i,j))/(maximum(i,j) - minimum(i,j))
     nvar = 1/(maximum(i,j) - minimum(i,j))**2*(var(i,j) + mean(i,j)**2 &
            - 2*minimum(i,j)*mean(i,j) + minimum(i,j)**2) - nmean**2
     alpha(i,j) = ((1-nmean)/nvar - (1/nmean))*nmean**2
     beta(i,j) = alpha(i,j)*(1/nmean - 1)
   endif
  enddo 
 enddo

end subroutine

subroutine compute_parameters_point(max_in,min_in,mean_in,var_in,&
                        probs,argsort,mapping,minimum,maximum,mean,var,&
                        alpha,beta,ncmax,np,nc,nc_all,nm)

 implicit none
 integer,intent(in) :: np,nc,nc_all,nm,ncmax
 integer,intent(in) :: argsort(nc,np),mapping(nm)
 real,intent(in) :: probs(nc,np),max_in(nc_all),min_in(nc_all),mean_in(nc_all),var_in(nc_all)
 real,intent(out) :: minimum(np),maximum(np),mean(np),alpha(np),beta(np),var(np)
 integer :: i,j,k,l,samples(nc),args(nc),count
 real :: probs_cell(nc),undef
 real :: var_cell(nc),mean_cell(nc),min_cell(nc),max_cell(nc)
 real :: nmean,nvar
 undef = -9999.0
 minimum = undef
 maximum = undef
 mean = undef
 alpha = undef
 beta = undef
 var = undef

 do i=1,np
   !Initialize array
   probs_cell = 0.0
   !Clean up the probabilities (If we don't have an estimate then it makes no
   !sense to use it...)
   args = argsort(:,i) + 1
   count = 0
   do k=1,nc
    if (args(k) .lt. 0)cycle
    if (mean_in(mapping(args(k))+1) .eq. -9999)then
     probs_cell(k) = 0.0
    else 
     probs_cell(k) = probs(k,i)
     count = count + 1
    endif
    if (count .eq. ncmax)exit
   enddo
   !Set to undef if the sum of probabilities is now 0....
   if (sum(probs_cell) .eq. 0.0)then
     maximum(i) = undef
     minimum(i) = undef
     var(i) = undef
     mean(i) = undef
     alpha(i) = undef
     beta(i) = undef
   else
     !Assemble the parameters
     min_cell = undef
     max_cell = undef 
     mean_cell = undef
     var_cell = undef
     do k=1,nc
      if (args(k) .lt. 0)cycle
      min_cell(k) = min_in(mapping(args(k))+1)
      max_cell(k) = max_in(mapping(args(k))+1)
      mean_cell(k) = mean_in(mapping(args(k))+1)
      var_cell(k) = var_in(mapping(args(k))+1)
     enddo
     probs_cell = probs_cell/sum(probs_cell)
     !Compute the parameters
     minimum(i) = minval(min_cell,mask=probs_cell.ne.0)
     maximum(i) = maxval(max_cell,mask=probs_cell.ne.0)
     mean(i) = sum(probs_cell*mean_cell)
     var(i) = sum(probs_cell*(var_cell + mean_cell**2)) - mean(i)**2
     !Normalize the mean and variance
     nmean = (mean(i) - minimum(i))/(maximum(i) - minimum(i))
     nvar = 1/(maximum(i) - minimum(i))**2*(var(i) + mean(i)**2 &
            - 2*minimum(i)*mean(i) + minimum(i)**2) - nmean**2
     alpha(i) = ((1-nmean)/nvar - (1/nmean))*nmean**2
     beta(i) = alpha(i)*(1/nmean - 1)
   endif
 enddo

end subroutine

subroutine compute_weighted_histogram(hist,probs,argsort,mapping,whist,ncmax,nx,ny,nc,nc_all,nm,nb)

 implicit none
 integer,intent(in) :: nx,ny,nc,nc_all,nm,ncmax,nb
 integer,intent(in) :: argsort(nc,nx,ny),mapping(nm)
 real,intent(in) :: probs(nc,nx,ny),hist(nc_all,nb)
 real,intent(out) :: whist(nx,ny,nb)
 integer :: i,j,k,l,args(nc),count
 real :: probs_cell(nc),undef
 real :: hist_cell(nc,nb)
 real :: nmean,nvar
 undef = -9999.0
 whist = undef

 do i=1,nx
  do j=1,ny
   !Initialize array
   probs_cell = 0.0
   !Clean up the probabilities (If we don't have an estimate then it makes no
   !sense to use it...)
   args = argsort(:,i,j) + 1
   count = 0
   do k=1,nc
    if (args(k) .lt. 0)cycle
    if (hist(mapping(args(k))+1,1) .eq. -9999)then
     probs_cell(k) = 0.0
    else 
     probs_cell(k) = probs(k,i,j)
     count = count + 1
    endif
    if (count .eq. ncmax)exit
   enddo
   !Set to undef if the sum of probabilities is now 0....
   if (sum(probs_cell) .eq. 0.0)then
     whist(i,j,:) = undef
   else
     !Assemble the parameters
     hist_cell = undef
     do k=1,nc
      if (args(k) .lt. 0)cycle
      hist_cell(k,:) = hist(mapping(args(k))+1,:)
     enddo
     probs_cell = probs_cell/sum(probs_cell)
     !Compute the weighted histogram
     do l=1,nb
      whist(i,j,l) = sum(probs_cell*hist_cell(:,l))
     enddo
   endif
  enddo 
 enddo

end subroutine
