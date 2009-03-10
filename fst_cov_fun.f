! def my_Gt_fun(tx,ty,scal_t,t_lim_corr,sin_frac,space_diff):
!     """
!     Converts two vectors of times, tx and ty, into a 
!     matrix whose i,j'th entry is gamma(abs(t[i]-t[j])),
!     gamma being Stein's 'valid variogram'. Parameters of
!     this variogram are:
!     - scal_t: Scales time.
!     - t_lim_corr: The limiting correlation between two points
!       as the distance between them grows. Note that this will be
!       multiplied by the overall 'amp' parameter.
!     - space_diff: The desired degree of differentiability of
!       the spatial margin.
!     - sin_frac: The fraction of the partial sill taken up by the first harmonic.
!     """
! 
!     k = t_lim_corr/space_diff
!     c = 1./space_diff-k
!     dt = np.asarray(abs(np.asmatrix(tx).T-ty))
!     Gt = 1./((np.exp(-dt/scal_t)*(1.-sin_frac) + sin_frac*np.cos(2.*np.pi*dt))*c+k)
!     return Gt, 1./(k+c)
! 

      SUBROUTINE my_gt_fun(D,x,y,nx,ny,st,tlc,sf,sd,
     *cmin,cmax,symm,origin_val)

cf2py intent(inplace) D
cf2py integer intent(optional) :: cmin=0
cf2py integer intent(optional) :: cmax=-1
cf2py logical intent(optional) :: symm=0
cf2py intent(out) origin_val
cf2py intent(hide) nx, ny
cf2py threadsafe

      INTEGER nx,ny,i,j,cmin,cmax
      DOUBLE PRECISION D(nx,ny), x(nx), y(ny)
      LOGICAL symm
      DOUBLE PRECISION st,tlc,sf,sd,k,c,dt,one,pi
      DOUBLE PRECISION origin_val
      
      PARAMETER (pi=3.141592653589793238462643d0)         

      one=1.0D0

      if (cmax.EQ.-1) then
          cmax = ny
      end if
      
      k=tlc/sd
      c=one/sd-k
      
      origin_val = one/(k+c)

      if(symm) then

        do j=cmin+1,cmax
          D(j,j) = one/(k+c)
          do i=1,j-1
            dt=dabs(x(i)-y(j))
            D(i,j) = 1/((dexp(-dt/st)*(one-sf) + 
     *        sf*dcos(2*pi*dt))*c+k)            
            if (D(i,j).LE.-one) then
                print *,'WARNING my_gt_fun has written a value <= -1!'
                return
            end if
          enddo
        enddo
      else
        do j=cmin+1,cmax
          do i=1,nx
              dt=dabs(x(i)-y(j))
              D(i,j) = 1/((dexp(-dt/st)*(one-sf) + 
     *          sf*dcos(2*pi*dt))*c+k)
             if (D(i,j).LE.-one) then
                 print *,'WARNING my_gt_fun has written a value <= -1!'
                 return
             end if
     
          enddo    
        enddo  
      endif
      RETURN
      END
