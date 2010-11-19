! Author: Anand Patil
! Date: 19 Nov 2010
! License: Creative Commons BY-NC-SA


! def gtf(tx,ty,scal_t,t_lim_corr,sin_frac,space_diff):
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
! # If parameter values are illegal, just return zeros.
! # This case will be caught by the Potential.
! sd = np.add.outer(diff_degree(x)*.5, diff_degree(y)*.5)
! k=kwds['tlc']/sd
! c=1./sd-k
! sf=kwds['sf']
! tlc=kwds['tlc']
! 
! if -sd >= 1./(-sf*(1-tlc)+tlc):
!     return np.zeros((nx,ny),order='F')

      SUBROUTINE gtf(D,x,y,ddx,ddy,st,tlc,sf,nx,ny,
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
      DOUBLE PRECISION ddx(nx), ddy(ny)
      LOGICAL symm
      DOUBLE PRECISION st,tlc,sf,sd,k,c,dt,one,pi
      DOUBLE PRECISION origin_val(nx,ny)
      
      PARAMETER (pi=3.141592653589793238462643d0)         

      one=1.0D0

      if (cmax.EQ.-1) then
          cmax = ny
      end if

      if(symm) then

        do j=cmin+1,cmax
          do i=1,j
            sd = (ddx(i)+ddy(j))*0.5D0
            k = tlc/sd
            c = one/sd-k
            origin_val(i,j) = one/(k+c)
            if ((-1.0D0*sd).GE.(one/(tlc-sf*(1.0D0-tlc)))) then
                D(i,j)=-1.0D0
            else
                dt=dabs(x(i)-y(j))
                D(i,j) = 1/((dexp(-dt/st)*(one-sf) + 
     *                   sf*dcos(2*pi*dt))*c+k)            
                if (D(i,j).LE.-one) then
                    print *,'WARNING my_gt_fun wrote <= -1!'
                    return
                end if
            end if
          enddo
        enddo
      else
        do j=cmin+1,cmax
          do i=1,nx
              sd = (ddx(i)+ddy(j))*0.5D0
              k = tlc/sd
              c = one/sd-k
              origin_val(i,j) = one/(k+c)
              if ((-1.0D0*sd).GE.(one/(tlc-sf*(1.0D0-tlc)))) then
                  D(i,j)=-1.0D0
              else
                  dt=dabs(x(i)-y(j))
                  D(i,j) = 1/((dexp(-dt/st)*(one-sf) + 
     *                     sf*dcos(2*pi*dt))*c+k)
                 if (D(i,j).LE.-one) then
                     print *,'WARNING my_gt_fun wrote <= -1!'
                     return
                 end if
              end if
          enddo    
        enddo  
      endif
      RETURN
      END
