
! boundary.f90
!-------------------------------------------------------------------

module bcmod

use datamod

implicit none

integer, parameter :: isd = 1 !inside
integer, parameter :: osd = 2 !outside

contains

!-------------------------------------------------------------------

subroutine bc_periodic(dirc)

integer,intent(in) :: dirc
integer :: i,j

select case(dirc)
	case(q1dirc)
		do j=q2imin,q2imax
	    	ctr(q1imin-1,j)%w = ctr(q1imax,j)%w
	    	ctr(q1imin-1,j)%h = ctr(q1imax,j)%h
	    	ctr(q1imin-1,j)%b = ctr(q1imax,j)%b
	    	ctr(q1imin-1,j)%sh = ctr(q1imax,j)%sh
		    ctr(q1imin-1,j)%sb = ctr(q1imax,j)%sb

		    ctr(q1imax+1,j)%w = ctr(q1imin,j)%w
		    ctr(q1imax+1,j)%h = ctr(q1imin,j)%h
		    ctr(q1imax+1,j)%b = ctr(q1imin,j)%b
		    ctr(q1imax+1,j)%sh = ctr(q1imin,j)%sh
		    ctr(q1imax+1,j)%sb = ctr(q1imin,j)%sb
		end do
	case(q2dirc)
	    do i=q1imin,q1imax
		    ctr(i,q2imin-1)%w = ctr(i,q2imax)%w
		    ctr(i,q2imin-1)%h = ctr(i,q2imax)%h
		    ctr(i,q2imin-1)%b = ctr(i,q2imax)%b
		    ctr(i,q2imin-1)%sh = ctr(i,q2imax)%sh
		    ctr(i,q2imin-1)%sb = ctr(i,q2imax)%sb

		    ctr(i,q2imax+1)%w = ctr(i,q2imin)%w
		    ctr(i,q2imax+1)%h = ctr(i,q2imin)%h
		    ctr(i,q2imax+1)%b = ctr(i,q2imin)%b
		    ctr(i,q2imax+1)%sh = ctr(i,q2imin)%sh
		    ctr(i,q2imax+1)%sb = ctr(i,q2imin)%sb
		end do
end select

end subroutine bc_periodic

!-------------------------------------------------------------------

subroutine bc_extrapolation(dirc)

integer,intent(in) :: dirc
integer :: i,j

select case(dirc)
	case(q1dirc)
		do j=q2imin,q2imax
			if(geometry_method==recgeometry) then
		    	ctr(q1imax+1,j)%w = ctr(q1imax,j)%w
		    	ctr(q1imax+1,j)%h = ctr(q1imax,j)%h
		    	ctr(q1imax+1,j)%b = ctr(q1imax,j)%b
		    	ctr(q1imax+1,j)%sh = ctr(q1imax,j)%sh
			    ctr(q1imax+1,j)%sb = ctr(q1imax,j)%sb

			    ctr(q1imax+2,j)%w = ctr(q1imax,j)%w
		    	ctr(q1imax+2,j)%h = ctr(q1imax,j)%h
		    	ctr(q1imax+2,j)%b = ctr(q1imax,j)%b
		    	ctr(q1imax+2,j)%sh = ctr(q1imax,j)%sh
			    ctr(q1imax+2,j)%sb = ctr(q1imax,j)%sb
		    else if(geometry_method==polargeometry) then
		    	if(cos(ctr(q1imax+1,j)%theta)>=0) then
		    		ctr(q1imax+1,j)%w = ctr(q1imax,j)%w
			    	ctr(q1imax+1,j)%h = ctr(q1imax,j)%h
			    	ctr(q1imax+1,j)%b = ctr(q1imax,j)%b
			    	ctr(q1imax+1,j)%sh = ctr(q1imax,j)%sh
				    ctr(q1imax+1,j)%sb = ctr(q1imax,j)%sb

				    ctr(q1imax+2,j)%w = ctr(q1imax,j)%w
			    	ctr(q1imax+2,j)%h = ctr(q1imax,j)%h
			    	ctr(q1imax+2,j)%b = ctr(q1imax,j)%b
			    	ctr(q1imax+2,j)%sh = ctr(q1imax,j)%sh
				    ctr(q1imax+2,j)%sb = ctr(q1imax,j)%sb
			    end if
	    	else
		    end if
		end do
	case(q2dirc)
	    do i=q1imin,q1imax
	    	if(geometry_method==recgeometry) then
			    ctr(i,q2imax+1)%w = ctr(i,q2imax)%w
			    ctr(i,q2imax+1)%h = ctr(i,q2imax)%h
			    ctr(i,q2imax+1)%b = ctr(i,q2imax)%b
			    ctr(i,q2imax+1)%sh = ctr(i,q2imax)%sh
			    ctr(i,q2imax+1)%sb = ctr(i,q2imax)%sb

			    ctr(i,q2imax+2)%w = ctr(i,q2imax)%w
			    ctr(i,q2imax+2)%h = ctr(i,q2imax)%h
			    ctr(i,q2imax+2)%b = ctr(i,q2imax)%b
			    ctr(i,q2imax+2)%sh = ctr(i,q2imax)%sh
			    ctr(i,q2imax+2)%sb = ctr(i,q2imax)%sb
		    end if
		end do
end select

end subroutine bc_extrapolation

!-------------------------------------------------------------------

subroutine bc_insulate(dirc,side)

integer,intent(in) :: dirc,side
integer :: i,j
integer :: k,l

if(dirc==q1dirc) then
	if(side==isd) then
		do j=q2imin,q2imax
	        ctr(q1imin-1,j)%w(1) = ctr(q1imin,j)%w(1)
	        ctr(q1imin-1,j)%w(2) = -ctr(q1imin,j)%w(2)
	        ctr(q1imin-1,j)%w(3) = -ctr(q1imin,j)%w(3)
	        ctr(q1imin-1,j)%w(4) = ctr(q1imin,j)%w(4)

	        do k=1,unum
	            do l=1,vnum
	                ctr(q1imin-1,j)%h(k,l) = ctr(q1imin,j)%h(unum-k+1,vnum-l+1)
	                ctr(q1imin-1,j)%b(k,l) = ctr(q1imin,j)%b(unum-k+1,vnum-l+1)
	                ctr(q1imin-1,j)%sh = 0.d0
	                ctr(q1imin-1,j)%sb = 0.d0
	            enddo
	        enddo

	        ctr(q1imin-2,j)%w(1) = ctr(q1imin+1,j)%w(1)
	        ctr(q1imin-2,j)%w(2) = -ctr(q1imin+1,j)%w(2)
	        ctr(q1imin-2,j)%w(3) = -ctr(q1imin+1,j)%w(3)
	        ctr(q1imin-2,j)%w(4) = ctr(q1imin+1,j)%w(4)

	        do k=1,unum
	            do l=1,vnum
	                ctr(q1imin-2,j)%h(k,l) = ctr(q1imin+1,j)%h(unum-k+1,vnum-l+1)
	                ctr(q1imin-2,j)%b(k,l) = ctr(q1imin+1,j)%b(unum-k+1,vnum-l+1)
	                ctr(q1imin-2,j)%sh = 0.d0
	                ctr(q1imin-2,j)%sb = 0.d0
	            enddo
	        enddo
    	enddo
	else if(side==osd) then
		do j=q2imin,q2imax
	        ctr(q1imax+1,j)%w(1) = ctr(q1imax,j)%w(1)
	        ctr(q1imax+1,j)%w(2) = -ctr(q1imax,j)%w(2)
	        ctr(q1imax+1,j)%w(3) = -ctr(q1imax,j)%w(3)
	        ctr(q1imax+1,j)%w(4) = ctr(q1imax,j)%w(4)

	        do k=1,unum
	            do l=1,vnum
	                ctr(q1imax+1,j)%h(k,l) = ctr(q1imax,j)%h(unum-k+1,vnum-l+1)
	                ctr(q1imax+1,j)%b(k,l) = ctr(q1imax,j)%b(unum-k+1,vnum-l+1)
	                ctr(q1imax+1,j)%sh = 0.d0
	                ctr(q1imax+1,j)%sb = 0.d0
	            enddo
	        enddo

	        ctr(q1imax+2,j)%w(1) = ctr(q1imax-1,j)%w(1)
	        ctr(q1imax+2,j)%w(2) = -ctr(q1imax-1,j)%w(2)
	        ctr(q1imax+2,j)%w(3) = -ctr(q1imax-1,j)%w(3)
	        ctr(q1imax+2,j)%w(4) = ctr(q1imax-1,j)%w(4)

	        do k=1,unum
	            do l=1,vnum
	                ctr(q1imax+2,j)%h(k,l) = ctr(q1imax-1,j)%h(unum-k+1,vnum-l+1)
	                ctr(q1imax+2,j)%b(k,l) = ctr(q1imax-1,j)%b(unum-k+1,vnum-l+1)
	                ctr(q1imax+2,j)%sh = 0.d0
	                ctr(q1imax+2,j)%sb = 0.d0
	            enddo
	        enddo
    	enddo
	end if
else if(dirc==q2dirc) then
	if(side==isd) then
	    do i=q1imin,q1imax
	        ctr(i,q2imin-1)%w(1) = ctr(i,q2imin)%w(1)
	        ctr(i,q2imin-1)%w(2) = -ctr(i,q2imin)%w(2)
	        ctr(i,q2imin-1)%w(3) = -ctr(i,q2imin)%w(3)
	        ctr(i,q2imin-1)%w(4) = ctr(i,q2imin)%w(4)

	        do k=1,unum
	            do l=1,vnum
	                ctr(i,q2imin-1)%h(k,l) = ctr(i,q2imin)%h(unum-k+1,vnum-l+1)
	                ctr(i,q2imin-1)%b(k,l) = ctr(i,q2imin)%b(unum-k+1,vnum-l+1)
	                ctr(i,q2imin-1)%sh = 0.d0
	                ctr(i,q2imin-1)%sb = 0.d0
	            enddo
	        enddo

	        ctr(i,q2imin-2)%w(1) = ctr(i,q2imin+1)%w(1)
	        ctr(i,q2imin-2)%w(2) = -ctr(i,q2imin+1)%w(2)
	        ctr(i,q2imin-2)%w(3) = -ctr(i,q2imin+1)%w(3)
	        ctr(i,q2imin-2)%w(4) = ctr(i,q2imin+1)%w(4)

	        do k=1,unum
	            do l=1,vnum
	                ctr(i,q2imin-2)%h(k,l) = ctr(i,q2imin+1)%h(unum-k+1,vnum-l+1)
	                ctr(i,q2imin-2)%b(k,l) = ctr(i,q2imin+1)%b(unum-k+1,vnum-l+1)
	                ctr(i,q2imin-2)%sh = 0.d0
	                ctr(i,q2imin-2)%sb = 0.d0
	            enddo
	        enddo
    	enddo
	else if(side==osd) then
		do i=q1imin,q1imax
	        ctr(i,q2imax+1)%w(1) = ctr(i,q2imax)%w(1)
	        ctr(i,q2imax+1)%w(2) = -ctr(i,q2imax)%w(1)
	        ctr(i,q2imax+1)%w(3) = -ctr(i,q2imax)%w(1)
	        ctr(i,q2imax+1)%w(4) = ctr(i,q2imax)%w(4)

	        do k=1,unum
	            do l=1,vnum
	                ctr(i,q2imax+1)%h(k,l) = ctr(i,q2imax)%h(unum-k+1,vnum-l+1)
	                ctr(i,q2imax+1)%b(k,l) = ctr(i,q2imax)%b(unum-k+1,vnum-l+1)
	                ctr(i,q2imax+1)%sh = 0.d0
	                ctr(i,q2imax+1)%sb = 0.d0
	            enddo
	        enddo

	        ctr(i,q2imax+2)%w(1) = ctr(i,q2imax-1)%w(1)
	        ctr(i,q2imax+2)%w(2) = -ctr(i,q2imax-1)%w(2)
	        ctr(i,q2imax+2)%w(3) = -ctr(i,q2imax-1)%w(3)
	        ctr(i,q2imax+2)%w(4) = ctr(i,q2imax-1)%w(4)

	        do k=1,unum
	            do l=1,vnum
	                ctr(i,q2imax+2)%h(k,l) = ctr(i,q2imax-1)%h(unum-k+1,vnum-l+1)
	                ctr(i,q2imax+2)%b(k,l) = ctr(i,q2imax-1)%b(unum-k+1,vnum-l+1)
	                ctr(i,q2imax+2)%sh = 0.d0
	                ctr(i,q2imax+2)%sb = 0.d0
	            enddo
	        enddo
    	enddo
	end if
end if

end subroutine bc_insulate

!-------------------------------------------------------------------

end module bcmod