!
! SUBROUTINES FOR READING COMMAND LINE SWITCHES.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine interprets switches (i.e., seperates filenames from
! switch arguments)
!

Subroutine readswitch(options,nopts)

Implicit None

Integer, intent(in) :: nopts
Character(len=*), dimension(nopts,2), intent(inout) :: options
Character*80 buffer
Integer nswitch,status,lastswitch,newswitch,i
Integer locate

nswitch = nargs()

If (nswitch.EQ.1) then
  Call help()
End if

lastswitch=-1
Do i=1,nswitch-1
  Call getarg(i,buffer,status)
  
  newswitch=locate(buffer,options(:,1),nopts)
  If (newswitch.NE.-1) then
    If (lastswitch.NE.-1) then
      Write(6,*) "ERROR: No value for switch "//options(lastswitch,1)
      Stop
    End if
  Else
    If (lastswitch.NE.-1) then
      options(lastswitch,2)=buffer
    Else
      ! Later - Non-switches are assumed to be filenames
      Write(6,*) "ERROR: No switch specified for value "//buffer
      Stop
    End if
  Endif
  lastswitch=newswitch
 
Enddo

If (newswitch.NE.-1) then
  Write(6,*) "ERROR: No value for switch "//options(newswitch,1)
  Stop
End If

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function locates the chosen switch in the list
!

Integer Function locate(ser,list,n)

Implicit None

Integer, intent(in) :: n
Character(len=*), intent(in) :: ser
Character(len=*), dimension(n), intent(in) :: list
Integer i

locate=-1
Do i=1,n
  If (ser.EQ.list(i)) then
    locate=i
  End if
Enddo

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function returns an option value
!

Character*80 Function returnoption(flag,options,nopts)

Implicit None

Integer, intent(in) :: nopts
Character(len=*), dimension(nopts,2), intent(in) :: options
Character(len=*), intent(in) :: flag
Integer pos
Integer locate

pos=locate(flag,options(:,1),nopts)
returnoption=trim(options(pos,2))

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads the namelist
!

Subroutine procnml(num,name,maxnum)

Implicit None

Integer, intent(in) :: maxnum
Integer, intent(out) :: num
Character(len=*), dimension(1:maxnum), intent(in) :: name
Integer i,j
Logical blank

blank=.FALSE.
i=0
Do While (.NOT.blank)
  i=i+1
  If ((Ichar(name(i)(1:1)).EQ.32).OR.(Ichar(name(i)(1:1)).EQ.0)) Then
    blank=.TRUE.
  End if
  
  If (i.EQ.maxnum) Then
    Write(6,*) "ERROR: maxnum reached in namelist"
    Stop
  End If

End Do

num=i-1

If (num.EQ.0) Then
  Write(6,*) "ERROR: Zero elements in namelist"
  Stop
End If

Return
End
