MODULE MOD_ISO_VARYING_STRING
  
! Written by J.L.Schonfelder 
! Incorporating suggestions by C.Tanasescu, C.Weber, J.Wagener and W.Walter,
! and corrections due to L.Moss, M.Cohen, P.Griffiths, B.T.Smith
! and many other members of the committee ISO/IEC JTC1/SC22/WG5

! Version produced (??-Jul-94)

!-----------------------------------------------------------------------------! 
! This module defines the interface and one possible implementation for a     ! 
! dynamic length character string facility in Fortran 90. The Fortran 90      !
! language is defined by the standard ISO/IEC 1539 : 1991.                    !
! The publicly accessible interface defined by this module is conformant      !
! with the auxilliary standard, ISO/IEC 1539-2 : 1994.                        !
! The detailed implementation may be considered as an informal definition of  !
! the required semantics, and may also be used as a guide to the production   !
! of a portable implementation.                                               !
! N.B. Although every care has been taken to produce valid Fortran code in    ! 
!      construction of this module no guarantee is given or implied that this ! 
!      code will work correctly without error on any specific processor, nor  !
!      is this implementation intended to be in any way optimal either in use !
!      of storage or CPU cycles.                                              !
!-----------------------------------------------------------------------------! 
  
PRIVATE 
  
!-----------------------------------------------------------------------------! 
! By default all entities declared or defined in this module are private to   !
! the module. Only those entities declared explicitly as being public are     ! 
! accessible to programs using the module. In particular, the procedures and  ! 
! operators defined herein are made accessible via their generic identifiers  ! 
! only; their specific names are private.                                     ! 
!-----------------------------------------------------------------------------! 
  
TYPE VARYING_STRING
 PRIVATE 
 CHARACTER,DIMENSION(:),POINTER :: chars 
ENDTYPE VARYING_STRING 
  
!-----------------------------------------------------------------------------! 
! The representation chosen for this definition of the module is of a string  ! 
! type consisting of a single component that is a pointer to a rank one array ! 
! of characters.                                                              ! 
! Note: this Module is defined only for characters of default kind. A similar ! 
! module could be defined for non-default characters if these are supported   ! 
! on a processor by adding a KIND parameter to the component in the type      ! 
! definition, and to all delarations of objects of CHARACTER type.            ! 
!-----------------------------------------------------------------------------! 
  
CHARACTER,PARAMETER :: blank = " " 
  
!----- GENERIC PROCEDURE INTERFACE DEFINITIONS -------------------------------! 
  
!----- LEN interface ---------------------------------------------------------! 
INTERFACE LEN 
  MODULE PROCEDURE len_s   ! length of string
ENDINTERFACE 
  
!----- Conversion procedure interfaces ---------------------------------------!
INTERFACE VAR_STR
  MODULE PROCEDURE c_to_s   ! character to string
ENDINTERFACE 
  
INTERFACE CHAR
  MODULE PROCEDURE s_to_c, &   ! string to character
                   s_to_fix_c  ! string to specified length character
ENDINTERFACE 
  
!----- ASSIGNMENT interfaces -------------------------------------------------! 
INTERFACE ASSIGNMENT(=) 
  MODULE PROCEDURE s_ass_s, &   ! string    = string
                   c_ass_s, &   ! character = string
                   s_ass_c      ! string    = character
ENDINTERFACE 
  
!----- Concatenation operator interfaces -------------------------------------!
INTERFACE OPERATOR(//) 
  MODULE PROCEDURE s_concat_s, &  ! string//string
                   s_concat_c, &  ! string//character
                   c_concat_s     ! character//string
ENDINTERFACE 
  
!----- Repeated Concatenation interface --------------------------------------! 
INTERFACE REPEAT 
  MODULE PROCEDURE repeat_s
ENDINTERFACE 
  
!------ Equality comparison operator interfaces-------------------------------!
INTERFACE OPERATOR(==) 
  MODULE PROCEDURE s_eq_s, &  ! string==string
                   s_eq_c, &  ! string==character
                   c_eq_s     ! character==string
ENDINTERFACE 
  
!----- not-equality comparison operator interfaces ---------------------------! 
INTERFACE OPERATOR(/=) 
  MODULE PROCEDURE s_ne_s, &  ! string/=string
                   s_ne_c, &  ! string/=character
                   c_ne_s     ! character/=string
ENDINTERFACE 
  
!----- less-than comparison operator interfaces ------------------------------!
INTERFACE OPERATOR(<) 
  MODULE PROCEDURE s_lt_s, &  ! string<string
                   s_lt_c, &  ! string<character
                   c_lt_s     ! character<string
ENDINTERFACE 
  
!----- less-than-or-equal comparison operator interfaces ---------------------! 
INTERFACE OPERATOR(<=) 
  MODULE PROCEDURE s_le_s, &  ! string<=string
                   s_le_c, &  ! string<=character
                   c_le_s     ! character<=string
ENDINTERFACE 
  
!----- greater-than-or-equal comparison operator interfaces ------------------!
INTERFACE OPERATOR(>=) 
  MODULE PROCEDURE s_ge_s, &  ! string>=string
                   s_ge_c, &  ! string>=character
                   c_ge_s     ! character>=string
ENDINTERFACE 
  
!----- greater-than comparison operator interfaces ---------------------------! 
INTERFACE OPERATOR(>) 
  MODULE PROCEDURE s_gt_s, &  ! string>string
                   s_gt_c, &  ! string>character
                   c_gt_s     ! character>string
ENDINTERFACE 
  
!----- LLT procedure interfaces ----------------------------------------------!
INTERFACE LLT 
  MODULE PROCEDURE s_llt_s, &  ! LLT(string,string)
                   s_llt_c, &  ! LLT(string,character)
                   c_llt_s     ! LLT(character,string)
ENDINTERFACE 
  
!----- LLE procedure interfaces ----------------------------------------------! 
INTERFACE LLE 
  MODULE PROCEDURE s_lle_s, &  ! LLE(string,string)
                   s_lle_c, &  ! LLE(string,character)
                   c_lle_s     ! LLE(character,string)
ENDINTERFACE 
  
!----- LGE procedure interfaces ----------------------------------------------!
INTERFACE LGE 
  MODULE PROCEDURE s_lge_s, &  ! LGE(string,string)
                   s_lge_c, &  ! LGE(string,character)
                   c_lge_s     ! LGE(character,string)
ENDINTERFACE 
  
!----- LGT procedure interfaces ----------------------------------------------! 
INTERFACE LGT 
  MODULE PROCEDURE s_lgt_s, &  ! LGT(string,string)
                   s_lgt_c, &  ! LGT(string,character)
                   c_lgt_s     ! LGT(character,string)
ENDINTERFACE 
  
!----- Input procedure interfaces --------------------------------------------!
INTERFACE GET
  MODULE PROCEDURE get_d_eor, &    ! default unit, EoR termination
                   get_u_eor, &    ! specified unit, EoR termination
                   get_d_tset_s, & ! default unit, string set termination
                   get_u_tset_s, & ! specified unit, string set termination
                   get_d_tset_c, & ! default unit, char set termination
                   get_u_tset_c    ! specified unit, char set termination
ENDINTERFACE 
  
!----- Output procedure interfaces -------------------------------------------!
INTERFACE PUT
  MODULE PROCEDURE put_d_s, & ! string to default unit
                   put_u_s, & ! string to specified unit
                   put_d_c, & ! char to default unit
                   put_u_c    ! char to specified unit
ENDINTERFACE 
  
INTERFACE PUT_LINE
  MODULE PROCEDURE putline_d_s, & ! string to default unit
                   putline_u_s, & ! string to specified unit
                   putline_d_c, & ! char to default unit
                   putline_u_c    ! char to specified unit
ENDINTERFACE 
  
!----- Insert procedure interfaces -------------------------------------------!
INTERFACE INSERT 
  MODULE PROCEDURE insert_ss, & ! string in string
                   insert_sc, & ! char in string
                   insert_cs, & ! string in char
                   insert_cc    ! char in char
ENDINTERFACE 

!----- Replace procedure interfaces ------------------------------------------!
INTERFACE REPLACE 
  MODULE PROCEDURE replace_ss, &   ! string by string, at specified
                   replace_sc, &   ! string by char  , starting
                   replace_cs, &   ! char by string  , point
                   replace_cc, &   ! char by char
                   replace_ss_sf,& ! string by string, between
                   replace_sc_sf,& ! string by char  , specified
                   replace_cs_sf,& ! char by string  , starting and
                   replace_cc_sf,& ! char by char    , finishing points
                   replace_sss, &  ! in string replace string by string
                   replace_ssc, &  ! in string replace string by char
                   replace_scs, &  ! in string replace char by string
                   replace_scc, &  ! in string replace char by char
                   replace_css, &  ! in char replace string by string
                   replace_csc, &  ! in char replace string by char
                   replace_ccs, &  ! in char replace char by string
                   replace_ccc     ! in char replace char by char
ENDINTERFACE 

!----- Remove procedure interface --------------------------------------------! 
INTERFACE REMOVE 
  MODULE PROCEDURE remove_s, & ! characters from string, between start
                   remove_c    ! characters from char  , and finish
ENDINTERFACE 
  
!----- Extract procedure interface -------------------------------------------!
INTERFACE EXTRACT 
  MODULE PROCEDURE extract_s, & ! from string extract string, between start
                   extract_c    ! from char   extract string, and finish
ENDINTERFACE 
  
!----- Split procedure interface ---------------------------------------------!
INTERFACE SPLIT
  MODULE PROCEDURE split_s, & ! split string at first occurance of
                   split_c    !   character in set
ENDINTERFACE

!----- Index procedure interfaces --------------------------------------------!
INTERFACE INDEX 
  MODULE PROCEDURE index_ss, index_sc, index_cs
ENDINTERFACE 
  
!----- Scan procedure interfaces ---------------------------------------------!
INTERFACE SCAN 
  MODULE PROCEDURE scan_ss, scan_sc, scan_cs
ENDINTERFACE 
  
!----- Verify procedure interfaces -------------------------------------------!
INTERFACE VERIFY 
  MODULE PROCEDURE verify_ss, verify_sc, verify_cs
ENDINTERFACE 
  
!----- Interfaces for remaining intrinsic function overloads -----------------!
INTERFACE LEN_TRIM 
  MODULE PROCEDURE len_trim_s
ENDINTERFACE 
  
INTERFACE TRIM 
  MODULE PROCEDURE trim_s
ENDINTERFACE 
  
INTERFACE IACHAR
  MODULE PROCEDURE iachar_s
ENDINTERFACE 
  
INTERFACE ICHAR 
  MODULE PROCEDURE ichar_s
ENDINTERFACE 
  
INTERFACE ADJUSTL 
  MODULE PROCEDURE adjustl_s
ENDINTERFACE 
  
INTERFACE ADJUSTR 
  MODULE PROCEDURE adjustr_s
ENDINTERFACE 
  
!----- specification of publically accessible entities -----------------------! 
PUBLIC :: VARYING_STRING,VAR_STR,CHAR,LEN,GET,PUT,PUT_LINE,INSERT,REPLACE,   &
          SPLIT,REMOVE,REPEAT,EXTRACT,INDEX,SCAN,VERIFY,LLT,LLE,LGE,LGT,     &
          ASSIGNMENT(=),OPERATOR(//),OPERATOR(==),OPERATOR(/=),OPERATOR(<),  &
          OPERATOR(<=),OPERATOR(>=),OPERATOR(>),LEN_TRIM,TRIM,IACHAR,ICHAR,  &
          ADJUSTL,ADJUSTR
  
CONTAINS 
  
!----- LEN Procedure ---------------------------------------------------------! 
 FUNCTION len_s(string) 
  type(VARYING_STRING),INTENT(IN) :: string 
  INTEGER                         :: len_s 
  ! returns the length of the string argument or zero if there is no current 
  ! string value 
  IF(.NOT.ASSOCIATED(string%chars))THEN 
    len_s = 0 
  ELSE 
    len_s = SIZE(string%chars) 
  ENDIF 
 ENDFUNCTION len_s 
  
!----- Conversion Procedures ------------------------------------------------! 
 FUNCTION c_to_s(chr) 
  type(VARYING_STRING)        :: c_to_s 
  CHARACTER(LEN=*),INTENT(IN) :: chr 
  ! returns the string consisting of the characters char 
  INTEGER                     :: lc,i 
  lc=LEN(chr) 
  ALLOCATE(c_to_s%chars(1:lc)) 
  DO i=1,lc 
    c_to_s%chars(i) = chr(i:i) 
  ENDDO 
 ENDFUNCTION c_to_s 
  
 FUNCTION s_to_c(string) 
  type(VARYING_STRING),INTENT(IN)   :: string 
  CHARACTER(LEN=SIZE(string%chars)) :: s_to_c 
  ! returns the characters of string as an automatically sized character 
  INTEGER                           :: lc,i 
  lc=SIZE(string%chars) 
  DO i=1,lc 
    s_to_c(i:i) = string%chars(i) 
  ENDDO 
 ENDFUNCTION s_to_c 
  
 FUNCTION s_to_fix_c(string,length)
  type(VARYING_STRING),INTENT(IN) :: string
  INTEGER,INTENT(IN)              :: length
  CHARACTER(LEN=length)           :: s_to_fix_c
  ! returns the character of fixed length, length, containing the characters
  ! of string either padded with blanks or truncated on the right to fit
  INTEGER                         :: lc,i
  lc=MIN(SIZE(string%chars),length)
  DO i=1,lc 
    s_to_fix_c(i:i) = string%chars(i)
  ENDDO 
  IF(lc < length)THEN  ! result longer than string padding needed
    s_to_fix_c(lc+1:length) = blank
  ENDIF
 ENDFUNCTION s_to_fix_c

!----- ASSIGNMENT Procedures -------------------------------------------------!
 SUBROUTINE s_ass_s(var,expr) 
  type(VARYING_STRING),INTENT(OUT) :: var 
  type(VARYING_STRING),INTENT(IN)  :: expr 
  !  assign a string value to a string variable overriding default assignement 
  !  reallocates string variable to size of string value and copies characters 
  ALLOCATE(var%chars(1:LEN(expr)))
  var%chars = expr%chars 
 ENDSUBROUTINE s_ass_s 
  
 SUBROUTINE c_ass_s(var,expr) 
  CHARACTER(LEN=*),INTENT(OUT)    :: var 
  type(VARYING_STRING),INTENT(IN) :: expr 
  ! assign a string value to a character variable 
  ! if the string is longer than the character truncate the string on the right 
  ! if the string is shorter the character is blank padded on the right 
  INTEGER                         :: lc,ls,i
  lc = LEN(var); ls = MIN(LEN(expr),lc) 
  DO i = 1,ls 
   var(i:i) = expr%chars(i) 
  ENDDO 
  DO i = ls+1,lc 
   var(i:i) = blank 
  ENDDO 
 ENDSUBROUTINE c_ass_s 
  
 SUBROUTINE s_ass_c(var,expr)
  type(VARYING_STRING),INTENT(OUT) :: var 
  CHARACTER(LEN=*),INTENT(IN)      :: expr 
  !  assign a character value to a string variable 
  !  disassociates the string variable from its current value, allocates new 
  !  space to hold the characters and copies them from the character value 
  !  into this space.
  INTEGER                          :: lc,i
  lc = LEN(expr) 
  ALLOCATE(var%chars(1:lc)) 
  DO i = 1,lc 
    var%chars(i) = expr(i:i) 
  ENDDO 
 ENDSUBROUTINE s_ass_c 
  
!----- Concatenation operator procedures ------------------------------------! 
 FUNCTION s_concat_s(string_a,string_b)  ! string//string 
  type(VARYING_STRING),INTENT(IN) :: string_a,string_b 
  type(VARYING_STRING)            :: s_concat_s 
  INTEGER                         :: la,lb 
  la = LEN(string_a); lb = LEN(string_b)
  ALLOCATE(s_concat_s%chars(1:la+lb)) 
  s_concat_s%chars(1:la) = string_a%chars 
  s_concat_s%chars(1+la:la+lb) = string_b%chars 
 ENDFUNCTION s_concat_s 
  
 FUNCTION s_concat_c(string_a,string_b)  ! string//character
  type(VARYING_STRING),INTENT(IN) :: string_a 
  CHARACTER(LEN=*),INTENT(IN)     :: string_b 
  type(VARYING_STRING)            :: s_concat_c 
  INTEGER                         :: la,lb,i 
  la = LEN(string_a); lb = LEN(string_b)
  ALLOCATE(s_concat_c%chars(1:la+lb)) 
  s_concat_c%chars(1:la) = string_a%chars 
  DO i = 1,lb
    s_concat_c%chars(la+i) = string_b(i:i)
  ENDDO 
 ENDFUNCTION s_concat_c 
  
 FUNCTION c_concat_s(string_a,string_b)  ! character//string 
  CHARACTER(LEN=*),INTENT(IN)     :: string_a 
  type(VARYING_STRING),INTENT(IN) :: string_b 
  type(VARYING_STRING)            :: c_concat_s 
  INTEGER                         :: la,lb,i 
  la = LEN(string_a); lb = LEN(string_b)
  ALLOCATE(c_concat_s%chars(1:la+lb)) 
  DO i = 1,la 
     c_concat_s%chars(i) = string_a(i:i) 
  ENDDO 
  c_concat_s%chars(1+la:la+lb) = string_b%chars 
 ENDFUNCTION c_concat_s 
  
!----- Reapeated concatenation procedures -----------------------------------! 
FUNCTION repeat_s(string,ncopies)                                     
 type(VARYING_STRING),INTENT(IN) :: string 
 INTEGER,INTENT(IN)              :: ncopies 
 type(VARYING_STRING)            :: repeat_s
 ! Returns a string produced by the concatenation of ncopies of the 
 ! argument string 
 INTEGER                         :: lr,ls,i 
 IF (ncopies < 0) THEN 
     WRITE(*,*) " Negative ncopies requested in REPEAT" 
     STOP 
 ENDIF 
 ls = LEN(string); lr = ls*ncopies 
 ALLOCATE(repeat_s%chars(1:lr))
 DO i = 1,ncopies 
    repeat_s%chars(1+(i-1)*ls:i*ls) = string%chars 
 ENDDO 
ENDFUNCTION repeat_s 
  
!------ Equality comparison operators ----------------------------------------! 
 FUNCTION s_eq_s(string_a,string_b)  ! string==string 
  type(VARYING_STRING),INTENT(IN) :: string_a,string_b 
  LOGICAL                         :: s_eq_s 
  INTEGER                         :: la,lb 
  la = LEN(string_a); lb = LEN(string_b)
  IF (la > lb) THEN
     s_eq_s = ALL(string_a%chars(1:lb) == string_b%chars) .AND. &
              ALL(string_a%chars(lb+1:la) == blank)
  ELSEIF (la < lb) THEN
     s_eq_s = ALL(string_a%chars == string_b%chars(1:la)) .AND. &
              ALL(blank == string_b%chars(la+1:lb))
  ELSE
     s_eq_s = ALL(string_a%chars == string_b%chars) 
  ENDIF 
 ENDFUNCTION s_eq_s 

 FUNCTION s_eq_c(string_a,string_b)  ! string==character 
  type(VARYING_STRING),INTENT(IN) :: string_a 
  CHARACTER(LEN=*),INTENT(IN)     :: string_b 
  LOGICAL                         :: s_eq_c 
  INTEGER                         :: la,lb,ls,i
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls
    IF( string_a%chars(i) /= string_b(i:i) )THEN
      s_eq_c = .FALSE.; RETURN
    ENDIF
  ENDDO
  IF( la > lb .AND. ANY( string_a%chars(lb+1:la) /= blank ) )THEN
    s_eq_c = .FALSE.; RETURN
  ELSEIF( la < lb .AND. blank /= string_b(la+1:lb) )THEN
    s_eq_c = .FALSE.; RETURN
  ENDIF
  s_eq_c = .TRUE.
 ENDFUNCTION s_eq_c 
  
 FUNCTION c_eq_s(string_a,string_b)  ! character==string 
  CHARACTER(LEN=*),INTENT(IN)     :: string_a 
  type(VARYING_STRING),INTENT(IN) :: string_b 
  LOGICAL                         :: c_eq_s 
  INTEGER                         :: la,lb,ls,i
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls
    IF( string_a(i:i) /= string_b%chars(i) )THEN
      c_eq_s = .FALSE.; RETURN
    ENDIF
  ENDDO
  IF( la > lb .AND. string_a(lb+1:la) /= blank )THEN
    c_eq_s = .FALSE.; RETURN
  ELSEIF( la < lb .AND. ANY( blank /= string_b%chars(la+1:lb) ) )THEN
    c_eq_s = .FALSE.; RETURN
  ENDIF
  c_eq_s = .TRUE.
 ENDFUNCTION c_eq_s 
  
!------ Non-equality operators -----------------------------------------------! 
 FUNCTION s_ne_s(string_a,string_b)  ! string/=string 
  type(VARYING_STRING),INTENT(IN) :: string_a,string_b 
  LOGICAL                         :: s_ne_s 
  INTEGER                         :: la,lb 
  la = LEN(string_a); lb = LEN(string_b)
  IF (la > lb) THEN
     s_ne_s = ANY(string_a%chars(1:lb) /= string_b%chars) .OR. &
              ANY(string_a%chars(lb+1:la) /= blank)
  ELSEIF (la < lb) THEN
     s_ne_s = ANY(string_a%chars /= string_b%chars(1:la)) .OR. &
              ANY(blank /= string_b%chars(la+1:lb))
  ELSE
     s_ne_s = ANY(string_a%chars /= string_b%chars)
  ENDIF 
 ENDFUNCTION s_ne_s 
  
 FUNCTION s_ne_c(string_a,string_b)  ! string/=character
  type(VARYING_STRING),INTENT(IN) :: string_a 
  CHARACTER(LEN=*),INTENT(IN)     :: string_b 
  LOGICAL                         :: s_ne_c 
  INTEGER                         :: la,lb,ls,i
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls
    IF( string_a%chars(i) /= string_b(i:i) )THEN
      s_ne_c = .TRUE.; RETURN
    ENDIF
  ENDDO
  IF( la > lb .AND. ANY( string_a%chars(lb+1:la) /= blank ) )THEN
    s_ne_c = .TRUE.; RETURN
  ELSEIF( la < lb .AND. blank /= string_b(la+1:lb) )THEN
    s_ne_c = .TRUE.; RETURN
  ENDIF
  s_ne_c = .FALSE.
 ENDFUNCTION s_ne_c 
  
 FUNCTION c_ne_s(string_a,string_b)  ! character/=string 
  CHARACTER(LEN=*),INTENT(IN)     :: string_a 
  type(VARYING_STRING),INTENT(IN) :: string_b 
  LOGICAL                         :: c_ne_s 
  INTEGER                         :: la,lb,ls,i
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls
    IF( string_a(i:i) /= string_b%chars(i) )THEN
      c_ne_s = .TRUE.; RETURN
    ENDIF
  ENDDO
  IF( la > lb .AND. string_a(lb+1:la) /= blank )THEN
    c_ne_s = .TRUE.; RETURN
  ELSEIF( la < lb .AND. ANY( blank /= string_b%chars(la+1:lb) ) )THEN
    c_ne_s = .TRUE.; RETURN
  ENDIF
  c_ne_s = .FALSE.
 ENDFUNCTION c_ne_s 
  
!------ Less-than operators --------------------------------------------------! 
 FUNCTION s_lt_s(string_a,string_b)  ! string<string
  type(VARYING_STRING),INTENT(IN) :: string_a,string_b 
  LOGICAL                         :: s_lt_s 
  INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( string_a%chars(i) < string_b%chars(i) )THEN 
      s_lt_s = .TRUE.; RETURN
    ELSEIF( string_a%chars(i) > string_b%chars(i) )THEN 
      s_lt_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( blank < string_b%chars(i) )THEN
        s_lt_s = .TRUE.; RETURN
      ELSEIF( blank > string_b%chars(i) )THEN
        s_lt_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( string_a%chars(i) < blank )THEN
        s_lt_s = .TRUE.; RETURN
      ELSEIF( string_a%chars(i) > blank )THEN
        s_lt_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_lt_s = .FALSE.
 ENDFUNCTION s_lt_s 
  
 FUNCTION s_lt_c(string_a,string_b)  ! string<character 
  type(VARYING_STRING),INTENT(IN) :: string_a 
  CHARACTER(LEN=*),INTENT(IN)     :: string_b 
  LOGICAL                         :: s_lt_c 
  INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( string_a%chars(i) < string_b(i:i) )THEN 
      s_lt_c = .TRUE.; RETURN
    ELSEIF( string_a%chars(i) > string_b(i:i) )THEN 
      s_lt_c = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    IF( blank < string_b(la+1:lb) )THEN
      s_lt_c = .TRUE.; RETURN
    ELSEIF( blank > string_b(la+1:lb) )THEN
      s_lt_c = .FALSE.; RETURN
    ENDIF
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( string_a%chars(i) < blank )THEN
        s_lt_c = .TRUE.; RETURN
      ELSEIF( string_a%chars(i) > blank )THEN
        s_lt_c = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_lt_c = .FALSE.
 ENDFUNCTION s_lt_c 
  
 FUNCTION c_lt_s(string_a,string_b)  ! character<string
  CHARACTER(LEN=*),INTENT(IN)     :: string_a 
  type(VARYING_STRING),INTENT(IN) :: string_b 
  LOGICAL                         :: c_lt_s 
  INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( string_a(i:i) < string_b%chars(i) )THEN 
      c_lt_s = .TRUE.; RETURN
    ELSEIF( string_a(i:i) > string_b%chars(i) )THEN 
      c_lt_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( blank < string_b%chars(i) )THEN
        c_lt_s = .TRUE.; RETURN
      ELSEIF( blank > string_b%chars(i) )THEN
        c_lt_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    IF( string_a(lb+1:la) < blank )THEN
      c_lt_s = .TRUE.; RETURN
    ELSEIF( string_a(lb+1:la) > blank )THEN
      c_lt_s = .FALSE.; RETURN
    ENDIF
  ENDIF
  c_lt_s = .FALSE.
 ENDFUNCTION c_lt_s 

!------ Less-than-or-equal-to operators --------------------------------------! 
 FUNCTION s_le_s(string_a,string_b)  ! string<=string 
  type(VARYING_STRING),INTENT(IN) :: string_a,string_b 
  LOGICAL                         :: s_le_s 
  INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( string_a%chars(i) < string_b%chars(i) )THEN 
      s_le_s = .TRUE.; RETURN
    ELSEIF( string_a%chars(i) > string_b%chars(i) )THEN 
      s_le_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( blank < string_b%chars(i) )THEN
        s_le_s = .TRUE.; RETURN
      ELSEIF( blank > string_b%chars(i) )THEN
        s_le_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( string_a%chars(i) < blank )THEN
        s_le_s = .TRUE.; RETURN
      ELSEIF( string_a%chars(i) > blank )THEN
        s_le_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_le_s = .TRUE.
 ENDFUNCTION s_le_s 
  
 FUNCTION s_le_c(string_a,string_b)  ! string<=character
  type(VARYING_STRING),INTENT(IN) :: string_a 
  CHARACTER(LEN=*),INTENT(IN)     :: string_b 
  LOGICAL                         :: s_le_c 
  INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( string_a%chars(i) < string_b(i:i) )THEN 
      s_le_c = .TRUE.; RETURN
    ELSEIF( string_a%chars(i) > string_b(i:i) )THEN 
      s_le_c = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    IF( blank < string_b(la+1:lb) )THEN
      s_le_c = .TRUE.; RETURN
    ELSEIF( blank > string_b(la+1:lb) )THEN
      s_le_c = .FALSE.; RETURN
    ENDIF
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( string_a%chars(i) < blank )THEN
        s_le_c = .TRUE.; RETURN
      ELSEIF( string_a%chars(i) > blank )THEN
        s_le_c = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_le_c = .TRUE.
 ENDFUNCTION s_le_c 
  
 FUNCTION c_le_s(string_a,string_b)  ! character<=string 
  CHARACTER(LEN=*),INTENT(IN)     :: string_a 
  type(VARYING_STRING),INTENT(IN) :: string_b 
  LOGICAL                         :: c_le_s 
  INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( string_a(i:i) < string_b%chars(i) )THEN 
      c_le_s = .TRUE.; RETURN
    ELSEIF( string_a(i:i) > string_b%chars(i) )THEN 
      c_le_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( blank < string_b%chars(i) )THEN
        c_le_s = .TRUE.; RETURN
      ELSEIF( blank > string_b%chars(i) )THEN
        c_le_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    IF( string_a(lb+1:la) < blank )THEN
      c_le_s = .TRUE.; RETURN
    ELSEIF( string_a(lb+1:la) > blank )THEN
      c_le_s = .FALSE.; RETURN
    ENDIF
  ENDIF
  c_le_s = .TRUE.
 ENDFUNCTION c_le_s 
  
!------ Greater-than-or-equal-to operators -----------------------------------! 
 FUNCTION s_ge_s(string_a,string_b)  ! string>=string
  type(VARYING_STRING),INTENT(IN) :: string_a,string_b 
  LOGICAL                         :: s_ge_s 
  INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( string_a%chars(i) > string_b%chars(i) )THEN
      s_ge_s = .TRUE.; RETURN
    ELSEIF( string_a%chars(i) < string_b%chars(i) )THEN
      s_ge_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( blank > string_b%chars(i) )THEN
        s_ge_s = .TRUE.; RETURN
      ELSEIF( blank < string_b%chars(i) )THEN
        s_ge_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( string_a%chars(i) > blank )THEN
        s_ge_s = .TRUE.; RETURN
      ELSEIF( string_a%chars(i) < blank )THEN
        s_ge_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_ge_s = .TRUE.
 ENDFUNCTION s_ge_s 
  
 FUNCTION s_ge_c(string_a,string_b)  ! string>=character 
  type(VARYING_STRING),INTENT(IN) :: string_a 
  CHARACTER(LEN=*),INTENT(IN)     :: string_b 
  LOGICAL                         :: s_ge_c 
  INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( string_a%chars(i) > string_b(i:i) )THEN
      s_ge_c = .TRUE.; RETURN
    ELSEIF( string_a%chars(i) < string_b(i:i) )THEN
      s_ge_c = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    IF( blank > string_b(la+1:lb) )THEN
      s_ge_c = .TRUE.; RETURN
    ELSEIF( blank < string_b(la+1:lb) )THEN
      s_ge_c = .FALSE.; RETURN
    ENDIF
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( string_a%chars(i) > blank )THEN
        s_ge_c = .TRUE.; RETURN
      ELSEIF( string_a%chars(i) < blank )THEN
        s_ge_c = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_ge_c = .TRUE.
 ENDFUNCTION s_ge_c 
  
 FUNCTION c_ge_s(string_a,string_b)  ! character>=string
  CHARACTER(LEN=*),INTENT(IN)     :: string_a 
  type(VARYING_STRING),INTENT(IN) :: string_b 
  LOGICAL                         :: c_ge_s 
  INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( string_a(i:i) > string_b%chars(i) )THEN
      c_ge_s = .TRUE.; RETURN
    ELSEIF( string_a(i:i) < string_b%chars(i) )THEN
      c_ge_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( blank > string_b%chars(i) )THEN
        c_ge_s = .TRUE.; RETURN
      ELSEIF( blank < string_b%chars(i) )THEN
        c_ge_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    IF( string_a(lb+1:la) > blank )THEN
      c_ge_s = .TRUE.; RETURN
    ELSEIF( string_a(lb+1:la) < blank )THEN
      c_ge_s = .FALSE.; RETURN
    ENDIF
  ENDIF
  c_ge_s = .TRUE.
 ENDFUNCTION c_ge_s 
  
!------ Greater-than operators -----------------------------------------------! 
 FUNCTION s_gt_s(string_a,string_b)  ! string>string 
  type(VARYING_STRING),INTENT(IN) :: string_a,string_b 
  LOGICAL                         :: s_gt_s 
  INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( string_a%chars(i) > string_b%chars(i) )THEN
      s_gt_s = .TRUE.; RETURN
    ELSEIF( string_a%chars(i) < string_b%chars(i) )THEN
      s_gt_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( blank > string_b%chars(i) )THEN
        s_gt_s = .TRUE.; RETURN
      ELSEIF( blank < string_b%chars(i) )THEN
        s_gt_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( string_a%chars(i) > blank )THEN
        s_gt_s = .TRUE.; RETURN
      ELSEIF( string_a%chars(i) < blank )THEN
        s_gt_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_gt_s = .FALSE.
 ENDFUNCTION s_gt_s 
  
 FUNCTION s_gt_c(string_a,string_b)  ! string>character
  type(VARYING_STRING),INTENT(IN) :: string_a 
  CHARACTER(LEN=*),INTENT(IN)     :: string_b 
  LOGICAL                         :: s_gt_c 
  INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( string_a%chars(i) > string_b(i:i) )THEN
      s_gt_c = .TRUE.; RETURN
    ELSEIF( string_a%chars(i) < string_b(i:i) )THEN
      s_gt_c = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    IF( blank > string_b(la+1:lb) )THEN
      s_gt_c = .TRUE.; RETURN
    ELSEIF( blank < string_b(la+1:lb) )THEN
      s_gt_c = .FALSE.; RETURN
    ENDIF
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( string_a%chars(i) > blank )THEN
        s_gt_c = .TRUE.; RETURN
      ELSEIF( string_a%chars(i) < blank )THEN
        s_gt_c = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_gt_c = .FALSE.
 ENDFUNCTION s_gt_c 
  
 FUNCTION c_gt_s(string_a,string_b)  ! character>string 
  CHARACTER(LEN=*),INTENT(IN)     :: string_a 
  type(VARYING_STRING),INTENT(IN) :: string_b 
  LOGICAL                         :: c_gt_s 
  INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( string_a(i:i) > string_b%chars(i) )THEN
      c_gt_s = .TRUE.; RETURN
    ELSEIF( string_a(i:i) < string_b%chars(i) )THEN
      c_gt_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( blank > string_b%chars(i) )THEN
        c_gt_s = .TRUE.; RETURN
      ELSEIF( blank < string_b%chars(i) )THEN
        c_gt_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    IF( string_a(lb+1:la) > blank )THEN
      c_gt_s = .TRUE.; RETURN
    ELSEIF( string_a(lb+1:la) < blank )THEN
      c_gt_s = .FALSE.; RETURN
    ENDIF
  ENDIF
  c_gt_s = .FALSE.
 ENDFUNCTION c_gt_s 
  
!----- LLT procedures -------------------------------------------------------!
FUNCTION s_llt_s(string_a,string_b)  ! string_a<string_b ISO-646 ordering 
 type(VARYING_STRING),INTENT(IN) :: string_a,string_b 
 LOGICAL                         :: s_llt_s 
 ! Returns TRUE if string_a preceeds string_b in the ISO 646 collating 
 ! sequence. Otherwise the result is FALSE. The result is FALSE if both 
 ! string_a and string_b are zero length. 
 INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( LLT(string_a%chars(i),string_b%chars(i)) )THEN
      s_llt_s = .TRUE.; RETURN
    ELSEIF( LGT(string_a%chars(i),string_b%chars(i)) )THEN
      s_llt_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( LLT(blank,string_b%chars(i)) )THEN
        s_llt_s = .TRUE.; RETURN
      ELSEIF( LGT(blank,string_b%chars(i)) )THEN
        s_llt_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( LLT(string_a%chars(i),blank) )THEN
        s_llt_s = .TRUE.; RETURN
      ELSEIF( LGT(string_a%chars(i),blank) )THEN
        s_llt_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_llt_s = .FALSE.
ENDFUNCTION s_llt_s 
  
FUNCTION s_llt_c(string_a,string_b)   ! string_a<string_b ISO-646 ordering 
 type(VARYING_STRING),INTENT(IN) :: string_a 
 CHARACTER(LEN=*),INTENT(IN)     :: string_b 
 LOGICAL                         :: s_llt_c 
 INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( LLT(string_a%chars(i),string_b(i:i)) )THEN
      s_llt_c = .TRUE.; RETURN
    ELSEIF( LGT(string_a%chars(i),string_b(i:i)) )THEN
      s_llt_c = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    IF( LLT(blank,string_b(la+1:lb)) )THEN
      s_llt_c = .TRUE.; RETURN
    ELSEIF( LGT(blank,string_b(la+1:lb)) )THEN
      s_llt_c = .FALSE.; RETURN
    ENDIF
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( LLT(string_a%chars(i),blank) )THEN
        s_llt_c = .TRUE.; RETURN
      ELSEIF( LGT(string_a%chars(i),blank) )THEN
        s_llt_c = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_llt_c = .FALSE.
ENDFUNCTION s_llt_c 
  
FUNCTION c_llt_s(string_a,string_b)  ! string_a,string_b ISO-646 ordering
 CHARACTER(LEN=*),INTENT(IN)     :: string_a
 type(VARYING_STRING),INTENT(IN) :: string_b
 LOGICAL                         :: c_llt_s
 INTEGER                         :: ls,la,lb,i
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( LLT(string_a(i:i),string_b%chars(i)) )THEN
      c_llt_s = .TRUE.; RETURN
    ELSEIF( LGT(string_a(i:i),string_b%chars(i)) )THEN
      c_llt_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( LLT(blank,string_b%chars(i)) )THEN
        c_llt_s = .TRUE.; RETURN
      ELSEIF( LGT(blank,string_b%chars(i)) )THEN
        c_llt_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    IF( LLT(string_a(lb+1:la),blank) )THEN
      c_llt_s = .TRUE.; RETURN
    ELSEIF( LGT(string_a(lb+1:la),blank) )THEN
      c_llt_s = .FALSE.; RETURN
    ENDIF
  ENDIF
  c_llt_s = .FALSE.
ENDFUNCTION c_llt_s 
  
!----- LLE procedures -------------------------------------------------------! 
FUNCTION s_lle_s(string_a,string_b)  ! string_a<=string_b ISO-646 ordering 
 type(VARYING_STRING),INTENT(IN) :: string_a,string_b 
 LOGICAL                         :: s_lle_s 
 ! Returns TRUE if strings are equal or if string_a preceeds string_b in the 
 ! ISO 646 collating sequence.  Otherwise the result is FALSE.
 INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( LLT(string_a%chars(i),string_b%chars(i)) )THEN
      s_lle_s = .TRUE.; RETURN
    ELSEIF( LGT(string_a%chars(i),string_b%chars(i)) )THEN
      s_lle_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( LLT(blank,string_b%chars(i)) )THEN
        s_lle_s = .TRUE.; RETURN
      ELSEIF( LGT(blank,string_b%chars(i)) )THEN
        s_lle_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( LLT(string_a%chars(i),blank) )THEN
        s_lle_s = .TRUE.; RETURN
      ELSEIF( LGT(string_a%chars(i),blank) )THEN
        s_lle_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_lle_s = .TRUE.
ENDFUNCTION s_lle_s 
  
FUNCTION s_lle_c(string_a,string_b)  ! strung_a<=string_b ISO-646 ordering
 type(VARYING_STRING),INTENT(IN) :: string_a 
 CHARACTER(LEN=*),INTENT(IN)     :: string_b 
 LOGICAL                         :: s_lle_c 
 INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( LLT(string_a%chars(i),string_b(i:i)) )THEN
      s_lle_c = .TRUE.; RETURN
    ELSEIF( LGT(string_a%chars(i),string_b(i:i)) )THEN
      s_lle_c = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    IF( LLT(blank,string_b(la+1:lb)) )THEN
      s_lle_c = .TRUE.; RETURN
    ELSEIF( LGT(blank,string_b(la+1:lb)) )THEN
      s_lle_c = .FALSE.; RETURN
    ENDIF
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( LLT(string_a%chars(i),blank) )THEN
        s_lle_c = .TRUE.; RETURN
      ELSEIF( LGT(string_a%chars(i),blank) )THEN
        s_lle_c = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_lle_c = .TRUE.
ENDFUNCTION s_lle_c 
  
FUNCTION c_lle_s(string_a,string_b)  ! string_a<=string_b ISO-646 ordering 
 CHARACTER(LEN=*),INTENT(IN)     :: string_a 
 type(VARYING_STRING),INTENT(IN) :: string_b 
 LOGICAL                         :: c_lle_s 
 INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( LLT(string_a(i:i),string_b%chars(i)) )THEN
      c_lle_s = .TRUE.; RETURN
    ELSEIF( LGT(string_a(i:i),string_b%chars(i)) )THEN
      c_lle_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( LLT(blank,string_b%chars(i)) )THEN
        c_lle_s = .TRUE.; RETURN
      ELSEIF( LGT(blank,string_b%chars(i)) )THEN
        c_lle_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    IF( LLT(string_a(lb+1:la),blank) )THEN
      c_lle_s = .TRUE.; RETURN
    ELSEIF( LGT(string_a(lb+1:la),blank) )THEN
      c_lle_s = .FALSE.; RETURN
    ENDIF
  ENDIF
  c_lle_s = .TRUE.
ENDFUNCTION c_lle_s 
  
!----- LGE procedures -------------------------------------------------------! 
FUNCTION s_lge_s(string_a,string_b)  ! string_a>=string_b ISO-646 ordering 
 type(VARYING_STRING),INTENT(IN) :: string_a,string_b 
 LOGICAL                         :: s_lge_s 
 ! Returns TRUE if strings are equal or if string_a follows string_b in the 
 ! ISO 646 collating sequence. Otherwise the result is FALSE.
 INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( LGT(string_a%chars(i),string_b%chars(i)) )THEN
      s_lge_s = .TRUE.; RETURN
    ELSEIF( LLT(string_a%chars(i),string_b%chars(i)) )THEN
      s_lge_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( LGT(blank,string_b%chars(i)) )THEN
        s_lge_s = .TRUE.; RETURN
      ELSEIF( LLT(blank,string_b%chars(i)) )THEN
        s_lge_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( LGT(string_a%chars(i),blank) )THEN
        s_lge_s = .TRUE.; RETURN
      ELSEIF( LLT(string_a%chars(i),blank) )THEN
        s_lge_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_lge_s = .TRUE.
ENDFUNCTION s_lge_s 
  
FUNCTION s_lge_c(string_a,string_b)  ! string_a>=string_b ISO-646 ordering
 type(VARYING_STRING),INTENT(IN) :: string_a 
 CHARACTER(LEN=*),INTENT(IN)     :: string_b 
 LOGICAL                         :: s_lge_c 
 INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( LGT(string_a%chars(i),string_b(i:i)) )THEN
      s_lge_c = .TRUE.; RETURN
    ELSEIF( LLT(string_a%chars(i),string_b(i:i)) )THEN
      s_lge_c = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    IF( LGT(blank,string_b(la+1:lb)) )THEN
      s_lge_c = .TRUE.; RETURN
    ELSEIF( LLT(blank,string_b(la+1:lb)) )THEN
      s_lge_c = .FALSE.; RETURN
    ENDIF
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( LGT(string_a%chars(i),blank) )THEN
        s_lge_c = .TRUE.; RETURN
      ELSEIF( LLT(string_a%chars(i),blank) )THEN
        s_lge_c = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_lge_c = .TRUE.
ENDFUNCTION s_lge_c 
  
FUNCTION c_lge_s(string_a,string_b)  ! string_a>=string_b ISO-646 ordering 
 CHARACTER(LEN=*),INTENT(IN)     :: string_a 
 type(VARYING_STRING),INTENT(IN) :: string_b 
 LOGICAL                         :: c_lge_s 
 INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( LGT(string_a(i:i),string_b%chars(i)) )THEN
      c_lge_s = .TRUE.; RETURN
    ELSEIF( LLT(string_a(i:i),string_b%chars(i)) )THEN
      c_lge_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( LGT(blank,string_b%chars(i)) )THEN
        c_lge_s = .TRUE.; RETURN
      ELSEIF( LLT(blank,string_b%chars(i)) )THEN
        c_lge_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    IF( LGT(string_a(lb+1:la),blank) )THEN
      c_lge_s = .TRUE.; RETURN
    ELSEIF( LLT(string_a(lb+1:la),blank) )THEN
      c_lge_s = .FALSE.; RETURN
    ENDIF
  ENDIF
  c_lge_s = .TRUE.
ENDFUNCTION c_lge_s 

!----- LGT procedures -------------------------------------------------------!
FUNCTION s_lgt_s(string_a,string_b)  ! string_a>string_b ISO-646 ordering 
 type(VARYING_STRING),INTENT(IN) :: string_a,string_b 
 LOGICAL                         :: s_lgt_s 
 ! Returns TRUE if string_a follows string_b in the ISO 646 collating sequence.
 ! Otherwise the result is FALSE. The result is FALSE if both string_a and 
 ! string_b are zero length. 
 INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( LGT(string_a%chars(i),string_b%chars(i)) )THEN
      s_lgt_s = .TRUE.; RETURN
    ELSEIF( LLT(string_a%chars(i),string_b%chars(i)) )THEN
      s_lgt_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( LGT(blank,string_b%chars(i)) )THEN
        s_lgt_s = .TRUE.; RETURN
      ELSEIF( LLT(blank,string_b%chars(i)) )THEN
        s_lgt_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( LGT(string_a%chars(i),blank) )THEN
        s_lgt_s = .TRUE.; RETURN
      ELSEIF( LLT(string_a%chars(i),blank) )THEN
        s_lgt_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_lgt_s = .FALSE.
ENDFUNCTION s_lgt_s 
  
FUNCTION s_lgt_c(string_a,string_b)  ! string_a>string_b ISO-646 ordering 
 type(VARYING_STRING),INTENT(IN) :: string_a 
 CHARACTER(LEN=*),INTENT(IN)     :: string_b 
 LOGICAL                         :: s_lgt_c 
 INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( LGT(string_a%chars(i),string_b(i:i)) )THEN
      s_lgt_c = .TRUE.; RETURN
    ELSEIF( LLT(string_a%chars(i),string_b(i:i)) )THEN
      s_lgt_c = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    IF( LGT(blank,string_b(la+1:lb)) )THEN
      s_lgt_c = .TRUE.; RETURN
    ELSEIF( LLT(blank,string_b(la+1:lb)) )THEN
      s_lgt_c = .FALSE.; RETURN
    ENDIF
  ELSEIF( la > lb )THEN
    DO i = lb+1,la
      IF( LGT(string_a%chars(i),blank) )THEN
        s_lgt_c = .TRUE.; RETURN
      ELSEIF( LLT(string_a%chars(i),blank) )THEN
        s_lgt_c = .FALSE.; RETURN
      ENDIF
    ENDDO
  ENDIF
  s_lgt_c = .FALSE.
ENDFUNCTION s_lgt_c 
  
FUNCTION c_lgt_s(string_a,string_b)  ! string_a>string_b ISO-646 ordering
 CHARACTER(LEN=*),INTENT(IN)     :: string_a 
 type(VARYING_STRING),INTENT(IN) :: string_b 
 LOGICAL                         :: c_lgt_s 
 INTEGER                         :: ls,la,lb,i 
  la = LEN(string_a); lb = LEN(string_b); ls = MIN(la,lb)
  DO i = 1,ls 
    IF( LGT(string_a(i:i),string_b%chars(i)) )THEN
      c_lgt_s = .TRUE.; RETURN
    ELSEIF( LLT(string_a(i:i),string_b%chars(i)) )THEN
      c_lgt_s = .FALSE.; RETURN
    ENDIF 
  ENDDO 
  IF( la < lb )THEN
    DO i = la+1,lb
      IF( LGT(blank,string_b%chars(i)) )THEN
        c_lgt_s = .TRUE.; RETURN
      ELSEIF( LLT(blank,string_b%chars(i)) )THEN
        c_lgt_s = .FALSE.; RETURN
      ENDIF
    ENDDO
  ELSEIF( la > lb )THEN
    IF( LGT(string_a(lb+1:la),blank) )THEN
      c_lgt_s = .TRUE.; RETURN
    ELSEIF( LLT(string_a(lb+1:la),blank) )THEN
      c_lgt_s = .FALSE.; RETURN
    ENDIF
  ENDIF
  c_lgt_s = .FALSE.
ENDFUNCTION c_lgt_s 
  
  
!----- Input string procedure -----------------------------------------------!
SUBROUTINE get_d_eor(string,maxlen,iostat)
 type(VARYING_STRING),INTENT(OUT) :: string
                                  ! the string variable to be filled with
                                  ! characters read from the
                                  ! file connected to the default unit
 INTEGER,INTENT(IN),OPTIONAL      :: maxlen
                                  ! if present indicates the maximum
                                  ! number of characters that will be
                                  ! read from the file
 INTEGER,INTENT(OUT),OPTIONAL     :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
! reads string from the default unit starting at next character in the file
! and terminating at the end of record or after maxlen characters.
 CHARACTER(LEN=80) :: buffer
 INTEGER           :: ist,nch,toread,nb
 IF(PRESENT(maxlen))THEN
   toread=maxlen
 ELSE
   toread=HUGE(1)
 ENDIF
 string = "" ! clears return string
 DO  ! repeatedly read buffer and add to string until EoR
     ! or maxlen reached
   IF(toread <= 0)EXIT
   nb=MIN(80,toread)
   READ(*,FMT='(A)',ADVANCE='NO',EOR=9999,SIZE=nch,IOSTAT=ist) buffer(1:nb)
   IF( ist /= 0 )THEN 
     IF(PRESENT(iostat)) THEN 
       iostat=ist 
       RETURN 
     ELSE 
       WRITE(*,*) " Error No.",ist, &
                  " during READ_STRING of varying string on default unit"
       STOP 
     ENDIF 
   ENDIF 
   string = string //buffer(1:nb)
   toread = toread - nb
 ENDDO
 IF(PRESENT(iostat)) iostat = 0
 RETURN
 9999 string = string //buffer(1:nch) 
 IF(PRESENT(iostat)) iostat = ist
ENDSUBROUTINE get_d_eor
  
SUBROUTINE get_u_eor(unit,string,maxlen,iostat)
 INTEGER,INTENT(IN)               :: unit
                                  ! identifies the input unit which must be
                                  ! connected for sequential formatted read
 type(VARYING_STRING),INTENT(OUT) :: string
                                  ! the string variable to be filled with
                                  ! characters read from the
                                  ! file connected to the unit
 INTEGER,INTENT(IN),OPTIONAL      :: maxlen
                                  ! if present indicates the maximum
                                  ! number of characters that will be
                                  ! read from the file
 INTEGER,INTENT(OUT),OPTIONAL     :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
! reads string from unit starting at next character in the file and 
! terminating at the end of record or after maxlen characters.
 CHARACTER(LEN=80) :: buffer
 INTEGER           :: ist,nch,toread,nb
 IF(PRESENT(maxlen))THEN
   toread=maxlen
 ELSE
   toread=HUGE(1)
 ENDIF
 string=""  ! clears return string
 DO  ! repeatedly read buffer and add to string until EoR
     ! or maxlen reached
   IF(toread <= 0)EXIT
   nb=MIN(80,toread)
   READ(unit,FMT='(A)',ADVANCE='NO',EOR=9999,SIZE=nch,IOSTAT=ist) buffer(1:nb)
   IF( ist /= 0 )THEN 
     IF(PRESENT(iostat)) THEN 
       iostat=ist 
       RETURN 
     ELSE 
       WRITE(*,*) " Error No.",ist, &
                  " during READ_STRING of varying string on UNIT ",unit
       STOP 
     ENDIF 
   ENDIF 
   string = string //buffer(1:nb)
   toread = toread - nb
 ENDDO
 IF(PRESENT(iostat)) iostat = 0
 RETURN
 9999 string = string //buffer(1:nch) 
 IF(PRESENT(iostat)) iostat = ist
ENDSUBROUTINE get_u_eor

SUBROUTINE get_d_tset_s(string,set,separator,maxlen,iostat)
 type(VARYING_STRING),INTENT(OUT) :: string
                                  ! the string variable to be filled with
                                  ! characters read from the
                                  ! file connected to the default unit
 type(VARYING_STRING),INTENT(IN)  :: set
                                  ! the set of characters which if found in
                                  ! the input terminate the read
 type(VARYING_STRING),INTENT(OUT),OPTIONAL :: separator
                                  ! the actual separator character from set
                                  ! found as the input string terminator
                                  ! returned as zero length if termination
                                  ! by maxlen or EOR
 INTEGER,INTENT(IN),OPTIONAL      :: maxlen
                                  ! if present indicates the maximum
                                  ! number of characters that will be
                                  ! read from the file
 INTEGER,INTENT(OUT),OPTIONAL     :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
! reads string from the default unit starting at next character in the file and
! terminating at the end of record, occurance of a character in set,
! or after reading maxlen characters.
 CHARACTER :: buffer  ! characters must be read one at a time to detect
                      ! first terminator character in set
 INTEGER   :: ist,toread,lenset,j
 lenset = LEN(set)
 IF(PRESENT(maxlen))THEN
   toread=maxlen
 ELSE
   toread=HUGE(1)
 ENDIF
 string = ""  ! clears return string
 IF(PRESENT(separator)) separator="" ! clear the separator
 readchar:DO  ! repeatedly read buffer and add to string
   IF(toread <= 0)EXIT readchar  ! maxlen reached
   READ(*,FMT='(A)',ADVANCE='NO',EOR=9999,IOSTAT=ist) buffer
   IF( ist /= 0 )THEN 
     IF(PRESENT(iostat)) THEN 
       iostat=ist 
       RETURN 
     ELSE 
       WRITE(*,*) " Error No.",ist, &
                  " during GET of varying string on default unit"
       STOP 
     ENDIF 
   ENDIF 
   ! check for occurance of set character in buffer
     DO j = 1,lenset
       IF(buffer == set%chars(j))THEN
         IF(PRESENT(separator)) separator=buffer
         EXIT readchar  ! separator terminator found
       ENDIF
     ENDDO
   string = string//buffer
   toread = toread - 1
  ENDDO readchar
 IF(PRESENT(iostat)) iostat = 0
 RETURN
 9999 CONTINUE  ! EOR terminator read
 IF(PRESENT(iostat)) iostat = ist
ENDSUBROUTINE get_d_tset_s

SUBROUTINE get_u_tset_s(unit,string,set,separator,maxlen,iostat)
 INTEGER,INTENT(IN)               :: unit
                                  ! identifies the input unit which must be
                                  ! connected for sequential formatted read
 type(VARYING_STRING),INTENT(OUT) :: string
                                  ! the string variable to be filled with
                                  ! characters read from the
                                  ! file connected to the unit
 type(VARYING_STRING),INTENT(IN)  :: set
                                  ! the set of characters which if found in
                                  ! the input terminate the read
 type(VARYING_STRING),INTENT(OUT),OPTIONAL :: separator
                                  ! the actual separator character from set
                                  ! found as the input string terminator
                                  ! returned as zero length if termination
                                  ! by maxlen or EOR
 INTEGER,INTENT(IN),OPTIONAL      :: maxlen
                                  ! if present indicates the maximum
                                  ! number of characters that will be
                                  ! read from the file
 INTEGER,INTENT(OUT),OPTIONAL     :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
! reads string from unit starting at next character in the file and 
! terminating at the end of record, occurance of a character in set,
! or after reading maxlen characters.
 CHARACTER :: buffer  ! characters must be read one at a time to detect
                      ! first terminator character in set
 INTEGER   :: ist,toread,lenset,j
 lenset = LEN(set)
 IF(PRESENT(maxlen))THEN
   toread=maxlen
 ELSE
   toread=HUGE(1)
 ENDIF
 string = ""  ! clears return string
 IF(PRESENT(separator)) separator="" ! clear the separator
 readchar:DO  ! repeatedly read buffer and add to string
   IF(toread <= 0)EXIT readchar ! maxlen reached
   READ(unit,FMT='(A)',ADVANCE='NO',EOR=9999,IOSTAT=ist) buffer
   IF( ist /= 0 )THEN 
     IF(PRESENT(iostat)) THEN 
       iostat=ist 
       RETURN 
     ELSE 
       WRITE(*,*) " Error No.",ist, &
                  " during GET of varying string on unit ",unit
       STOP 
     ENDIF 
   ENDIF 
   ! check for occurance of set character in buffer
     DO j = 1,lenset
       IF(buffer == set%chars(j))THEN
         IF(PRESENT(separator)) separator=buffer
         EXIT readchar ! separator terminator found
       ENDIF
     ENDDO
   string = string//buffer
   toread = toread - 1
 ENDDO readchar
 IF(PRESENT(iostat)) iostat = 0
 RETURN
 9999 CONTINUE  ! EOR terminator found
 IF(PRESENT(iostat)) iostat = ist
ENDSUBROUTINE get_u_tset_s

SUBROUTINE get_d_tset_c(string,set,separator,maxlen,iostat)
 type(VARYING_STRING),INTENT(OUT) :: string
                                  ! the string variable to be filled with
                                  ! characters read from the
                                  ! file connected to the default unit
 CHARACTER(LEN=*),INTENT(IN)      :: set
                                  ! the set of characters which if found in
                                  ! the input terminate the read
 type(VARYING_STRING),INTENT(OUT),OPTIONAL :: separator
                                  ! the actual separator character from set
                                  ! found as the input string terminator
                                  ! returned as zero length if termination
                                  ! by maxlen or EOR
 INTEGER,INTENT(IN),OPTIONAL      :: maxlen
                                  ! if present indicates the maximum
                                  ! number of characters that will be
                                  ! read from the file
 INTEGER,INTENT(OUT),OPTIONAL     :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
! reads string from the default unit starting at next character in the file and
! terminating at the end of record, occurance of a character in set,
! or after reading maxlen characters.
 CHARACTER :: buffer  ! characters must be read one at a time to detect
                      ! first terminator character in set
 INTEGER   :: ist,toread,lenset,j
 lenset = LEN(set)
 IF(PRESENT(maxlen))THEN
   toread=maxlen
 ELSE
   toread=HUGE(1)
 ENDIF
 string = ""  ! clears return string
 IF(PRESENT(separator)) separator=""  ! clear separator
 readchar:DO  ! repeatedly read buffer and add to string
   IF(toread <= 0)EXIT readchar ! maxlen reached
   READ(*,FMT='(A)',ADVANCE='NO',EOR=9999,IOSTAT=ist) buffer
   IF( ist /= 0 )THEN 
     IF(PRESENT(iostat)) THEN 
       iostat=ist 
       RETURN 
     ELSE 
       WRITE(*,*) " Error No.",ist, &
                  " during GET of varying string on default unit"
       STOP 
     ENDIF 
   ENDIF 
   ! check for occurance of set character in buffer
     DO j = 1,lenset
       IF(buffer == set(j:j))THEN
         IF(PRESENT(separator)) separator=buffer
         EXIT readchar
       ENDIF
     ENDDO
   string = string//buffer
   toread = toread - 1
 ENDDO readchar
 IF(PRESENT(iostat)) iostat = 0
 RETURN
 9999 CONTINUE ! EOR terminator read
 IF(PRESENT(iostat)) iostat = ist
ENDSUBROUTINE get_d_tset_c

SUBROUTINE get_u_tset_c(unit,string,set,separator,maxlen,iostat)
 INTEGER,INTENT(IN)               :: unit
                                  ! identifies the input unit which must be
                                  ! connected for sequential formatted read
 type(VARYING_STRING),INTENT(OUT) :: string
                                  ! the string variable to be filled with
                                  ! characters read from the
                                  ! file connected to the unit
 CHARACTER(LEN=*),INTENT(IN)      :: set
                                  ! the set of characters which if found in
                                  ! the input terminate the read
 type(VARYING_STRING),INTENT(OUT),OPTIONAL :: separator
                                  ! the actual separator character from set
                                  ! found as the input string terminator
                                  ! returned as zero length if termination
                                  ! by maxlen or EOR
 INTEGER,INTENT(IN),OPTIONAL      :: maxlen
                                  ! if present indicates the maximum
                                  ! number of characters that will be
                                  ! read from the file
 INTEGER,INTENT(OUT),OPTIONAL     :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
! reads string from unit starting at next character in the file and
! terminating at the end of record, occurance of a character in set,
! or after reading maxlen characters.
 CHARACTER :: buffer  ! characters must be read one at a time to detect
                      ! first terminator character in set
 INTEGER   :: ist,toread,lenset,j
 lenset = LEN(set)
 IF(PRESENT(maxlen))THEN
   toread=maxlen
 ELSE
   toread=HUGE(1)
 ENDIF
 string = ""  ! clears return string
 IF(PRESENT(separator)) separator="" ! clear separator
 readchar:DO  ! repeatedly read buffer and add to string
   IF(toread <= 0)EXIT readchar ! maxlen reached
   READ(unit,FMT='(A)',ADVANCE='NO',EOR=9999,IOSTAT=ist) buffer
   IF( ist /= 0 )THEN 
     IF(PRESENT(iostat)) THEN 
       iostat=ist 
       RETURN 
     ELSE 
       WRITE(*,*) " Error No.",ist, &
                  " during GET of varying string on unit ",unit
       STOP 
     ENDIF 
   ENDIF 
   ! check for occurance of set character in buffer
     DO j = 1,lenset
       IF(buffer == set(j:j))THEN
         IF(PRESENT(separator)) separator=buffer
         EXIT readchar ! separator terminator found
       ENDIF
     ENDDO
   string = string//buffer
   toread = toread - 1
 ENDDO readchar
 IF(PRESENT(iostat)) iostat = 0
 RETURN
 9999 CONTINUE ! EOR terminator read
 IF(PRESENT(iostat)) iostat = ist
ENDSUBROUTINE get_u_tset_c

!----- Output string procedures ----------------------------------------------!
SUBROUTINE put_d_s(string,iostat)
 type(VARYING_STRING),INTENT(IN) :: string
                                  ! the string variable to be appended to
                                  ! the current record or to the start of
                                  ! the next record if there is no
                                  ! current record
                                  ! uses the default unit
 INTEGER,INTENT(OUT),OPTIONAL    :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
 INTEGER           :: ist
 WRITE(*,FMT='(A)',ADVANCE='NO',IOSTAT=ist) CHAR(string)
 IF( ist /= 0 )THEN
   IF(PRESENT(iostat))THEN
     iostat = ist
     RETURN
   ELSE
     WRITE(*,*) " Error No.",ist, &
                "  during PUT of varying string on default unit"
     STOP
   ENDIF
 ENDIF
 IF(PRESENT(iostat)) iostat=0
ENDSUBROUTINE put_d_s

SUBROUTINE put_u_s(unit,string,iostat)
 INTEGER,INTENT(IN)              :: unit
                                  ! identifies the output unit which must
                                  ! be connected for sequential formatted
                                  ! write
 type(VARYING_STRING),INTENT(IN) :: string
                                  ! the string variable to be appended to
                                  ! the current record or to the start of
                                  ! the next record if there is no
                                  ! current record
 INTEGER,INTENT(OUT),OPTIONAL    :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
  INTEGER :: ist
  WRITE(unit,FMT='(A)',ADVANCE='NO',IOSTAT=ist) CHAR(string)
  IF( ist /= 0 )THEN
   IF(PRESENT(iostat))THEN
    iostat = ist
    RETURN
   ELSE
    WRITE(*,*) " Error No.",ist, &
               "  during PUT of varying string on UNIT ",unit
    STOP
   ENDIF
  ENDIF
 IF(PRESENT(iostat)) iostat=0
ENDSUBROUTINE put_u_s
  
SUBROUTINE put_d_c(string,iostat)
 CHARACTER(LEN=*),INTENT(IN)     :: string
                                  ! the character variable to be appended to
                                  ! the current record or to the start of
                                  ! the next record if there is no
                                  ! current record
                                  ! uses the default unit
 INTEGER,INTENT(OUT),OPTIONAL    :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
 INTEGER :: ist
 WRITE(*,FMT='(A)',ADVANCE='NO',IOSTAT=ist) string
 IF( ist /= 0 )THEN
  IF(PRESENT(iostat))THEN
   iostat = ist
   RETURN
  ELSE
   WRITE(*,*) " Error No.",ist, &
              " during PUT of character on default unit"
   STOP
  ENDIF
 ENDIF
 IF(PRESENT(iostat)) iostat=0
ENDSUBROUTINE put_d_c

SUBROUTINE put_u_c(unit,string,iostat)
 INTEGER,INTENT(IN)              :: unit
                                  ! identifies the output unit which must
                                  ! be connected for sequential formatted
                                  ! write
 CHARACTER(LEN=*),INTENT(IN)     :: string
                                  ! the character variable to be appended to
                                  ! the current record or to the start of
                                  ! the next record if there is no
                                  ! current record
 INTEGER,INTENT(OUT),OPTIONAL    :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
 INTEGER :: ist
 WRITE(unit,FMT='(A)',ADVANCE='NO',IOSTAT=ist) string
 IF( ist /= 0 )THEN
  IF(PRESENT(iostat))THEN
   iostat = ist
   RETURN
  ELSE
   WRITE(*,*) " Error No.",ist," during PUT of character on UNIT ",unit
   STOP
  ENDIF
 ENDIF
 IF(PRESENT(iostat)) iostat=0
ENDSUBROUTINE put_u_c

SUBROUTINE putline_d_s(string,iostat)
 type(VARYING_STRING),INTENT(IN) :: string
                                  ! the string variable to be appended to
                                  ! the current record or to the start of
                                  ! the next record if there is no
                                  ! current record
                                  ! uses the default unit
 INTEGER,INTENT(OUT),OPTIONAL    :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
! appends the string to the current record and then ends the record
! leaves the file positioned after the record just completed which then
! becomes the previous and last record in the file.
 INTEGER           :: ist
  WRITE(*,FMT='(A,/)',ADVANCE='NO',IOSTAT=ist) CHAR(string)
  IF( ist /= 0 )THEN
   IF(PRESENT(iostat))THEN
    iostat = ist; RETURN
   ELSE
    WRITE(*,*) " Error No.",ist, &
               " during PUT_LINE of varying string on default unit"
    STOP
   ENDIF
  ENDIF
 IF(PRESENT(iostat)) iostat=0
ENDSUBROUTINE putline_d_s
  
SUBROUTINE putline_u_s(unit,string,iostat)
 INTEGER,INTENT(IN)              :: unit
                                  ! identifies the output unit which must
                                  ! be connected for sequential formatted
                                  ! write
 type(VARYING_STRING),INTENT(IN) :: string
                                  ! the string variable to be appended to
                                  ! the current record or to the start of
                                  ! the next record if there is no
                                  ! current record
 INTEGER,INTENT(OUT),OPTIONAL    :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
! appends the string to the current record and then ends the record
! leaves the file positioned after the record just completed which then
! becomes the previous and last record in the file.
 INTEGER  :: ist
  WRITE(unit,FMT='(A,/)',ADVANCE='NO',IOSTAT=ist) CHAR(string)
  IF( ist /= 0 )THEN
   IF(PRESENT(iostat))THEN
    iostat = ist; RETURN
   ELSE
    WRITE(*,*) " Error No.",ist, &
               " during PUT_LINE of varying string on UNIT",unit
    STOP
   ENDIF
  ENDIF
 IF(PRESENT(iostat)) iostat=0
ENDSUBROUTINE putline_u_s

SUBROUTINE putline_d_c(string,iostat)
 CHARACTER(LEN=*),INTENT(IN)     :: string
                                  ! the character variable to be appended to
                                  ! the current record or to the start of
                                  ! the next record if there is no
                                  ! current record
                                  ! uses the default unit
 INTEGER,INTENT(OUT),OPTIONAL    :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
! appends the string to the current record and then ends the record
! leaves the file positioned after the record just completed which then
! becomes the previous and last record in the file.
 INTEGER :: ist
 WRITE(*,FMT='(A,/)',ADVANCE='NO',IOSTAT=ist) string
 IF(PRESENT(iostat))THEN
  iostat = ist
  RETURN
 ELSEIF( ist /= 0 )THEN
  WRITE(*,*) " Error No.",ist, &
              " during PUT_LINE of character on default unit"
  STOP
 ENDIF
ENDSUBROUTINE putline_d_c
  
SUBROUTINE putline_u_c(unit,string,iostat)
 INTEGER,INTENT(IN)              :: unit
                                  ! identifies the output unit which must
                                  ! be connected for sequential formatted
                                  ! write
 CHARACTER(LEN=*),INTENT(IN)     :: string
                                  ! the character variable to be appended to
                                  ! the current record or to the start of
                                  ! the next record if there is no
                                  ! current record
 INTEGER,INTENT(OUT),OPTIONAL    :: iostat
                                  ! if present used to return the status
                                  ! of the data transfer
                                  ! if absent errors cause termination
! appends the string to the current record and then ends the record
! leaves the file positioned after the record just completed which then
! becomes the previous and last record in the file.
 INTEGER :: ist
 WRITE(unit,FMT='(A,/)',ADVANCE='NO',IOSTAT=ist) string
 IF(PRESENT(iostat))THEN
  iostat = ist
  RETURN
 ELSEIF( ist /= 0 )THEN
  WRITE(*,*) " Error No.",ist, &
              " during WRITE_LINE of character on UNIT",unit
  STOP
 ENDIF
ENDSUBROUTINE putline_u_c
  
!----- Insert procedures ----------------------------------------------------!
 FUNCTION insert_ss(string,start,substring)
  type(VARYING_STRING)            :: insert_ss
  type(VARYING_STRING),INTENT(IN) :: string
  INTEGER,INTENT(IN)              :: start
  type(VARYING_STRING),INTENT(IN) :: substring
  ! calculates result string by inserting the substring into string
  ! beginning at position start pushing the remainder of the string
  ! to the right and enlarging it accordingly,
  ! if start is greater than LEN(string) the substring is simply appended
  ! to string by concatenation. if start is less than 1
  ! substring is inserted before string, ie. start is treated as if it were 1 
  CHARACTER,POINTER,DIMENSION(:) :: work 
  INTEGER                        :: ip,is,lsub,ls 
  lsub = LEN(substring); ls = LEN(string)
  is = MAX(start,1) 
  ip = MIN(ls+1,is) 
  ALLOCATE(work(1:lsub+ls))
  work(1:ip-1) = string%chars(1:ip-1) 
  work(ip:ip+lsub-1) =substring%chars
  work(ip+lsub:lsub+ls) = string%chars(ip:ls)
  insert_ss%chars => work
 ENDFUNCTION insert_ss
  
 FUNCTION insert_sc(string,start,substring)
  type(VARYING_STRING)            :: insert_sc
  type(VARYING_STRING),INTENT(IN) :: string
  INTEGER,INTENT(IN)              :: start
  CHARACTER(LEN=*),INTENT(IN) :: substring
  ! calculates result string by inserting the substring into string
  ! beginning at position start pushing the remainder of the string
  ! to the right and enlarging it accordingly,
  ! if start is greater than LEN(string) the substring is simply appended
  ! to string by concatenation. if start is less than 1
  ! substring is inserted before string, ie. start is treated as if it were 1 
  CHARACTER,POINTER,DIMENSION(:) :: work 
  INTEGER                        :: ip,is,lsub,ls,i 
  lsub = LEN(substring); ls = LEN(string)
  is = MAX(start,1) 
  ip = MIN(ls+1,is) 
  ALLOCATE(work(1:lsub+ls))
  work(1:ip-1) = string%chars(1:ip-1) 
  DO i = 1,lsub 
   work(ip-1+i) = substring(i:i) 
  ENDDO 
  work(ip+lsub:lsub+ls) = string%chars(ip:ls)
  insert_sc%chars => work
 ENDFUNCTION insert_sc

 FUNCTION insert_cs(string,start,substring)
  type(VARYING_STRING)            :: insert_cs
  CHARACTER(LEN=*),INTENT(IN)     :: string
  INTEGER,INTENT(IN)              :: start
  type(VARYING_STRING),INTENT(IN) :: substring
  ! calculates result string by inserting the substring into string
  ! beginning at position start pushing the remainder of the string
  ! to the right and enlarging it accordingly,
  ! if start is greater than LEN(string) the substring is simply appended
  ! to string by concatenation. if start is less than 1
  ! substring is inserted before string, ie. start is treated as if it were 1 
  CHARACTER,POINTER,DIMENSION(:) :: work 
  INTEGER                        :: ip,is,lsub,ls,i 
  lsub = LEN(substring); ls = LEN(string)
  is = MAX(start,1) 
  ip = MIN(ls+1,is) 
  ALLOCATE(work(1:lsub+ls))
  DO i=1,ip-1
    work(i) = string(i:i)
  ENDDO
  work(ip:ip+lsub-1) =substring%chars
  DO i=ip,ls
    work(i+lsub) = string(i:i)
  ENDDO
  insert_cs%chars => work
 ENDFUNCTION insert_cs
  
 FUNCTION insert_cc(string,start,substring)
  type(VARYING_STRING)        :: insert_cc
  CHARACTER(LEN=*),INTENT(IN) :: string
  INTEGER,INTENT(IN)          :: start
  CHARACTER(LEN=*),INTENT(IN) :: substring
  ! calculates result string by inserting the substring into string
  ! beginning at position start pushing the remainder of the string
  ! to the right and enlarging it accordingly,
  ! if start is greater than LEN(string) the substring is simply appended
  ! to string by concatenation. if start is less than 1
  ! substring is inserted before string, ie. start is treated as if it were 1 
  CHARACTER,POINTER,DIMENSION(:) :: work 
  INTEGER                        :: ip,is,lsub,ls,i 
  lsub = LEN(substring); ls = LEN(string)
  is = MAX(start,1) 
  ip = MIN(ls+1,is) 
  ALLOCATE(work(1:lsub+ls))
  DO i=1,ip-1
    work(i) = string(i:i)
  ENDDO
  DO i = 1,lsub 
   work(ip-1+i) = substring(i:i) 
  ENDDO 
  DO i=ip,ls
    work(i+lsub) = string(i:i)
  ENDDO
  insert_cc%chars => work
 ENDFUNCTION insert_cc

!----- Replace procedures ---------------------------------------------------!
FUNCTION replace_ss(string,start,substring)
 type(VARYING_STRING)            :: replace_ss
 type(VARYING_STRING),INTENT(IN) :: string
 INTEGER,INTENT(IN)              :: start
 type(VARYING_STRING),INTENT(IN) :: substring
 !  calculates the result string by the following actions:
 !  inserts the substring into string beginning at position 
 !  start replacing the following LEN(substring) characters of the string 
 !  and enlarging string if necessary. if start is greater than LEN(string) 
 !  substring is simply appended to string by concatenation. If start is less 
 !  than 1, substring replaces characters in string starting at 1
 CHARACTER,POINTER,DIMENSION(:) :: work
 INTEGER                        :: ip,is,nw,lsub,ls
 lsub = LEN(substring); ls = LEN(string)
 is = MAX(start,1)
 ip = MIN(ls+1,is)
 nw = MAX(ls,ip+lsub-1)
 ALLOCATE(work(1:nw))
 work(1:ip-1) = string%chars(1:ip-1)
 work(ip:ip+lsub-1) = substring%chars
 work(ip+lsub:nw) = string%chars(ip+lsub:ls)
 replace_ss%chars => work
ENDFUNCTION replace_ss
  
FUNCTION replace_ss_sf(string,start,finish,substring)
 type(VARYING_STRING)            :: replace_ss_sf
 type(VARYING_STRING),INTENT(IN) :: string
 INTEGER,INTENT(IN)              :: start,finish
 type(VARYING_STRING),INTENT(IN) :: substring
 !  calculates the result string by the following actions:
 !  inserts the substring into string beginning at position 
 !  start replacing the following finish-start+1 characters of the string
 !  and enlarging or shrinking the string if necessary.
 !  If start is greater than LEN(string) substring is simply appended to string
 !  by concatenation. If start is less than 1, start = 1 is used
 !  If finish is greater than LEN(string), finish = LEN(string) is used
 !  If finish is less than start, substring is inserted before start
 CHARACTER,POINTER,DIMENSION(:) :: work
 INTEGER                        :: ip,is,if,nw,lsub,ls
 lsub = LEN(substring); ls = LEN(string)
 is = MAX(start,1)
 ip = MIN(ls+1,is)
 if = MAX(ip-1,MIN(finish,ls))
 nw = lsub + ls - if+ip-1
 ALLOCATE(work(1:nw))
 work(1:ip-1) = string%chars(1:ip-1)
 work(ip:ip+lsub-1) = substring%chars
 work(ip+lsub:nw) = string%chars(if+1:ls)
 replace_ss_sf%chars => work
ENDFUNCTION replace_ss_sf

FUNCTION replace_sc(string,start,substring)
 type(VARYING_STRING)            :: replace_sc
 type(VARYING_STRING),INTENT(IN) :: string
 INTEGER,INTENT(IN)              :: start
 CHARACTER(LEN=*),INTENT(IN)     :: substring
 !  calculates the result string by the following actions:
 !  inserts the characters from substring into string beginning at position 
 !  start replacing the following LEN(substring) characters of the string 
 !  and enlarging string if necessary. If start is greater than LEN(string) 
 !  substring is simply appended to string by concatenation. If start is less 
 !  than 1, substring replaces characters in string starting at 1
 CHARACTER,POINTER,DIMENSION(:) :: work
 INTEGER                        :: ip,is,nw,lsub,ls,i
 lsub = LEN(substring); ls = LEN(string)
 is = MAX(start,1)
 ip = MIN(ls+1,is)
 nw = MAX(ls,ip+lsub-1)
 ALLOCATE(work(1:nw))
 work(1:ip-1) = string%chars(1:ip-1)
 DO i = 1,lsub
   work(ip-1+i) = substring(i:i)
 ENDDO
 work(ip+lsub:nw) = string%chars(ip+lsub:ls)
 replace_sc%chars => work
ENDFUNCTION replace_sc
  
FUNCTION replace_sc_sf(string,start,finish,substring)
 type(VARYING_STRING)            :: replace_sc_sf
 type(VARYING_STRING),INTENT(IN) :: string
 INTEGER,INTENT(IN)              :: start,finish
 CHARACTER(LEN=*),INTENT(IN)     :: substring
 !  calculates the result string by the following actions:
 !  inserts the substring into string beginning at position 
 !  start replacing the following finish-start+1 characters of the string
 !  and enlarging or shrinking the string if necessary.
 !  If start is greater than LEN(string) substring is simply appended to string
 !  by concatenation. If start is less than 1, start = 1 is used
 !  If finish is greater than LEN(string), finish = LEN(string) is used
 !  If finish is less than start, substring is inserted before start
 CHARACTER,POINTER,DIMENSION(:) :: work
 INTEGER                        :: ip,is,if,nw,lsub,ls,i
 lsub = LEN(substring); ls = LEN(string)
 is = MAX(start,1)
 ip = MIN(ls+1,is)
 if = MAX(ip-1,MIN(finish,ls))
 nw = lsub + ls - if+ip-1
 ALLOCATE(work(1:nw))
 work(1:ip-1) = string%chars(1:ip-1)
 DO i = 1,lsub
   work(ip-1+i) = substring(i:i)
 ENDDO
 work(ip+lsub:nw) = string%chars(if+1:ls)
 replace_sc_sf%chars => work
ENDFUNCTION replace_sc_sf

FUNCTION replace_cs(string,start,substring)
 type(VARYING_STRING)            :: replace_cs
 CHARACTER(LEN=*),INTENT(IN)     :: string
 INTEGER,INTENT(IN)              :: start
 type(VARYING_STRING),INTENT(IN) :: substring
 !  calculates the result string by the following actions:
 !  inserts the substring into string beginning at position 
 !  start replacing the following LEN(substring) characters of the string 
 !  and enlarging string if necessary. if start is greater than LEN(string) 
 !  substring is simply appended to string by concatenation. If start is less 
 !  than 1, substring replaces characters in string starting at 1
 CHARACTER,POINTER,DIMENSION(:) :: work
 INTEGER                        :: ip,is,nw,lsub,ls,i
 lsub = LEN(substring); ls = LEN(string)
 is = MAX(start,1)
 ip = MIN(ls+1,is)
 nw = MAX(ls,ip+lsub-1)
 ALLOCATE(work(1:nw))
 DO i=1,ip-1
   work(i) = string(i:i)
 ENDDO
 work(ip:ip+lsub-1) = substring%chars
 DO i=ip+lsub,nw
   work(i) = string(i:i)
 ENDDO
 replace_cs%chars => work
ENDFUNCTION replace_cs
  
FUNCTION replace_cs_sf(string,start,finish,substring)
 type(VARYING_STRING)            :: replace_cs_sf
 CHARACTER(LEN=*),INTENT(IN)     :: string
 INTEGER,INTENT(IN)              :: start,finish
 type(VARYING_STRING),INTENT(IN) :: substring
 !  calculates the result string by the following actions:
 !  inserts the substring into string beginning at position 
 !  start replacing the following finish-start+1 characters of the string
 !  and enlarging or shrinking the string if necessary.
 !  If start is greater than LEN(string) substring is simply appended to string
 !  by concatenation. If start is less than 1, start = 1 is used
 !  If finish is greater than LEN(string), finish = LEN(string) is used
 !  If finish is less than start, substring is inserted before start
 CHARACTER,POINTER,DIMENSION(:) :: work
 INTEGER                        :: ip,is,if,nw,lsub,ls,i
 lsub = LEN(substring); ls = LEN(string)
 is = MAX(start,1)
 ip = MIN(ls+1,is)
 if = MAX(ip-1,MIN(finish,ls))
 nw = lsub + ls - if+ip-1
 ALLOCATE(work(1:nw))
 DO i=1,ip-1
   work(i) = string(i:i)
 ENDDO
 work(ip:ip+lsub-1) = substring%chars
 DO i=1,nw-ip-lsub+1
   work(i+ip+lsub-1) = string(if+i:if+i)
 ENDDO
 replace_cs_sf%chars => work
ENDFUNCTION replace_cs_sf

FUNCTION replace_cc(string,start,substring)
 type(VARYING_STRING)            :: replace_cc
 CHARACTER(LEN=*),INTENT(IN)     :: string
 INTEGER,INTENT(IN)              :: start
 CHARACTER(LEN=*),INTENT(IN)     :: substring
 !  calculates the result string by the following actions:
 !  inserts the characters from substring into string beginning at position 
 !  start replacing the following LEN(substring) characters of the string 
 !  and enlarging string if necessary. If start is greater than LEN(string) 
 !  substring is simply appended to string by concatenation. If start is less 
 !  than 1, substring replaces characters in string starting at 1
 CHARACTER,POINTER,DIMENSION(:) :: work
 INTEGER                        :: ip,is,nw,lsub,ls,i
 lsub = LEN(substring); ls = LEN(string)
 is = MAX(start,1)
 ip = MIN(ls+1,is)
 nw = MAX(ls,ip+lsub-1)
 ALLOCATE(work(1:nw))
 DO i=1,ip-1
   work(i) = string(i:i)
 ENDDO
 DO i=1,lsub
   work(ip-1+i) = substring(i:i)
 ENDDO
 DO i=ip+lsub,nw
   work(i) = string(i:i)
 ENDDO
 replace_cc%chars => work
ENDFUNCTION replace_cc
  
FUNCTION replace_cc_sf(string,start,finish,substring)
 type(VARYING_STRING)            :: replace_cc_sf
 CHARACTER(LEN=*),INTENT(IN)     :: string
 INTEGER,INTENT(IN)              :: start,finish
 CHARACTER(LEN=*),INTENT(IN)     :: substring
 !  calculates the result string by the following actions:
 !  inserts the substring into string beginning at position 
 !  start replacing the following finish-start+1 characters of the string
 !  and enlarging or shrinking the string if necessary.
 !  If start is greater than LEN(string) substring is simply appended to string
 !  by concatenation. If start is less than 1, start = 1 is used
 !  If finish is greater than LEN(string), finish = LEN(string) is used
 !  If finish is less than start, substring is inserted before start
 CHARACTER,POINTER,DIMENSION(:) :: work
 INTEGER                        :: ip,is,if,nw,lsub,ls,i
 lsub = LEN(substring); ls = LEN(string)
 is = MAX(start,1)
 ip = MIN(ls+1,is)
 if = MAX(ip-1,MIN(finish,ls))
 nw = lsub + ls - if+ip-1
 ALLOCATE(work(1:nw))
 DO i=1,ip-1
   work(i) = string(i:i)
 ENDDO
 DO i=1,lsub
   work(i+ip-1) = substring(i:i)
 ENDDO
 DO i=1,nw-ip-lsub+1
   work(i+ip+lsub-1) = string(if+i:if+i)
 ENDDO
 replace_cc_sf%chars => work
ENDFUNCTION replace_cc_sf

FUNCTION replace_sss(string,target,substring,every,back)
 type(VARYING_STRING)            :: replace_sss
 type(VARYING_STRING),INTENT(IN) :: string,target,substring
 LOGICAL,INTENT(IN),OPTIONAL     :: every,back
 !  calculates the result string by the following actions:
 !  searches for occurences of target in string, and replaces these with
 !  substring. if back present with value true search is backward otherwise
 !  search is done forward. if every present with value true all occurences
 !  of target in string are replaced, otherwise only the first found is
 !  replaced. if target is not found the result is the same as string.
 LOGICAL                        :: dir_switch, rep_search
 CHARACTER,DIMENSION(:),POINTER :: work,temp
 INTEGER                        :: ls,lt,lsub,ipos,ipow
 ls = LEN(string); lt = LEN(target); lsub = LEN(substring)
 IF(lt==0)THEN
   WRITE(*,*) " Zero length target in REPLACE"
   STOP
 ENDIF
 ALLOCATE(work(1:ls)); work = string%chars
 IF( PRESENT(back) )THEN
   dir_switch = back
 ELSE
   dir_switch = .FALSE.
 ENDIF
 IF( PRESENT(every) )THEN
   rep_search = every
 ELSE
   rep_search = .FALSE.
 ENDIF
 IF( dir_switch )THEN ! backwards search
   ipos = ls-lt+1
   DO
     IF( ipos < 1 )EXIT ! search past start of string
     ! test for occurance of target in string at this position
     IF( ALL(string%chars(ipos:ipos+lt-1) == target%chars) )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipos-1) = work(1:ipos-1)
       temp(ipos:ipos+lsub-1) = substring%chars
       temp(ipos+lsub:) = work(ipos+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos-lt+1
     ENDIF
     ipos=ipos-1
   ENDDO
 ELSE ! forward search
   ipos = 1; ipow = 1
   DO
     IF( ipos > ls-lt+1 )EXIT ! search past end of string
     ! test for occurance of target in string at this position
     IF( ALL(string%chars(ipos:ipos+lt-1) == target%chars) )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipow-1) = work(1:ipow-1)
       temp(ipow:ipow+lsub-1) = substring%chars
       temp(ipow+lsub:) = work(ipow+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos+lt-1; ipow = ipow+lsub-1
     ENDIF
     ipos=ipos+1; ipow=ipow+1
   ENDDO
 ENDIF
 replace_sss%chars => work
ENDFUNCTION replace_sss

FUNCTION replace_ssc(string,target,substring,every,back)
 type(VARYING_STRING)            :: replace_ssc
 type(VARYING_STRING),INTENT(IN) :: string,target
 CHARACTER(LEN=*),INTENT(IN)     :: substring
 LOGICAL,INTENT(IN),OPTIONAL     :: every,back
 !  calculates the result string by the following actions:
 !  searches for occurences of target in string, and replaces these with
 !  substring. if back present with value true search is backward otherwise
 !  search is done forward. if every present with value true all occurences
 !  of target in string are replaced, otherwise only the first found is
 !  replaced. if target is not found the result is the same as string.
 LOGICAL                        :: dir_switch, rep_search
 CHARACTER,DIMENSION(:),POINTER :: work,temp
 INTEGER                        :: ls,lt,lsub,ipos,ipow,i
 ls = LEN(string); lt = LEN(target); lsub = LEN(substring)
 IF(lt==0)THEN
   WRITE(*,*) " Zero length target in REPLACE"
   STOP
 ENDIF
 ALLOCATE(work(1:ls)); work = string%chars
 IF( PRESENT(back) )THEN
   dir_switch = back
 ELSE
   dir_switch = .FALSE.
 ENDIF
 IF( PRESENT(every) )THEN
   rep_search = every
 ELSE
   rep_search = .FALSE.
 ENDIF
 IF( dir_switch )THEN ! backwards search
   ipos = ls-lt+1
   DO
     IF( ipos < 1 )EXIT ! search past start of string
     ! test for occurance of target in string at this position
     IF( ALL(string%chars(ipos:ipos+lt-1) == target%chars) )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipos-1) = work(1:ipos-1)
       DO i=1,lsub
         temp(i+ipos-1) = substring(i:i)
       ENDDO
       temp(ipos+lsub:) = work(ipos+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos-lt+1
     ENDIF
     ipos=ipos-1
   ENDDO
 ELSE ! forward search
   ipos = 1; ipow = 1
   DO
     IF( ipos > ls-lt+1 )EXIT ! search past end of string
     ! test for occurance of target in string at this position
     IF( ALL(string%chars(ipos:ipos+lt-1) == target%chars) )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipow-1) = work(1:ipow-1)
       DO i=1,lsub
         temp(i+ipow-1) = substring(i:i)
       ENDDO
       temp(ipow+lsub:) = work(ipow+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos+lt-1; ipow = ipow+lsub-1
     ENDIF
     ipos=ipos+1; ipow=ipow+1
   ENDDO
 ENDIF
 replace_ssc%chars => work
ENDFUNCTION replace_ssc

FUNCTION replace_scs(string,target,substring,every,back)
 type(VARYING_STRING)            :: replace_scs
 type(VARYING_STRING),INTENT(IN) :: string,substring
 CHARACTER(LEN=*),INTENT(IN)     :: target
 LOGICAL,INTENT(IN),OPTIONAL     :: every,back
 !  calculates the result string by the following actions:
 !  searches for occurences of target in string, and replaces these with
 !  substring. if back present with value true search is backward otherwise
 !  search is done forward. if every present with value true all accurences
 !  of target in string are replaced, otherwise only the first found is
 !  replaced. if target is not found the result is the same as string.
 LOGICAL                        :: dir_switch, rep_search
 CHARACTER,DIMENSION(:),POINTER :: work,temp,tget
 INTEGER                        :: ls,lt,lsub,ipos,ipow,i
 ls = LEN(string); lt = LEN(target); lsub = LEN(substring)
 IF(lt==0)THEN
   WRITE(*,*) " Zero length target in REPLACE"
   STOP
 ENDIF
 ALLOCATE(work(1:ls)); work = string%chars
 ALLOCATE(tget(1:lt))
 DO i=1,lt
   tget(i) = target(i:i)
 ENDDO
 IF( PRESENT(back) )THEN
   dir_switch = back
 ELSE
   dir_switch = .FALSE.
 ENDIF
 IF( PRESENT(every) )THEN
   rep_search = every
 ELSE
   rep_search = .FALSE.
 ENDIF
 IF( dir_switch )THEN ! backwards search
   ipos = ls-lt+1
   DO
     IF( ipos < 1 )EXIT ! search past start of string
     ! test for occurance of target in string at this position
     IF( ALL(string%chars(ipos:ipos+lt-1) == tget) )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipos-1) = work(1:ipos-1)
       temp(ipos:ipos+lsub-1) = substring%chars
       temp(ipos+lsub:) = work(ipos+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos-lt+1
     ENDIF
     ipos=ipos-1
   ENDDO
 ELSE ! forward search
   ipos = 1; ipow = 1
   DO
     IF( ipos > ls-lt+1 )EXIT ! search past end of string
     ! test for occurance of target in string at this position
     IF( ALL(string%chars(ipos:ipos+lt-1) == tget) )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipow-1) = work(1:ipow-1)
       temp(ipow:ipow+lsub-1) = substring%chars
       temp(ipow+lsub:) = work(ipow+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos+lt-1; ipow = ipow+lsub-1
     ENDIF
     ipos=ipos+1; ipow=ipow+1
   ENDDO
 ENDIF
 replace_scs%chars => work
ENDFUNCTION replace_scs

FUNCTION replace_scc(string,target,substring,every,back)
 type(VARYING_STRING)            :: replace_scc
 type(VARYING_STRING),INTENT(IN) :: string
 CHARACTER(LEN=*),INTENT(IN)     :: target,substring
 LOGICAL,INTENT(IN),OPTIONAL     :: every,back
 !  calculates the result string by the following actions:
 !  searches for occurences of target in string, and replaces these with
 !  substring. if back present with value true search is backward otherwise
 !  search is done forward. if every present with value true all accurences
 !  of target in string are replaced, otherwise only the first found is
 !  replaced. if target is not found the result is the same as string.
 LOGICAL                        :: dir_switch, rep_search
 CHARACTER,DIMENSION(:),POINTER :: work,temp,tget
 INTEGER                        :: ls,lt,lsub,ipos,ipow,i
 ls = LEN(string); lt = LEN(target); lsub = LEN(substring)
 IF(lt==0)THEN
   WRITE(*,*) " Zero length target in REPLACE"
   STOP
 ENDIF
 ALLOCATE(work(1:ls)); work = string%chars
 ALLOCATE(tget(1:lt))
 DO i=1,lt
   tget(i) = target(i:i)
 ENDDO
 IF( PRESENT(back) )THEN
   dir_switch = back
 ELSE
   dir_switch = .FALSE.
 ENDIF
 IF( PRESENT(every) )THEN
   rep_search = every
 ELSE
   rep_search = .FALSE.
 ENDIF
 IF( dir_switch )THEN ! backwards search
   ipos = ls-lt+1
   DO
     IF( ipos < 1 )EXIT ! search past start of string
     ! test for occurance of target in string at this position
     IF( ALL(string%chars(ipos:ipos+lt-1) == tget) )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipos-1) = work(1:ipos-1)
       DO i=1,lsub
         temp(i+ipos-1) = substring(i:i)
       ENDDO
       temp(ipos+lsub:) = work(ipos+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos-lt+1
     ENDIF
     ipos=ipos-1
   ENDDO
 ELSE ! forward search
   ipos = 1; ipow = 1
   DO
     IF( ipos > ls-lt+1 )EXIT ! search past end of string
     ! test for occurance of target in string at this position
     IF( ALL(string%chars(ipos:ipos+lt-1) == tget) )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipow-1) = work(1:ipow-1)
       DO i=1,lsub
         temp(i+ipow-1) = substring(i:i)
       ENDDO
       temp(ipow+lsub:) = work(ipow+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos+lt-1; ipow = ipow+lsub-1
     ENDIF
     ipos=ipos+1; ipow=ipow+1
   ENDDO
 ENDIF
 replace_scc%chars => work
ENDFUNCTION replace_scc

FUNCTION replace_css(string,target,substring,every,back)
 type(VARYING_STRING)            :: replace_css
 CHARACTER(LEN=*),INTENT(IN)     :: string
 type(VARYING_STRING),INTENT(IN) :: target,substring
 LOGICAL,INTENT(IN),OPTIONAL     :: every,back
 ! calculates the result string by the following actions:
 ! searches for occurences of target in string, and replaces these with
 ! substring. if back present with value true search is backward otherwise
 ! search is done forward. if every present with value true all accurences
 ! of target in string are replaced, otherwise only the first found is
 ! replaced. if target is not found the result is the same as string.
 LOGICAL                        :: dir_switch, rep_search
 CHARACTER,DIMENSION(:),POINTER :: work,temp,str
 INTEGER                        :: ls,lt,lsub,ipos,ipow,i
 ls = LEN(string); lt = LEN(target); lsub = LEN(substring)
 IF(lt==0)THEN
   WRITE(*,*) " Zero length target in REPLACE"
   STOP
 ENDIF
 ALLOCATE(work(1:ls)); ALLOCATE(str(1:ls))
 DO i=1,ls
   str(i) = string(i:i)
 ENDDO
 work = str
 IF( PRESENT(back) )THEN
   dir_switch = back
 ELSE
   dir_switch = .FALSE.
 ENDIF
 IF( PRESENT(every) )THEN
   rep_search = every
 ELSE
   rep_search = .FALSE.
 ENDIF
 IF( dir_switch )THEN ! backwards search
   ipos = ls-lt+1
   DO
     IF( ipos < 1 )EXIT ! search past start of string
     ! test for occurance of target in string at this position
     IF( ALL(str(ipos:ipos+lt-1) == target%chars) )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipos-1) = work(1:ipos-1)
       temp(ipos:ipos+lsub-1) = substring%chars
       temp(ipos+lsub:) = work(ipos+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos-lt+1
     ENDIF
     ipos=ipos-1
   ENDDO
 ELSE ! forward search
   ipos = 1; ipow = 1
   DO
     IF( ipos > ls-lt+1 )EXIT ! search past end of string
     ! test for occurance of target in string at this position
     IF( ALL(str(ipos:ipos+lt-1) == target%chars) )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipow-1) = work(1:ipow-1)
       temp(ipow:ipow+lsub-1) = substring%chars
       temp(ipow+lsub:) = work(ipow+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos+lt-1; ipow = ipow+lsub-1
     ENDIF
     ipos=ipos+1; ipow=ipow+1
   ENDDO
 ENDIF
 replace_css%chars => work
ENDFUNCTION replace_css

FUNCTION replace_csc(string,target,substring,every,back)
 type(VARYING_STRING)            :: replace_csc
 type(VARYING_STRING),INTENT(IN) :: target
 CHARACTER(LEN=*),INTENT(IN)     :: string,substring
 LOGICAL,INTENT(IN),OPTIONAL     :: every,back
 !  calculates the result string by the following actions:
 !  searches for occurences of target in string, and replaces these with
 !  substring. if back present with value true search is backward otherwise
 !  search is done forward. if every present with value true all accurences
 !  of target in string are replaced, otherwise only the first found is
 !  replaced. if target is not found the result is the same as string.
 LOGICAL                        :: dir_switch, rep_search
 CHARACTER,DIMENSION(:),POINTER :: work,temp,str
 INTEGER                        :: ls,lt,lsub,ipos,ipow,i
 ls = LEN(string); lt = LEN(target); lsub = LEN(substring)
 IF(lt==0)THEN
   WRITE(*,*) " Zero length target in REPLACE"
   STOP
 ENDIF
 ALLOCATE(work(1:ls)); ALLOCATE(str(1:ls))
 DO i=1,ls
   str(i) = string(i:i)
 ENDDO
 work = str
 IF( PRESENT(back) )THEN
   dir_switch = back
 ELSE
   dir_switch = .FALSE.
 ENDIF
 IF( PRESENT(every) )THEN
   rep_search = every
 ELSE
   rep_search = .FALSE.
 ENDIF
 IF( dir_switch )THEN ! backwards search
   ipos = ls-lt+1
   DO
     IF( ipos < 1 )EXIT ! search past start of string
     ! test for occurance of target in string at this position
     IF( ALL(str(ipos:ipos+lt-1) == target%chars) )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipos-1) = work(1:ipos-1)
       DO i=1,lsub
         temp(i+ipos-1) = substring(i:i)
       ENDDO
       temp(ipos+lsub:) = work(ipos+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos-lt+1
     ENDIF
     ipos=ipos-1
   ENDDO
 ELSE ! forward search
   ipos = 1; ipow = 1
   DO
     IF( ipos > ls-lt+1 )EXIT ! search past end of string
     ! test for occurance of target in string at this position
     IF( ALL(str(ipos:ipos+lt-1) == target%chars) )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipow-1) = work(1:ipow-1)
       DO i=1,lsub
         temp(i+ipow-1) = substring(i:i)
       ENDDO
       temp(ipow+lsub:) = work(ipow+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos+lt-1; ipow = ipow+lsub-1
     ENDIF
     ipos=ipos+1; ipow=ipow+1
   ENDDO
 ENDIF
 replace_csc%chars => work
ENDFUNCTION replace_csc

FUNCTION replace_ccs(string,target,substring,every,back)
 type(VARYING_STRING)            :: replace_ccs
 type(VARYING_STRING),INTENT(IN) :: substring
 CHARACTER(LEN=*),INTENT(IN)     :: string,target
 LOGICAL,INTENT(IN),OPTIONAL     :: every,back
 !  calculates the result string by the following actions:
 !  searches for occurences of target in string, and replaces these with
 !  substring. if back present with value true search is backward otherwise
 !  search is done forward. if every present with value true all accurences
 !  of target in string are replaced, otherwise only the first found is
 !  replaced. if target is not found the result is the same as string.
 LOGICAL                        :: dir_switch, rep_search
 CHARACTER,DIMENSION(:),POINTER :: work,temp
 INTEGER                        :: ls,lt,lsub,ipos,ipow,i
 ls = LEN(string); lt = LEN(target); lsub = LEN(substring)
 IF(lt==0)THEN
   WRITE(*,*) " Zero length target in REPLACE"
   STOP
 ENDIF
 ALLOCATE(work(1:ls))
 DO i=1,ls
   work(i) = string(i:i)
 ENDDO
 IF( PRESENT(back) )THEN
   dir_switch = back
 ELSE
   dir_switch = .FALSE.
 ENDIF
 IF( PRESENT(every) )THEN
   rep_search = every
 ELSE
   rep_search = .FALSE.
 ENDIF
 IF( dir_switch )THEN ! backwards search
   ipos = ls-lt+1
   DO
     IF( ipos < 1 )EXIT ! search past start of string
     ! test for occurance of target in string at this position
     IF( string(ipos:ipos+lt-1) == target )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipos-1) = work(1:ipos-1)
       temp(ipos:ipos+lsub-1) = substring%chars
       temp(ipos+lsub:) = work(ipos+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos-lt+1
     ENDIF
     ipos=ipos-1
   ENDDO
 ELSE ! forward search
   ipos = 1; ipow = 1
   DO
     IF( ipos > ls-lt+1 )EXIT ! search past end of string
     ! test for occurance of target in string at this position
     IF( string(ipos:ipos+lt-1) == target )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipow-1) = work(1:ipow-1)
       temp(ipow:ipow+lsub-1) = substring%chars
       temp(ipow+lsub:) = work(ipow+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos+lt-1; ipow = ipow+lsub-1
     ENDIF
     ipos=ipos+1; ipow=ipow+1
   ENDDO
 ENDIF
 replace_ccs%chars => work
ENDFUNCTION replace_ccs

FUNCTION replace_ccc(string,target,substring,every,back)
 type(VARYING_STRING)            :: replace_ccc
 CHARACTER(LEN=*),INTENT(IN)     :: string,target,substring
 LOGICAL,INTENT(IN),OPTIONAL     :: every,back
 !  calculates the result string by the following actions:
 !  searches for occurences of target in string, and replaces these with
 !  substring. if back present with value true search is backward otherwise
 !  search is done forward. if every present with value true all accurences
 !  of target in string are replaced, otherwise only the first found is
 !  replaced. if target is not found the result is the same as string.
 LOGICAL                        :: dir_switch, rep_search
 CHARACTER,DIMENSION(:),POINTER :: work,temp
 INTEGER                        :: ls,lt,lsub,ipos,ipow,i
 ls = LEN(string); lt = LEN(target); lsub = LEN(substring)
 IF(lt==0)THEN
   WRITE(*,*) " Zero length target in REPLACE"
   STOP
 ENDIF
 ALLOCATE(work(1:ls))
 DO i=1,ls
   work(i) = string(i:i)
 ENDDO
 IF( PRESENT(back) )THEN
   dir_switch = back
 ELSE
   dir_switch = .FALSE.
 ENDIF
 IF( PRESENT(every) )THEN
   rep_search = every
 ELSE
   rep_search = .FALSE.
 ENDIF
 IF( dir_switch )THEN ! backwards search
   ipos = ls-lt+1
   DO
     IF( ipos < 1 )EXIT ! search past start of string
     ! test for occurance of target in string at this position
     IF( string(ipos:ipos+lt-1) == target )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipos-1) = work(1:ipos-1)
       DO i=1,lsub
         temp(i+ipos-1) = substring(i:i)
       ENDDO
       temp(ipos+lsub:) = work(ipos+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos-lt+1
     ENDIF
     ipos=ipos-1
   ENDDO
 ELSE ! forward search
   ipos = 1; ipow = 1
   DO
     IF( ipos > ls-lt+1 )EXIT ! search past end of string
     ! test for occurance of target in string at this position
     IF( string(ipos:ipos+lt-1) == target )THEN
       ! match found allocate space for string with this occurance of
       ! target replaced by substring
       ALLOCATE(temp(1:SIZE(work)+lsub-lt))
       ! copy work into temp replacing this occurance of target by
       ! substring
       temp(1:ipow-1) = work(1:ipow-1)
       DO i=1,lsub
         temp(i+ipow-1) = substring(i:i)
       ENDDO
       temp(ipow+lsub:) = work(ipow+lt:)
       work => temp ! make new version of work point at the temp space
       IF(.NOT.rep_search)EXIT ! exit if only first replacement wanted
       ! move search and replacement positions over the effected positions
       ipos = ipos+lt-1; ipow = ipow+lsub-1
     ENDIF
     ipos=ipos+1; ipow=ipow+1
   ENDDO
 ENDIF
 replace_ccc%chars => work
ENDFUNCTION replace_ccc

!----- Remove procedures ----------------------------------------------------!
FUNCTION remove_s(string,start,finish)
 type(VARYING_STRING)            :: remove_s
 type(VARYING_STRING),INTENT(IN) :: string
 INTEGER,INTENT(IN),OPTIONAL     :: start
 INTEGER,INTENT(IN),OPTIONAL     :: finish
 !  returns as result the string produced by the actions
 !  removes the characters between start and finish from string reducing it in 
 !  size by MAX(0,ABS(finish-start+1)) 
 !  if start < 1 or is missing then assumes start=1 
 !  if finish > LEN(string) or is missing then assumes finish=LEN(string) 
 CHARACTER,DIMENSION(:),POINTER :: arg_str
 INTEGER                        :: is,if,ls
 ls = LEN(string)
 IF (PRESENT(start)) THEN
   is = MAX(1,start)
 ELSE
   is = 1
 ENDIF
 IF (PRESENT(finish)) THEN
   if = MIN(ls,finish)
 ELSE
   if = ls
 ENDIF
 IF( if < is ) THEN  ! zero characters to be removed, string is unchanged
   ALLOCATE(arg_str(1:ls))
   arg_str = string%chars
 ELSE
   ALLOCATE(arg_str(1:ls-if+is-1) )
   arg_str(1:is-1) = string%chars(1:is-1)
   arg_str(is:) = string%chars(if+1:)
 ENDIF
 remove_s%chars => arg_str
ENDFUNCTION remove_s
  
FUNCTION remove_c(string,start,finish)
 type(VARYING_STRING)        :: remove_c
 CHARACTER(LEN=*),INTENT(IN) :: string
 INTEGER,INTENT(IN),OPTIONAL :: start
 INTEGER,INTENT(IN),OPTIONAL :: finish
 !  returns as result the string produced by the actions
 !  removes the characters between start and finish from string reducing it in 
 !  size by MAX(0,ABS(finish-start+1)) 
 !  if start < 1 or is missing then assumes start=1 
 !  if finish > LEN(string) or is missing then assumes finish=LEN(string) 
 CHARACTER,DIMENSION(:),POINTER :: arg_str
 INTEGER                        :: is,if,ls,i
 ls = LEN(string)
 IF (PRESENT(start)) THEN
   is = MAX(1,start)
 ELSE
   is = 1
 ENDIF
 IF (PRESENT(finish)) THEN
   if = MIN(ls,finish)
 ELSE
   if = ls
 ENDIF
 IF( if < is ) THEN  ! zero characters to be removed, string is unchanged
   ALLOCATE(arg_str(1:ls))
   DO i=1,ls
     arg_str(i) = string(i:i)
   ENDDO
 ELSE
   ALLOCATE(arg_str(1:ls-if+is-1) )
   DO i=1,is-1
     arg_str(i) = string(i:i)
   ENDDO
   DO i=is,ls-if+is-1
     arg_str(i) = string(i-is+if+1:i-is+if+1)
   ENDDO
 ENDIF
 remove_c%chars => arg_str
ENDFUNCTION remove_c
  
!----- Extract procedures ---------------------------------------------------!
 FUNCTION extract_s(string,start,finish) 
  type(VARYING_STRING),INTENT(IN) :: string 
  INTEGER,INTENT(IN),OPTIONAL     :: start     
  INTEGER,INTENT(IN),OPTIONAL     :: finish     
  type(VARYING_STRING)            :: extract_s 
  ! extracts the characters between start and finish from string  and 
  ! delivers these as the result of the function, string is unchanged 
  ! if start < 1 or is missing then it is treated as 1 
  ! if finish > LEN(string) or is missing then it is treated as LEN(string) 
  INTEGER                         :: is,if 
  IF (PRESENT(start)) THEN  
     is = MAX(1,start) 
  ELSE 
     is = 1 
  ENDIF 
  IF (PRESENT(finish)) THEN  
     if = MIN(LEN(string),finish) 
  ELSE 
     if = LEN(string) 
  ENDIF 
  ALLOCATE(extract_s%chars(1:if-is+1)) 
  extract_s%chars = string%chars(is:if)
 ENDFUNCTION extract_s 
  
 FUNCTION extract_c(string,start,finish)
  CHARACTER(LEN=*),INTENT(IN) :: string 
  INTEGER,INTENT(IN),OPTIONAL :: start   
  INTEGER,INTENT(IN),OPTIONAL :: finish  
  type(VARYING_STRING)        :: extract_c 
  ! extracts the characters between start and finish from character string and 
  ! delivers these as the result of the function, string is unchanged 
  ! if start < 1 or is missing then it is treated as 1 
  ! if finish > LEN(string) or is missing then it is treated as LEN(string) 
  INTEGER                      :: is,if,i 
  IF (PRESENT(start)) THEN    
     is = MAX(1,start) 
  ELSE 
     is = 1 
  ENDIF 
  IF (PRESENT(finish)) THEN  
     if = MIN(LEN(string),finish) 
  ELSE 
     if = LEN(string) 
  ENDIF 
  ALLOCATE(extract_c%chars(1:if-is+1)) 
  DO i=is,if 
    extract_c%chars(i-is+1) = string(i:i) 
  ENDDO 
 ENDFUNCTION extract_c 

!----- Split procedures ------------------------------------------------------!
 SUBROUTINE split_s(string,word,set,separator,back)
  type(VARYING_STRING),INTENT(INOUT)        :: string
  type(VARYING_STRING),INTENT(OUT)          :: word
  type(VARYING_STRING),INTENT(IN)           :: set
  type(VARYING_STRING),INTENT(OUT),OPTIONAL :: separator
  LOGICAL,INTENT(IN),OPTIONAL               :: back
  ! splits the input string at the first(last) character in set
  ! returns the leading(trailing) substring in word and the trailing(leading)
  ! substring in string. The search is done in the forward or backward
  ! direction depending on back. If separator is present, the actual separator
  ! character found is returned in separator.
  ! If no character in set is found string and separator are returned as
  ! zero length and the whole input string is returned in word.
  LOGICAL                         :: dir_switch 
  INTEGER                         :: ls,tpos
  ls = LEN(string)
  IF( PRESENT(back) )THEN 
    dir_switch = back 
  ELSE 
    dir_switch = .FALSE. 
  ENDIF 
  IF(dir_switch)THEN ! backwards search 
    DO tpos = ls,1,-1
       IF(ANY(string%chars(tpos) == set%chars))EXIT
    ENDDO
    word%chars => string%chars(tpos+1:ls)
    IF(PRESENT(separator))THEN
      IF(tpos==0)THEN
        separator = ""
      ELSE
        separator%chars => string%chars(tpos:tpos)
      ENDIF
    ENDIF
    string%chars => string%chars(1:tpos-1)
  ELSE ! forwards search
    DO tpos =1,ls
       IF(ANY(string%chars(tpos) == set%chars))EXIT
    ENDDO
    word%chars => string%chars(1:tpos-1)
    IF(PRESENT(separator))THEN
      IF(tpos==ls+1)THEN
        separator = ""
      ELSE
        separator%chars => string%chars(tpos:tpos)
      ENDIF
    ENDIF
    string%chars => string%chars(tpos+1:ls)
  ENDIF
 ENDSUBROUTINE split_s

 SUBROUTINE split_c(string,word,set,separator,back)
  type(VARYING_STRING),INTENT(INOUT)        :: string
  type(VARYING_STRING),INTENT(OUT)          :: word
  CHARACTER(LEN=*),INTENT(IN)               :: set
  type(VARYING_STRING),INTENT(OUT),OPTIONAL :: separator
  LOGICAL,INTENT(IN),OPTIONAL               :: back
  ! splits the input string at the first(last) character in set
  ! returns the leading(trailing) substring in word and the trailing(leading)
  ! substring in string. The search is done in the forward or backward
  ! direction depending on back. If separator is present, the actual separator
  ! character found is returned in separator.
  ! If no character in set is found string and separator are returned as
  ! zero length and the whole input string is returned in word.
  LOGICAL                         :: dir_switch 
  INTEGER                         :: ls,tpos,lset,i
  ls = LEN(string); lset = LEN(set)
  IF( PRESENT(back) )THEN 
    dir_switch = back 
  ELSE 
    dir_switch = .FALSE. 
  ENDIF 
  IF(dir_switch)THEN ! backwards search 
    BSEARCH:DO tpos = ls,1,-1
       DO i=1,lset
         IF(string%chars(tpos) == set(i:i))EXIT BSEARCH
       ENDDO
    ENDDO BSEARCH
    word%chars => string%chars(tpos+1:ls)
    IF(PRESENT(separator))THEN
      IF(tpos==0)THEN
        separator = ""
      ELSE
        separator%chars => string%chars(tpos:tpos)
      ENDIF
    ENDIF
    string%chars => string%chars(1:tpos-1)
  ELSE ! forwards search
    FSEARCH:DO tpos =1,ls
       DO i=1,lset
         IF(string%chars(tpos) == set(i:i))EXIT FSEARCH
       ENDDO
    ENDDO FSEARCH
    word%chars => string%chars(1:tpos-1)
    IF(PRESENT(separator))THEN
      IF(tpos==ls+1)THEN
        separator = ""
      ELSE
        separator%chars => string%chars(tpos:tpos)
      ENDIF
    ENDIF
    string%chars => string%chars(tpos+1:ls)
  ENDIF
 ENDSUBROUTINE split_c

!----- INDEX procedures ------------------------------------------------------!
 FUNCTION index_ss(string,substring,back) 
  type(VARYING_STRING),INTENT(IN) :: string,substring 
  LOGICAL,INTENT(IN),OPTIONAL     :: back 
  INTEGER                         :: index_ss 
  ! returns the starting position in string of the substring 
  ! scanning from the front or back depending on the logical argument back 
  LOGICAL                         :: dir_switch 
  INTEGER                         :: ls,lsub,i 
  ls = LEN(string); lsub = LEN(substring) 
  IF( PRESENT(back) )THEN 
    dir_switch = back 
  ELSE 
    dir_switch = .FALSE. 
  ENDIF 
  IF(dir_switch)THEN ! backwards search 
    DO i = ls-lsub+1,1,-1 
      IF( ALL(string%chars(i:i+lsub-1) == substring%chars) )THEN 
        index_ss = i 
        RETURN 
      ENDIF 
    ENDDO 
    index_ss = 0 
  ELSE ! forward search 
    DO i = 1,ls-lsub+1 
      IF( ALL(string%chars(i:i+lsub-1) == substring%chars) )THEN 
        index_ss = i 
        RETURN 
      ENDIF 
    ENDDO 
    index_ss = 0 
  ENDIF 
 ENDFUNCTION index_ss 
  
 FUNCTION index_sc(string,substring,back)
  type(VARYING_STRING),INTENT(IN) :: string 
  CHARACTER(LEN=*),INTENT(IN)     :: substring 
  LOGICAL,INTENT(IN),OPTIONAL     :: back 
  INTEGER                         :: index_sc 
  ! returns the starting position in string of the substring 
  ! scanning from the front or back depending on the logical argument back 
  LOGICAL                         :: dir_switch,matched 
  INTEGER                         :: ls,lsub,i,j 
  ls = LEN(string); lsub = LEN(substring) 
  IF( PRESENT(back) )THEN 
    dir_switch = back 
  ELSE 
    dir_switch = .FALSE. 
  ENDIF 
  IF (dir_switch) THEN ! backwards search 
    DO i = ls-lsub+1,1,-1 
      matched = .TRUE. 
      DO j = 1,lsub 
        IF( string%chars(i+j-1) /= substring(j:j) )THEN 
          matched = .FALSE. 
          EXIT 
        ENDIF 
      ENDDO 
      IF( matched )THEN 
        index_sc = i 
        RETURN 
      ENDIF 
    ENDDO 
    index_sc = 0 
  ELSE ! forward search 
    DO i = 1,ls-lsub+1 
      matched = .TRUE. 
      DO j = 1,lsub 
        IF( string%chars(i+j-1) /= substring(j:j) )THEN 
          matched = .FALSE. 
          EXIT 
        ENDIF 
      ENDDO 
      IF( matched )THEN 
        index_sc = i 
        RETURN 
      ENDIF 
    ENDDO 
    index_sc = 0 
  ENDIF 
 ENDFUNCTION index_sc 
  
 FUNCTION index_cs(string,substring,back)
  CHARACTER(LEN=*),INTENT(IN)     :: string 
  type(VARYING_STRING),INTENT(IN) :: substring 
  LOGICAL,INTENT(IN),OPTIONAL     :: back 
  INTEGER                         :: index_cs 
  ! returns the starting position in string of the substring 
  ! scanning from the front or back depending on the logical argument back 
  LOGICAL                         :: dir_switch,matched 
  INTEGER                         :: ls,lsub,i,j 
  ls = LEN(string); lsub = LEN(substring) 
  IF( PRESENT(back) )THEN 
    dir_switch = back 
  ELSE 
    dir_switch = .FALSE. 
  ENDIF 
  IF(dir_switch)THEN ! backwards search 
    DO i = ls-lsub+1,1,-1 
      matched = .TRUE. 
      DO j = 1,lsub 
        IF( string(i+j-1:i+j-1) /= substring%chars(j) )THEN 
          matched = .FALSE. 
          EXIT 
        ENDIF 
      ENDDO 
      IF( matched )THEN 
        index_cs = i 
        RETURN 
      ENDIF 
    ENDDO 
    index_cs = 0 
  ELSE ! forward search 
    DO i = 1,ls-lsub+1 
      matched = .TRUE. 
      DO j = 1,lsub 
        IF( string(i+j-1:i+j-1) /= substring%chars(j) )THEN 
          matched = .FALSE. 
          EXIT 
        ENDIF 
      ENDDO 
      IF( matched )THEN 
        index_cs = i 
        RETURN 
      ENDIF 
    ENDDO 
    index_cs = 0 
  ENDIF 
 ENDFUNCTION index_cs 
  
!----- SCAN procedures ------------------------------------------------------!
 FUNCTION scan_ss(string,set,back) 
  type(VARYING_STRING),INTENT(IN) :: string,set 
  LOGICAL,INTENT(IN),OPTIONAL     :: back 
  INTEGER                         :: scan_ss 
  ! returns the first position in string occupied by a character from 
  ! the characters in set, scanning is forward or backwards depending on back 
  LOGICAL                         :: dir_switch 
  INTEGER                         :: ls,i 
  ls = LEN(string) 
  IF( PRESENT(back) )THEN 
    dir_switch = back 
  ELSE 
    dir_switch = .FALSE. 
  ENDIF 
  IF(dir_switch)THEN ! backwards search 
    DO i = ls,1,-1 
      IF( ANY( set%chars == string%chars(i) ) )THEN 
        scan_ss = i 
        RETURN 
      ENDIF 
    ENDDO 
    scan_ss = 0 
  ELSE ! forward search 
    DO i = 1,ls 
      IF( ANY( set%chars == string%chars(i) ) )THEN 
        scan_ss = i 
        RETURN 
      ENDIF 
    ENDDO 
    scan_ss = 0 
  ENDIF 
 ENDFUNCTION scan_ss 
  
 FUNCTION scan_sc(string,set,back)
  type(VARYING_STRING),INTENT(IN) :: string 
  CHARACTER(LEN=*),INTENT(IN)     :: set 
  LOGICAL,INTENT(IN),OPTIONAL     :: back 
  INTEGER                         :: scan_sc 
  ! returns the first position in string occupied by a character from 
  ! the characters in set, scanning is forward or backwards depending on back 
  LOGICAL                         :: dir_switch,matched 
  INTEGER                         :: ls,i,j 
  ls = LEN(string) 
  IF( PRESENT(back) )THEN 
    dir_switch = back 
  ELSE 
    dir_switch = .FALSE. 
  ENDIF 
  IF(dir_switch)THEN ! backwards search 
    DO i = ls,1,-1 
      matched = .FALSE. 
      DO j = 1,LEN(set) 
        IF( string%chars(i) == set(j:j) )THEN 
          matched = .TRUE. 
          EXIT 
        ENDIF 
      ENDDO 
      IF( matched )THEN 
        scan_sc = i 
        RETURN 
      ENDIF 
    ENDDO 
    scan_sc = 0 
  ELSE ! forward search 
    DO i = 1,ls 
      matched = .FALSE. 
      DO j = 1,LEN(set) 
        IF( string%chars(i) == set(j:j) )THEN 
          matched = .TRUE. 
          EXIT 
        ENDIF 
      ENDDO 
      IF( matched )THEN 
        scan_sc = i 
        RETURN 
      ENDIF 
    ENDDO 
    scan_sc = 0 
  ENDIF 
 ENDFUNCTION scan_sc 
  
 FUNCTION scan_cs(string,set,back)
  CHARACTER(LEN=*),INTENT(IN)     :: string 
  type(VARYING_STRING),INTENT(IN) :: set 
  LOGICAL,INTENT(IN),OPTIONAL     :: back 
  INTEGER                         :: scan_cs 
  ! returns the first position in character string occupied by a character from 
  ! the characters in set, scanning is forward or backwards depending on back 
  LOGICAL                         :: dir_switch,matched 
  INTEGER                         :: ls,i,j 
  ls = LEN(string) 
  IF( PRESENT(back) )THEN 
    dir_switch = back 
  ELSE 
    dir_switch = .FALSE. 
  ENDIF 
  IF(dir_switch)THEN ! backwards search 
    DO i = ls,1,-1 
      matched = .FALSE. 
      DO j = 1,LEN(set) 
        IF( string(i:i) == set%chars(j) )THEN 
          matched = .TRUE. 
          EXIT 
        ENDIF 
      ENDDO 
      IF( matched )THEN 
        scan_cs = i 
        RETURN 
      ENDIF 
    ENDDO 
    scan_cs = 0 
  ELSE ! forward search 
    DO i = 1,ls 
      matched = .FALSE. 
      DO j = 1,LEN(set) 
        IF( string(i:i) == set%chars(j) )THEN 
          matched = .TRUE. 
          EXIT 
        ENDIF 
      ENDDO 
      IF( matched )THEN 
        scan_cs = i 
        RETURN 
      ENDIF 
    ENDDO 
    scan_cs = 0 
  ENDIF 
 ENDFUNCTION scan_cs 
  
!----- VERIFY procedures ----------------------------------------------------!
 FUNCTION verify_ss(string,set,back) 
  type(VARYING_STRING),INTENT(IN) :: string,set 
  LOGICAL,INTENT(IN),OPTIONAL     :: back 
  INTEGER                         :: verify_ss 
  ! returns the first position in string not occupied by a character from 
  ! the characters in set, scanning is forward or backwards depending on back 
  LOGICAL                     :: dir_switch 
  INTEGER                     :: ls,i 
  ls = LEN(string) 
  IF( PRESENT(back) )THEN 
    dir_switch = back 
  ELSE 
    dir_switch = .FALSE. 
  ENDIF 
  IF(dir_switch)THEN ! backwards search 
    DO i = ls,1,-1 
      IF( .NOT.(ANY( set%chars == string%chars(i) )) )THEN 
        verify_ss = i 
        RETURN 
      ENDIF 
    ENDDO 
    verify_ss = 0 
  ELSE ! forward search 
    DO i = 1,ls 
      IF( .NOT.(ANY( set%chars == string%chars(i) )) )THEN 
        verify_ss = i 
        RETURN 
      ENDIF 
    ENDDO 
    verify_ss = 0 
  ENDIF 
 ENDFUNCTION verify_ss 
  
 FUNCTION verify_sc(string,set,back)
  type(VARYING_STRING),INTENT(IN) :: string 
  CHARACTER(LEN=*),INTENT(IN)     :: set 
  LOGICAL,INTENT(IN),OPTIONAL     :: back 
  INTEGER                         :: verify_sc 
  ! returns the first position in string not occupied by a character from 
  ! the characters in set, scanning is forward or backwards depending on back 
  LOGICAL                     :: dir_switch
  INTEGER                     :: ls,i,j 
  ls = LEN(string) 
  IF( PRESENT(back) )THEN 
    dir_switch = back 
  ELSE 
    dir_switch = .FALSE. 
  ENDIF 
  IF(dir_switch)THEN ! backwards search 
    back_string_search:DO i = ls,1,-1 
      DO j = 1,LEN(set) 
        IF( string%chars(i) == set(j:j) )CYCLE back_string_search
        ! cycle string search if string character found in set
      ENDDO 
      ! string character not found in set index i is result
        verify_sc = i 
        RETURN 
    ENDDO back_string_search
    ! each string character found in set
    verify_sc = 0 
  ELSE ! forward search 
    frwd_string_search:DO i = 1,ls 
      DO j = 1,LEN(set) 
        IF( string%chars(i) == set(j:j) )CYCLE frwd_string_search
      ENDDO 
        verify_sc = i 
        RETURN 
    ENDDO frwd_string_search
    verify_sc = 0 
  ENDIF 
 ENDFUNCTION verify_sc 
  
 FUNCTION verify_cs(string,set,back)
  CHARACTER(LEN=*),INTENT(IN)     :: string 
  type(VARYING_STRING),INTENT(IN) :: set 
  LOGICAL,INTENT(IN),OPTIONAL     :: back 
  INTEGER                         :: verify_cs 
  ! returns the first position in icharacter string not occupied by a character 
  ! from the characters in set, scanning is forward or backwards depending on 
  ! back 
  LOGICAL                     :: dir_switch
  INTEGER                     :: ls,i,j 
  ls = LEN(string) 
  IF( PRESENT(back) )THEN 
    dir_switch = back 
  ELSE 
    dir_switch = .FALSE. 
  ENDIF 
  IF(dir_switch)THEN ! backwards search 
    back_string_search:DO i = ls,1,-1 
      DO j = 1,LEN(set) 
        IF( string(i:i) == set%chars(j) )CYCLE back_string_search
      ENDDO 
        verify_cs = i 
        RETURN 
    ENDDO back_string_search
    verify_cs = 0 
  ELSE ! forward search 
    frwd_string_search:DO i = 1,ls 
      DO j = 1,LEN(set) 
        IF( string(i:i) == set%chars(j) )CYCLE frwd_string_search
      ENDDO 
        verify_cs = i 
        RETURN 
    ENDDO frwd_string_search
    verify_cs = 0 
  ENDIF 
 ENDFUNCTION verify_cs 
    
!----- LEN_TRIM procedure ----------------------------------------------------!
FUNCTION len_trim_s(string) 
 type(VARYING_STRING),INTENT(IN) :: string 
 INTEGER                         :: len_trim_s 
 ! Returns the length of the string without counting trailing blanks 
 INTEGER                         :: ls,i 
 ls=LEN(string) 
 len_trim_s = 0 
 DO i = ls,1,-1 
    IF (string%chars(i) /= BLANK) THEN 
       len_trim_s = i 
       EXIT 
    ENDIF 
 ENDDO 
ENDFUNCTION len_trim_s 
  
!----- TRIM procedure -------------------------------------------------------! 
FUNCTION trim_s(string) 
 type(VARYING_STRING),INTENT(IN)  :: string 
 type(VARYING_STRING)             :: trim_s 
 ! Returns the argument string with trailing blanks removed 
 INTEGER                      :: ls,pos,i 
 ls=LEN(string) 
 pos=0 
 DO i = ls,1,-1 
    IF(string%chars(i) /= BLANK) THEN 
       pos=i 
       EXIT 
    ENDIF 
 ENDDO 
 ALLOCATE(trim_s%chars(1:pos))
 trim_s%chars(1:pos) = string%chars(1:pos) 
ENDFUNCTION trim_s 
  
!----- IACHAR procedure ------------------------------------------------------! 
FUNCTION iachar_s(string) 
 type(VARYING_STRING),INTENT(IN) :: string 
 INTEGER                         :: iachar_s 
 ! returns the position of the character string in the ISO 646 
 ! collating sequence. 
 ! string must be of length one 
 IF (LEN(string) /= 1) THEN 
    WRITE(*,*) " ERROR, argument in IACHAR not of length one" 
    STOP 
 ENDIF 
 iachar_s = IACHAR(string%chars(1)) 
ENDFUNCTION iachar_s 

!----- ICHAR procedure ------------------------------------------------------!
FUNCTION ichar_s(string) 
 type(VARYING_STRING),INTENT(IN) :: string 
 INTEGER                         :: ichar_s 
 ! returns the position of character from string in the processor collating 
 ! sequence. 
 ! string must be of length one 
 IF (LEN(string) /= 1) THEN 
    WRITE(*,*) " Argument string in ICHAR has to be of length one" 
    STOP 
 ENDIF 
 ichar_s = ICHAR(string%chars(1)) 
ENDFUNCTION ichar_s 
  
!----- ADJUSTL procedure ----------------------------------------------------! 
FUNCTION adjustl_s(string) 
 type(VARYING_STRING),INTENT(IN) :: string 
 type(VARYING_STRING)            :: adjustl_s 
 ! Returns the string adjusted to the left, removing leading blanks and 
 ! inserting trailing blanks 
 INTEGER                         :: ls,pos 
 ls=LEN(string) 
 DO pos = 1,ls 
    IF(string%chars(pos) /= blank) EXIT 
 ENDDO 
 ! pos now holds the position of the first non-blank character 
 ! or ls+1 if all characters are blank
 ALLOCATE(adjustl_s%chars(1:ls)) 
 adjustl_s%chars(1:ls-pos+1) = string%chars(pos:ls)
 adjustl_s%chars(ls-pos+2:ls) = blank
ENDFUNCTION adjustl_s 
  
!----- ADJUSTR procedure ----------------------------------------------------! 
FUNCTION adjustr_s(string) 
 type(VARYING_STRING),INTENT(IN) :: string 
 type(VARYING_STRING)            :: adjustr_s 
 ! Returns the string adjusted to the right, removing trailing blanks 
 ! and inserting leading blanks 
 INTEGER                         :: ls,pos 
 ls=LEN(string) 
 DO pos = ls,1,-1 
    IF(string%chars(pos) /= blank) EXIT 
 ENDDO 
 ! pos now holds the position of the last non-blank character
 ! or zero if all characters are blank
 ALLOCATE(adjustr_s%chars(1:ls)) 
 adjustr_s%chars(ls-pos+1:ls) = string%chars(1:pos) 
 adjustr_s%chars(1:ls-pos) = blank 
ENDFUNCTION adjustr_s 
  
ENDMODULE MOD_ISO_VARYING_STRING
