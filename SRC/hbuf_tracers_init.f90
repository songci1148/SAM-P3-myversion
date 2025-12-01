! Wrapper call to avoid circular dependency between hbuffer and tracers modules

subroutine hbuf_tracers_init(namelist,deflist,unitlist,status,average_type,count,trcount)
   use tracers
   implicit none
   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*),average_type(*),count,trcount

   call tracers_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)

end
