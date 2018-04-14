include_directories(.)
file(GLOB bowstring_sources *.c io/*.c)
file(GLOB domlib_sources ${CMAKE_SOURCE_DIR}/${DOMLIB_PATH}/*.c)

# library                                         
add_library(bowstring ${BOWSTRING_LIBRARY_TYPE} 
  ${bowstring_sources}
  ${domlib_sources}
) 
install(TARGETS bowstring 
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# binary
add_executable(bowstring_bin bowstring_bin.c)
set_target_properties(bowstring_bin PROPERTIES OUTPUT_NAME bowstring)
target_link_libraries(bowstring_bin bowstring ${REQ_LIBS})
install(TARGETS bowstring_bin
  RUNTIME DESTINATION bin
)


if (NOT DEFINED NO_MPI OR NO_MPI EQUAL 0)
  add_subdirectory("mpi")
endif()
