include_directories(.)
include_directories(..)

# binaries
add_executable(pushrelabel pushrelabel.c)
set_target_properties(pushrelabel PROPERTIES OUTPUT_NAME pushrelabel)
target_link_libraries(pushrelabel ${MPI_LIBRARIES} bowstring z m)
install(TARGETS pushrelabel RUNTIME DESTINATION bin)


