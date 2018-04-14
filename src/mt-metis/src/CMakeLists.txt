include_directories(.)
file(GLOB mtmetis_sources *.c)
file(GLOB domlib_sources 
  ${CMAKE_SOURCE_DIR}/${DOMLIB_PATH}/*.c)
file(GLOB bowstring_sources 
  ${BOWSTRING_LIB_PATH}/*.c 
  ${BOWSTRING_LIB_PATH}/io/*.c)
file(GLOB metis_sources 
  ${METIS_LIB_PATH}/*.c
  ${METIS_LIB_PATH}/../GKlib/memory.c
  ${METIS_LIB_PATH}/../GKlib/mcore.c
  ${METIS_LIB_PATH}/../GKlib/timers.c
  ${METIS_LIB_PATH}/../GKlib/random.c
  ${METIS_LIB_PATH}/../GKlib/fs.c
  ${METIS_LIB_PATH}/../GKlib/error.c
  ${METIS_LIB_PATH}/../GKlib/io.c
  ${METIS_LIB_PATH}/../GKlib/gk_util.c
  ${METIS_LIB_PATH}/../GKlib/string.c
)
list(REMOVE_ITEM bowstring_sources 
  ${BOWSTRING_LIB_PATH}/bowstring_bin.c)
# libraries
add_library(mtmetis ${MTMETIS_LIBRARY_TYPE} 
  ${mtmetis_sources}
  ${metis_sources}
  ${bowstring_sources}
  ${domlib_sources}
) 
install(TARGETS mtmetis 
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
)



# binary
add_executable(mtmetis_bin mtmetis_bin.c)
set_target_properties(mtmetis_bin PROPERTIES OUTPUT_NAME mtmetis)
target_link_libraries(mtmetis_bin mtmetis ${REQ_LIBS})
install(TARGETS mtmetis_bin
  RUNTIME DESTINATION bin
)
