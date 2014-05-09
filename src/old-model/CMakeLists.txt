if (INCLUDE_OLD_MODEL)

    find_package (GSL)
    
    include_directories ("${GSL_INCLUDE_DIRS}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GSL_DEFINITIONS}")
    
    include_directories("${PROJECT_SOURCE_DIR}/src"
        "${PROJECT_SOURCE_DIR}/abacus/src"
        "${PROJECT_SOURCE_DIR}/src/hashtab"
        "${PROJECT_SOURCE_DIR}/src/old-model"
        )
    
    find_library (MATH_LIBRARY m)
    
    add_executable (msprior
    	msprior.c
    	msprior.h
    	setup.c
    	setup.h
    	../initvars.c
    	../initvars.h
    	../whiteSpaces.c
    	../whiteSpaces.h
    	../hashtab/hashtab.c
    	../hashtab/hashtab.h
    	../stringUtils.c
    	../stringUtils.h
    	../abacus/src/array_utils.c
    	../abacus/src/array_utils.h
    	../abacus/src/probability.c
    	../abacus/src/probability.h
    	../abacus/src/parsing.c
    	../abacus/src/parsing.h
    	../abacus/src/partition_combinatorics.c
    	../abacus/src/partition_combinatorics.h
    	../abacus/src/partition_combinatorics_random.c
    	../abacus/src/partition_combinatorics_random.h
        )
    target_link_libraries(msprior
        "${GSL_LIBRARIES}"
        "${MATH_LIBRARY}"
        )
    
    install (TARGETS msprior DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
    install (FILES msbayes.pl DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

endif (INCLUDE_OLD_MODEL)
