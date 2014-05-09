find_package (GSL)

include_directories ("${GSL_INCLUDE_DIRS}")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GSL_DEFINITIONS}")

include_directories("${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/abacus/src"
    "${PROJECT_SOURCE_DIR}/src/hashtab"
    )

find_library (MATH_LIBRARY m)

add_executable (sumstatsvector
	msStatsDQH.c
	sumStatsVector.c
	sumStatsVector.h
	hashtab/hashtab.c
	hashtab/hashtab.h
	whiteSpaces.c
	whiteSpaces.h
	stringUtils.c
	stringUtils.h
    )
target_link_libraries(sumstatsvector
    "${GSL_LIBRARIES}"
    "${MATH_LIBRARY}"
    )

add_executable(msDQH
	msQHarbpop.c
	msQHarbpop.h
	QHarbpop-streec.c
	QHsubs.c
	QHsubs.h
	QHarbpop-getpars.c
	rand2new.c
    )
target_link_libraries(msDQH
    "${GSL_LIBRARIES}"
    "${MATH_LIBRARY}"
    )

add_executable(dpp-msprior
	msprior.c
	msprior.h
	setup.c
	setup.h
	initvars.c
	initvars.h
	whiteSpaces.c
	whiteSpaces.h
	hashtab/hashtab.c
	hashtab/hashtab.h
	stringUtils.c
	stringUtils.h
	abacus/src/array_utils.c
	abacus/src/array_utils.h
	abacus/src/probability.c
	abacus/src/probability.h
	abacus/src/parsing.c
	abacus/src/parsing.h
	abacus/src/partition_combinatorics.c
	abacus/src/partition_combinatorics.h
	abacus/src/partition_combinatorics_random.c
	abacus/src/partition_combinatorics_random.h
    )
target_link_libraries(dpp-msprior
    "${GSL_LIBRARIES}"
    "${MATH_LIBRARY}"
    )

add_executable(msReject
    rejectinC.c
    )
target_link_libraries(msReject
    "${MATH_LIBRARY}"
    )

install(TARGETS
    sumstatsvector
    msDQH
    dpp-msprior
    msReject
    DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

install (FILES
 	acceptRej.pl
 	acceptRej_flex.pl
 	convertIM.pl
 	dataSummary.pl
 	dpp-msbayes.pl
 	msCombModels.pl
 	obsSumStats.pl
 	regress_cli.r
    DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

install (FILES
 	acceptRej_flex.r
 	acceptRej.r
 	calmod.r
 	loc2plot.r
 	make_pd2005.r
    DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/msbayes")
