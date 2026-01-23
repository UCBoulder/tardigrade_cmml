if(${project_name_upper}_BUILD_AS_SHARED)
    add_library(${PROJECT_NAME} SHARED ${PROJECT_SOURCE_FILES})
    set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${PROJECT_NAME}.h" SUFFIX ".so")
    target_compile_options(${PROJECT_NAME} PUBLIC)
    target_link_libraries(
        ${PROJECT_NAME}
        PUBLIC
            "tardigrade_error_tools"
            ${ADDITIONAL_HEADER_ONLY_UPSTREAM_LIBRARIES}
            ${SHARED_OBJECT_MATERIALS_INTERFACES}
    )
    if(NOT cmake_build_type_lower STREQUAL "release")
        foreach(additional_lib ${ADDITIONAL_HEADER_ONLY_UPSTREAM_LIBRARIES})
            target_include_directories(${PROJECT_NAME} PUBLIC "${${additional_lib}_SOURCE_DIR}/${CPP_SRC_PATH}")
        endforeach(additional_lib)
    endif()
    foreach(user_interface ${SHARED_OBJECT_MATERIALS_INTERFACES})
        target_link_options(
            ${PROJECT_NAME}
            PUBLIC "SHELL:-Wl,--whole-archive $<TARGET_FILE:${user_interface}> -Wl,--no-whole-archive"
        )
    endforeach(user_interface)
else()
    add_library(${PROJECT_NAME} INTERFACE ${PROJECT_SOURCE_FILES})
    target_include_directories(
        ${PROJECT_NAME}
        INTERFACE
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/${CPP_SRC_PATH}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
endif()

install(FILES ${PROJECT_PRIVATE_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

foreach(package ${PROJECT_NAME})
    install(
        TARGETS ${package}
        EXPORT ${package}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
    get_target_property(target_sources ${package} SOURCES)
    foreach(source ${target_sources})
        set_property(GLOBAL APPEND PROPERTY CLANG_FORMAT_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
    endforeach(source)
endforeach(package)

install(FILES ${PROJECT_SOURCE_FILES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

foreach(package ${ADDITIONAL_HEADER_ONLY_LIBRARIES})
    install(
        TARGETS ${package}
        EXPORT ${package}_Targets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
    install(FILES ${package}.h ${package}.cpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    get_target_property(target_sources ${package} SOURCES)
    foreach(source ${target_sources})
        set_property(GLOBAL APPEND PROPERTY CLANG_FORMAT_SOURCE_FILES "${PROJECT_SOURCE_DIR}/${source}")
    endforeach(source)
endforeach(package)

get_target_property(target_sources ${PROJECT_NAME} SOURCES)
foreach(source ${target_sources})
    set_property(GLOBAL APPEND PROPERTY CLANG_FORMAT_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
endforeach(source)
