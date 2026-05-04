set(This TorsionRigidModel_test)

find_package(Sofa.Testing REQUIRED)

set(HEADER_FILES
        ../math/Lie.h
        ../core/RobotParameters.h
        ../core/ForwardKinematics.h
        ../core/InverseKinematics.h
        )

set(SOURCE_FILES
        test_lie.cpp
        test_fk.cpp
        test_ik.cpp
        test_sofa_engine.cpp)

add_executable(${This} ${SOURCE_FILES} ${HEADER_FILES})

target_link_libraries(${This}
        Sofa.Testing
        TorsionRigidModel
)

target_include_directories(${This}
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>"
    )

add_test(
    NAME ${This}
    COMMAND ${This}
)