cmake_minimum_required(VERSION 3.21)

project(
  PENKNIFE
  VERSION 0.0.1
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# ##############################################################################
# Find dependencies
# ##############################################################################

find_package(NESO REQUIRED)
find_package(NESO-Particles REQUIRED)
find_package(NESO-RNG-Toolkit REQUIRED)
find_package(VANTAGE-Reactions REQUIRED)

# Identify source files

set(PENKNIFEClosureTypes Braginskii Zhdanov)

# The combined volume flux uses the same explicit registration as FluxVolDiff.
# Retain the generator below for future additional volume-flux implementations.
set(PENKNIFEVolFluxes)
foreach(relPath IN LISTS PENKNIFEVolFluxes)
  set(dir ${CMAKE_CURRENT_SOURCE_DIR}/${relPath})
  get_filename_component(name "${relPath}" NAME)
  set(NAME_SOURCES ${name}_SOURCES)
  set(OPERATOR ${name}) # Operator parameter for boiler plate code

  # These operator have an additional argument for closure and SolverTypes which
  # requires more declarations.
  set(CONFIG_FILE ${CMAKE_SOURCE_DIR}/VolumeFluxFactoryDec.cpp.in)

  foreach(ClosureType IN LISTS PENKNIFEClosureTypes)
    # Loop through each possible execution space.
    foreach(ExecSpace IN LISTS ExecSpaces)

      # Set file extension.
      set(HEADER_EXT hpp)
      if("${ExecSpace}" STREQUAL "Device" AND NEKTAR_ENABLE_DEVICE STREQUAL
                                              "CUDA")
        set(SOURCE_EXT cu)
      elseif("${ExecSpace}" STREQUAL "Device" AND NEKTAR_ENABLE_DEVICE STREQUAL
                                                  "HIP")
        set(SOURCE_EXT hip)
      else()
        set(SOURCE_EXT cpp)
      endif()

      # Loop through each possible data type.
      foreach(TData IN LISTS DataTypes)
        # Set up the name used for the .cpp declaration file.
        set(FactoryDeclName
            src/${name}${SolverType}${ClosureType}${ExecSpace}${TData}.${SOURCE_EXT}
        )

        # Get the source and header files for this operator, execution space,
        # and implementation.
        set(CURRENT_HEADERS ${name}_${ExecSpace}_HEADERS)
        set(CURRENT_SOURCES ${name}_${ExecSpace}_SOURCES)

        # Get the headers matching the name and execution space
        file(
          GLOB ${CURRENT_HEADERS}
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${dir}/${name}${SolverType}${ExecSpace}*.hpp)

        file(RELATIVE_PATH include_file ${CMAKE_CURRENT_SOURCE_DIR}
             "${dir}/${name}OpImpl.hpp")
        set(IMPL_HEADER "#include \"${include_file}\"")

        file(RELATIVE_PATH include_file ${CMAKE_CURRENT_SOURCE_DIR}
             "${dir}/SolverKernels.hpp")
        set(KERNEL_HEADER "#include \"${include_file}\"")

        # Anchored at the solver root rather than relative to the operator, so
        # an operator in a subdirectory still finds the equation of state
        # headers.
        file(RELATIVE_PATH include_file ${CMAKE_CURRENT_SOURCE_DIR}
             "${CMAKE_CURRENT_SOURCE_DIR}/Closures/${ClosureType}.hpp")
        set(CLOSURE_HEADER "#include \"${include_file}\"")
        set(CLOSURE_TYPE ${ClosureType})

        # Create the cpp file from the implementation file.
        configure_file(${CONFIG_FILE} ${FactoryDeclName})

        # Add the respective cpp file.
        set(${CURRENT_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/${FactoryDeclName})
        set(${NAME_SOURCES} ${${NAME_SOURCES}} ${${CURRENT_SOURCES}})

      endforeach()
    endforeach()
  endforeach()

  # Add this operator's implemenations to the global operator source and header
  # files.
  set(PENKNIFE_OPERATORS_SOURCES ${PENKNIFE_OPERATORS_SOURCES}
                                 ${${NAME_SOURCES}})

  # Make list visible to subdirectories
  set(VOLUME_FACTORY_SOURCES
      "${PENKNIFE_OPERATORS_SOURCES}"
      CACHE INTERNAL "Generated Volume flux factory sources")
endforeach()

# ##############################################################################

# Trace fluxes use explicit registration, matching the volume operators.

# Register PENKNIFE operators without Nektar++'s internal ADD_OPERATOR macro.
list(
  APPEND
  PENKNIFE_OPERATORS_SOURCES
  src/BndCond/BndCondBohm/BndCondBohmFactory.cpp
  src/MagneticField/MagneticFieldFactory.cpp
  src/Pointwise/PlasmaPointwiseFactory.cpp
  src/FluxVolDiff/FluxVolDiffFactory.cpp
  src/FluxVolAdv/FluxVolAdvFactory.cpp
  src/FluxVolAdvDiff/FluxVolAdvDiffFactory.cpp
  src/FluxTraceAdvDiv/FluxTraceAdvDiffFactory.cpp)

set(SRC_FILES
    src/EquationSystems/PlasmaSystem.cpp
    src/EquationSystems/NullSystem.cpp
    src/EquationSystems/AnisotropicDiffusion.cpp
    src/EquationSystems/ReducedBraginskii.cpp
    src/EquationSystems/ElectrostaticTurbulence.cpp
    src/ParticleSystems/ParticleSystem.cpp
    src/ParticleSystems/ReactionSystem.cpp
    src/ParticleSystems/Reactions.cpp
    ${PENKNIFE_OPERATORS_SOURCES})

# ##############################################################################
# Build object library; used by solver exec and tests
# ##############################################################################
set(OBJLIB_NAME ObjLib)
set(SOLVER_LIBS
    ${SOLVER_LIBS} ${OBJLIB_NAME}
    CACHE INTERNAL "")
add_library(${OBJLIB_NAME} OBJECT ${SRC_FILES})
target_compile_options(${OBJLIB_NAME} PRIVATE ${BUILD_TYPE_COMPILE_FLAGS})
target_link_libraries(
  ${OBJLIB_NAME}
  PRIVATE NESO-Particles::NESO-Particles Nektar++::nektar++ NESO::nesolib
          NESO-RNG-Toolkit::NESO-RNG-Toolkit
          VANTAGE-Reactions::VANTAGE-Reactions)

target_include_directories(${OBJLIB_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR}
                                                 ${CMAKE_CURRENT_LIST_DIR}/src)

add_sycl_to_target(TARGET ${OBJLIB_NAME} SOURCES ${SRC_FILES})

# ##############################################################################
# Build solver exec
# ##############################################################################
set(EXEC_NAME penknife)

add_executable(${EXEC_NAME} src/PENKNIFE.cpp $<TARGET_OBJECTS:${OBJLIB_NAME}>)

# Compile options
target_compile_options(${EXEC_NAME} PRIVATE ${BUILD_TYPE_COMPILE_FLAGS})

target_include_directories(${EXEC_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR})

# Linker options, target libs
target_link_options(${EXEC_NAME} PRIVATE ${BUILD_TYPE_LINK_FLAGS})
target_link_libraries(
  ${EXEC_NAME}
  PRIVATE NESO-Particles::NESO-Particles Nektar++::nektar++ NESO::nesolib
          NESO-RNG-Toolkit::NESO-RNG-Toolkit
          VANTAGE-Reactions::VANTAGE-Reactions)

add_sycl_to_target(TARGET ${EXEC_NAME} SOURCES src/PENKNIFE.cpp)

# ##############################################################################
# Install exec
# ##############################################################################
# Set default install path if not provided on the command-line
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX
      ${CMAKE_BINARY_DIR}/dist
      CACHE PATH "" FORCE)
endif()

# Install location
install(TARGETS ${EXEC_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# Exec rpath
set_target_properties(${EXEC_NAME} PROPERTIES INSTALL_RPATH ${INSTALL_RPATH}
                                              ${NEKTAR++_LIBRARY_DIRS})

option(PENKNIFE_BUILD_TESTS "Build PENKNIFE equation-system checks" OFF)
if(PENKNIFE_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests/EquationSystems)
endif()
