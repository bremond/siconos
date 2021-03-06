# ============================================================
# Prepare  config to generate siconos documentation :
#
# - automatic doxygen documentation, from source code
# - sphinx documentation for manuals
# 
# ============================================================
include(FindPythonModule)
include(FindSphinxModule)

set(DOC_ROOT_DIR ${CMAKE_BINARY_DIR}/Docs/build/html/ CACHE INTERNAL "Root path to documentation - Corresponds to 'htdocs' on gforge page.")

#  ====== doxygen documentation setup ======
# This setup is required if one at least of
# the following is needed:
# - doxygen warnings (see WITH_DOXYGEN_WARNINGS)
# - python docstrings from doxygen (xml) (see WITH_DOXY2SWIG)
# - doxygen documentation generated from source files (see WITH_DOCUMENTATION)
if(USE_DOXYGEN)
  find_package(Doxygen REQUIRED)

  # Output path for doxygen documentation
  if(NOT DOXYGEN_OUTPUT)
    set(DOXYGEN_OUTPUT ${DOC_ROOT_DIR}/doxygen CACHE INTERNAL "Path to doxygen outputs.")
  endif()
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/tags)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/html)
  
  file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT}/man)
  
  # --- Config file and its options. --
  # --> one for doxygen warnings (from doxy_warnings.config.in)
  # --> one for documentation and python docstrings (from doxy.config.in)
  
  # == set configuration file for doxygen warnings ==
  #  - Results in binary_dir/doxygen_warnings
  #  - input config from config/doxy_warnings.config
  #  - only xml output.
  if(WITH_DOXYGEN_WARNINGS)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doxygen_warnings)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doxygen_warnings/man)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doxygen_warnings/xml)
    set(DOXY_WARNINGS_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/config/doxy_warnings.config.in CACHE INTERNAL "config file for doxygen warnings output.")
    # All the configure process will be done for each component
    # thanks to 'include(doxygen_warnings)'. See LibraryProjectSetup.
  endif()

  # == set configuration file for doxygen doc from sources ==
  #  - Results in binary_dir/Docs/
  #  - input config from config/doxy.config
  #  - only html output (add option for latex?)
  # 'inputs' are updated by each component, during call to update_xml_doxy_config_file macro.
  # Doc will be built when 'make doxygen' is called.
  if(WITH_DOCUMENTATION)
    # config file name
    set(DOXY_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/config/doxy.config" CACHE INTERNAL "Doxygen configuration file for doc from sources.")

    # The generation of the config and the creation of the target will be done later,
    # after the update of inputs by each component, with macro 'finalize_doxygen'
    
  endif()
    
  # == set configuration file for doxy2swig (use doxygen ouput to create docstrings) ==
  #  - Results in binary_dir/Docs/
  #  - input config from config/<component_name>doxy.config
  #  - only xml output
  # 'inputs' are updated by each component, during call to update_xml_doxy_config_file macro.
  # The xml doc is built for each component, during call to build_doc_xml macro.
  # Here, we just need to create a target to create xml files for a given component
  # for example --> make numerics_xml to create xml files for numerics.
  if(WITH_DOXY2SWIG)
    foreach(comp ${COMPONENTS})
      set(DOXY_CONFIG_XML ${CMAKE_CURRENT_BINARY_DIR}/config/${comp}doxy.config.xml)
      add_custom_target(${comp}_xml ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG_XML} WORKING_DIRECTORY ${DOC_ROOT_DIR}
	COMMENT "Build xml documentation files.")
    endforeach()
  endif()

endif()

  
#  ====== Create a Target to generate sphinx documentation ======
if(USE_SPHINX)
  # find_python_module(breathe REQUIRED)
  find_python_module(sphinx REQUIRED)
  find_sphinx_module(sphinxcontrib bibtex REQUIRED)
  find_sphinx_module(sphinxcontrib doxylink REQUIRED)
  
  # --- conf.py describes sphinx setup ---
  configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/conf.py.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/conf.py" @ONLY)

  configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/index.rst.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/index.rst" @ONLY)

  configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/reference/siconos_api_reference.rst.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/reference/siconos_api_reference.rst" @ONLY)
  
  # --- Build a target for sphinx ---
  # --> use make html or make latex ... to build documentation.
  set(SPHINX_EXECUTABLE sphinx-build)

  # html doc
  set(SPHINX_PARAMETERS
    -b html # html output
    -d build/doctrees # path to doctree files
    ${CMAKE_CURRENT_SOURCE_DIR}/sphinx  # path to rst source
    )

  add_custom_target(html
    COMMAND ${SPHINX_EXECUTABLE} ${SPHINX_PARAMETERS}
    ${CMAKE_CURRENT_BINARY_DIR}/build/html
    )
  
  # latex doc
  set(SPHINX_LATEX_PARAMETERS
    -b latex # latex output
    -d build/doctrees # path to doctree files
    ${CMAKE_CURRENT_SOURCE_DIR}/sphinx  # path to rst source
    )
  
  add_custom_target(latex
    COMMAND ${SPHINX_EXECUTABLE} ${SPHINX_LATEX_PARAMETERS}
    ${CMAKE_CURRENT_BINARY_DIR}/build/latex
    )

  # python api auto doc
  if(WITH_PYTHON_WRAPPER)

    set(SPHINX_APIDOC sphinx-apidoc)
    set(GENERATED_DOC ${CMAKE_CURRENT_BINARY_DIR}/sphinx/reference)
    set(APIDOC_PARAMETERS -e -P -o ${GENERATED_DOC})
    list(GET COMPONENTS -1 _last_comp)
    add_custom_target(apidoc
      COMMAND ${SPHINX_APIDOC} ${APIDOC_PARAMETERS}
      ${CMAKE_BINARY_DIR}/wrap/${SICONOS_PYTHON_PACKAGE}
      DEPENDS ${_last_comp}
      )
  endif()

  add_custom_target(doc DEPENDS html)

endif()
