cmake_minimum_required(VERSION 3.13)
project(checkVar)
set(CMAKE_CXX_STANDARD 17)
find_package(NAIA REQUIRED)
find_package(NSL REQUIRED)

add_executable(checkVar src/checkVar.cpp)
target_link_libraries(checkVar NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(checkVar PRIVATE include/)
set_target_properties(checkVar PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

# Add TrTrackEffSel.cpp to the project sources
set(SOURCES_selectEv
    src/selectEv.cpp
    include/Efficiencies/Tracker.cpp  # Add this line
    include/Efficiencies/Trigger.cpp
    include/Efficiencies/InnerTracker.cpp
    include/Efficiencies/TofStandalone.cpp
    include/Efficiencies/TRD.cpp
    include/Efficiencies/RTI.cpp
    include/Efficiencies/MySel.cpp
)

set(SOURCES_fragment
    src/fragmentation.cpp
    include/Efficiencies/Tracker.cpp  # Add this line
    include/Efficiencies/Trigger.cpp
    include/Efficiencies/InnerTracker.cpp
    include/Efficiencies/TofStandalone.cpp
    include/Efficiencies/TRD.cpp
    include/Efficiencies/RTI.cpp
    include/Efficiencies/MySel.cpp
    include/Templates/BuildTemplatesSel.cpp
)

add_executable(selectEv ${SOURCES_selectEv})
target_link_libraries(selectEv NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(selectEv PRIVATE include/)
set_target_properties(selectEv PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(hplot src/hplot.cpp)
target_link_libraries(hplot NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(hplot PRIVATE include/)
set_target_properties(hplot PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(eff src/eff.cpp)
target_link_libraries(eff NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(eff PRIVATE include/)
set_target_properties(eff PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(fragmentation ${SOURCES_fragment})
target_link_libraries(fragmentation NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(fragmentation PRIVATE include/)
set_target_properties(fragmentation PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CreateTemplatesBelowL1 src/CreateTemplatesBelowL1.cpp)
target_link_libraries(CreateTemplatesBelowL1 NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CreateTemplatesBelowL1 PRIVATE include/)
set_target_properties(CreateTemplatesBelowL1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)