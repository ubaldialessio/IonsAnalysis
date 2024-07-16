cmake_minimum_required(VERSION 3.13)
project(IonsAnalysis)
set(CMAKE_CXX_STANDARD 14)

find_package(NAIA REQUIRED)
find_package(NSL REQUIRED)
find_package(Python3 COMPONENTS Interpreter REQUIRED)


# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS RooFit RooFitCore RooStats)
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

# Add TrTrackEffSel.cpp to the project sources
set(SOURCES_selectEv
    src/selectEv.cpp
    include/Efficiencies/Tracker.cpp
    include/Efficiencies/Trigger.cpp
    include/Efficiencies/InnerTracker.cpp
    include/Efficiencies/TofStandalone.cpp
    include/Efficiencies/TRD.cpp
    include/Efficiencies/RTI.cpp
    include/Efficiencies/MySel.cpp
)

set(SOURCES_fragment
    src/fragmentation.cpp
    include/Efficiencies/Tracker.cpp
    include/Efficiencies/Trigger.cpp
    include/Efficiencies/InnerTracker.cpp
    include/Efficiencies/TofStandalone.cpp
    include/Efficiencies/TRD.cpp
    include/Efficiencies/RTI.cpp
    include/Efficiencies/MySel.cpp
    include/Templates/BuildTemplatesSel.cpp
)

add_executable(checkVar src/checkVar.cpp)
target_link_libraries(checkVar NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(checkVar PRIVATE include/)
set_target_properties(checkVar PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

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

add_executable(BuildFlux src/BuildFlux.cpp)
target_link_libraries(BuildFlux NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(BuildFlux PRIVATE include/)
set_target_properties(BuildFlux PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(fitTemplates src/fitTemplates.cpp)
target_link_libraries(fitTemplates NAIA::NAIAChain NSL::NSLSelections RooFitCore RooFitMore)
target_include_directories(fitTemplates PRIVATE include/)
set_target_properties(fitTemplates PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareEff src/CompareEff.cpp)
target_link_libraries(CompareEff NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareEff PRIVATE include/)
set_target_properties(CompareEff PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(printEverything src/printEverything.cpp)
target_link_libraries(printEverything NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(printEverything PRIVATE include/)
set_target_properties(printEverything PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareFlux src/CompareFlux.cpp)
target_link_libraries(CompareFlux NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareFlux PRIVATE include/)
set_target_properties(CompareFlux PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

