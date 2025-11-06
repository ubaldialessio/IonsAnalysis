cmake_minimum_required(VERSION 3.13)
project(IonsAnalysis)
set(CMAKE_CXX_STANDARD 17)

find_package(NAIA REQUIRED)
find_package(NSL REQUIRED)
find_package(Python3 COMPONENTS Interpreter REQUIRED)


# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS RooFit RooFitCore RooStats)
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

# Add TrTrackEffSel.cpp to the project sources
set(SOURCES_selectEv
    src/Selection/selectEv.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_SelectPassedTree
    src/Selection/SelectPassedTree.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_SelectCounts
    src/Selection/SelectCounts.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_RawAcceptance
    src/Selection/RawAcceptance.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_AboveL1
    src/Fragmentation/AboveL1.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_L1Eff
    src/Efficiencies/L1Eff.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_L1EffUnb
    src/Efficiencies/L1EffUnb.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_TrackEff
    src/Efficiencies/TrackEff.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_TrackChEff
    src/Efficiencies/TrackChEff.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/ChargeAnalysisManager.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_TofEff
    src/Efficiencies/TofEff.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_TriggEff
    src/Efficiencies/TriggEff.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/ChargeAnalysisManager.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_DaqEff
    src/Efficiencies/DaqEff.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_fragment
    src/Fragmentation/fragmentation.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_l1temp
    src/Fragmentation/BuildL1Template.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_RigidityMap
    src/RigidityMap.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_SecondTrackVsFirst
    src/ControlPlots/SecondTrackVsFirst.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_ClustersNumbers
    src/ControlPlots/ClustersNumbers.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_VariableController
    src/VariableController.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_CreateTemplatesBelowL1
    src/Fragmentation/CreateTemplatesBelowL1.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_FitTemplates
    src/Fragmentation/fitTemplates.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_ShiftTemplates
    src/Fragmentation/ShiftTemplates.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)


set(SOURCES_EvaluateFragments
    src/Fragmentation/EvaluateFragments.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_BuildFlux
    src/BuildFlux.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_BuildFluxForJack
    src/BuildFluxForJack.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_Livetime
    src/Selection/Livetime.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_InterpolateCorrection
    src/Efficiencies/InterpolateCorrections.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_BuildAcceptance
    src/Efficiencies/BuildAcceptance.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_UnfoldAcceptance
    src/Unfolding/UnfoldAcceptance.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_FitPurity
    src/Fragmentation/FitPurity.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_TestFluxFit
    src/Plotting/TestFluxFit.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_PurityMC
    src/Fragmentation/PurityMc.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_DataSelector
    src/Selection/DataSelector.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/ChargeAnalysisManager.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_DistributionsSelector
    src/Selection/DistributionsSelector.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/ChargeAnalysisManager.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_McSelector
    src/Selection/McSelector.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/ChargeAnalysisManager.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)

set(SOURCES_TemplatesSelector
    src/Selection/TemplatesSelector.cpp
    src/SelectionEfficiency/Tracker.cpp
    src/SelectionEfficiency/Trigger.cpp
    src/SelectionEfficiency/InnerTracker.cpp
    src/SelectionEfficiency/TofStandalone.cpp
    src/SelectionEfficiency/TRD.cpp
    src/SelectionEfficiency/RTI.cpp
    src/SelectionEfficiency/MySel.cpp
    src/ChargeAnalysisManager.cpp
    src/SelectionEfficiency/BuildTemplatesSel.cpp
)


add_executable(VariableController ${SOURCES_VariableController})
target_link_libraries(VariableController NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(VariableController PRIVATE include/)
set_target_properties(VariableController PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(selectEv ${SOURCES_selectEv})
target_link_libraries(selectEv NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(selectEv PRIVATE include/)
set_target_properties(selectEv PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(hplot src/Plotting/hplot.cpp)
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

add_executable(CreateTemplatesBelowL1 ${SOURCES_CreateTemplatesBelowL1})
target_link_libraries(CreateTemplatesBelowL1 NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CreateTemplatesBelowL1 PRIVATE include/)
set_target_properties(CreateTemplatesBelowL1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(BuildFlux ${SOURCES_BuildFlux})
target_link_libraries(BuildFlux NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(BuildFlux PRIVATE include/)
set_target_properties(BuildFlux PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(BuildFluxForJack ${SOURCES_BuildFluxForJack})
target_link_libraries(BuildFluxForJack NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(BuildFluxForJack PRIVATE include/)
set_target_properties(BuildFluxForJack PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(fitTemplates ${SOURCES_FitTemplates})
target_link_libraries(fitTemplates NAIA::NAIAChain NSL::NSLSelections RooFitCore RooFit)
target_include_directories(fitTemplates PRIVATE include/)
set_target_properties(fitTemplates PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareEff src/Plotting/CompareEff.cpp)
target_link_libraries(CompareEff NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareEff PRIVATE include/)
set_target_properties(CompareEff PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareRatios src/Plotting/CompareRatios.cpp)
target_link_libraries(CompareRatios NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareRatios PRIVATE include/)
set_target_properties(CompareRatios PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareRatiosAll src/Plotting/CompareRatiosAll.cpp)
target_link_libraries(CompareRatiosAll NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareRatiosAll PRIVATE include/)
set_target_properties(CompareRatiosAll PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(printEverything src/Plotting/printEverything.cpp)
target_link_libraries(printEverything NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(printEverything PRIVATE include/)
set_target_properties(printEverything PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareFlux src/Plotting/CompareFlux.cpp)
target_link_libraries(CompareFlux NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareFlux PRIVATE include/)
set_target_properties(CompareFlux PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(BuildL1Template ${SOURCES_l1temp})
target_link_libraries(BuildL1Template NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(BuildL1Template PRIVATE include/)
set_target_properties(BuildL1Template PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(RigidityMap ${SOURCES_RigidityMap})
target_link_libraries(RigidityMap NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(RigidityMap PRIVATE include/)
set_target_properties(RigidityMap PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(L1Eff ${SOURCES_L1Eff})
target_link_libraries(L1Eff NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(L1Eff PRIVATE include/)
set_target_properties(L1Eff PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(L1EffUnb ${SOURCES_L1EffUnb})
target_link_libraries(L1EffUnb NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(L1EffUnb PRIVATE include/)
set_target_properties(L1EffUnb PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(TofEff ${SOURCES_TofEff})
target_link_libraries(TofEff NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(TofEff PRIVATE include/)
set_target_properties(TofEff PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(TrackEff ${SOURCES_TrackEff})
target_link_libraries(TrackEff NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(TrackEff PRIVATE include/)
set_target_properties(TrackEff PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(TrackChEff ${SOURCES_TrackChEff})
target_link_libraries(TrackChEff NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(TrackChEff PRIVATE include/)
set_target_properties(TrackChEff PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(TriggEff ${SOURCES_TriggEff})
target_link_libraries(TriggEff NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(TriggEff PRIVATE include/)
set_target_properties(TriggEff PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(DaqEff ${SOURCES_DaqEff})
target_link_libraries(DaqEff NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(DaqEff PRIVATE include/)
set_target_properties(DaqEff PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareWithPublished src/Plotting/CompareWithPublished.cpp)
target_link_libraries(CompareWithPublished NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareWithPublished PRIVATE include/)
set_target_properties(CompareWithPublished PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(InterpolateCorrections ${SOURCES_InterpolateCorrection})
target_link_libraries(InterpolateCorrections NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(InterpolateCorrections PRIVATE include/)
set_target_properties(InterpolateCorrections PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(BuildAcceptance ${SOURCES_BuildAcceptance})
target_link_libraries(BuildAcceptance NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(BuildAcceptance PRIVATE include/)
set_target_properties(BuildAcceptance PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareCorrections src/Plotting/CompareCorrections.cpp)
target_link_libraries(CompareCorrections NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareCorrections PRIVATE include/)
set_target_properties(CompareCorrections PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(RawAcceptance ${SOURCES_RawAcceptance})
target_link_libraries(RawAcceptance  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(RawAcceptance PRIVATE include/)
set_target_properties(RawAcceptance  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(AboveL1 ${SOURCES_AboveL1})
target_link_libraries(AboveL1  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(AboveL1 PRIVATE include/)
set_target_properties(AboveL1  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(Livetime ${SOURCES_Livetime})
target_link_libraries(Livetime  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(Livetime PRIVATE include/)
set_target_properties(Livetime  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(SelectCounts ${SOURCES_SelectCounts})
target_link_libraries(SelectCounts  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(SelectCounts PRIVATE include/)
set_target_properties(SelectCounts PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(SelectPassedTree ${SOURCES_SelectPassedTree})
target_link_libraries(SelectPassedTree  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(SelectPassedTree PRIVATE include/)
set_target_properties(SelectPassedTree  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(UnfoldAcceptance ${SOURCES_UnfoldAcceptance})
target_link_libraries(UnfoldAcceptance NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(UnfoldAcceptance PRIVATE include/)
set_target_properties(UnfoldAcceptance PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(BuildFileToPrint src/Plotting/BuildFileToPrint.cpp)
target_link_libraries(BuildFileToPrint NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(BuildFileToPrint PRIVATE include/)
set_target_properties(BuildFileToPrint PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(TotalCorrection  src/Plotting/TotalCorrection.cpp)
target_link_libraries(TotalCorrection  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(TotalCorrection  PRIVATE include/)
set_target_properties(TotalCorrection  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareUnfoldIteration  src/Plotting/CompareUnfoldIteration.cpp)
target_link_libraries(CompareUnfoldIteration  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareUnfoldIteration  PRIVATE include/)
set_target_properties(CompareUnfoldIteration  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(AcceptanceMC  src/Plotting/AcceptanceMC.cpp)
target_link_libraries(AcceptanceMC  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(AcceptanceMC  PRIVATE include/)
set_target_properties(AcceptanceMC  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(EvaluateFragments  ${SOURCES_EvaluateFragments})
target_link_libraries(EvaluateFragments  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(EvaluateFragments  PRIVATE include/)
set_target_properties(EvaluateFragments  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(InnerChargeDistribution  src/ControlPlots/InnerChargeDistribution.cpp)
target_link_libraries(InnerChargeDistribution  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(InnerChargeDistribution  PRIVATE include/)
set_target_properties(InnerChargeDistribution  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(L1Distribution  src/ControlPlots/L1Distribution.cpp)
target_link_libraries(L1Distribution  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(L1Distribution  PRIVATE include/)
set_target_properties(L1Distribution  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(L1UnbDistribution  src/ControlPlots/L1UnbDistribution.cpp)
target_link_libraries(L1UnbDistribution  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(L1UnbDistribution  PRIVATE include/)
set_target_properties(L1UnbDistribution  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(L9UnbDistribution  src/ControlPlots/L9UnbDistribution.cpp)
target_link_libraries(L9UnbDistribution  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(L9UnbDistribution  PRIVATE include/)
set_target_properties(L9UnbDistribution  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(TofDistribution  src/ControlPlots/TofDistribution.cpp)
target_link_libraries(TofDistribution  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(TofDistribution  PRIVATE include/)
set_target_properties(TofDistribution  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareGlobalMc src/Plotting/CompareGlobalMc.cpp)
target_link_libraries(CompareGlobalMc NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareGlobalMc PRIVATE include/)
set_target_properties(CompareGlobalMc PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareRates src/Plotting/CompareRates.cpp)
target_link_libraries(CompareRates NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareRates PRIVATE include/)
set_target_properties(CompareRates PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(PlotUnfoldingFactor src/Plotting/PlotUnfoldingFactor.cpp)
target_link_libraries(PlotUnfoldingFactor NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(PlotUnfoldingFactor PRIVATE include/)
set_target_properties(PlotUnfoldingFactor PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareAcceptanceMC src/Plotting/CompareAcceptanceMC.cpp)
target_link_libraries(CompareAcceptanceMC NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareAcceptanceMC PRIVATE include/)
set_target_properties(CompareAcceptanceMC PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareEffectiveAcc src/Plotting/CompareEffectiveAcc.cpp)
target_link_libraries(CompareEffectiveAcc NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareEffectiveAcc PRIVATE include/)
set_target_properties(CompareEffectiveAcc PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(SecondTrackVsFirst ${SOURCES_SecondTrackVsFirst})
target_link_libraries(SecondTrackVsFirst  NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(SecondTrackVsFirst  PRIVATE include/)
set_target_properties(SecondTrackVsFirst  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(ShiftTemplates ${SOURCES_ShiftTemplates})
target_link_libraries(ShiftTemplates NAIA::NAIAChain NSL::NSLSelections RooFitCore RooFit)
target_include_directories(ShiftTemplates PRIVATE include/)
set_target_properties(ShiftTemplates PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(CompareL1ChargeSecTrack src/Plotting/CompareL1ChargeSecTrack.cpp)
target_link_libraries(CompareL1ChargeSecTrack NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(CompareL1ChargeSecTrack PRIVATE include/)
set_target_properties(CompareL1ChargeSecTrack PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(ClustersNumbers ${SOURCES_ClustersNumbers})
target_link_libraries(ClustersNumbers NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(ClustersNumbers PRIVATE include/)
set_target_properties(ClustersNumbers PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)


add_executable(FitPurity ${SOURCES_FitPurity})
target_link_libraries(FitPurity NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(FitPurity PRIVATE include/)
set_target_properties(FitPurity PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(TestFluxFit ${SOURCES_TestFluxFit})
target_link_libraries(TestFluxFit NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(TestFluxFit PRIVATE include/)
set_target_properties(TestFluxFit PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(PurityMC ${SOURCES_PurityMC})
target_link_libraries(PurityMC NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(PurityMC PRIVATE include/)
set_target_properties(PurityMC PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(DataSelector ${SOURCES_DataSelector})
target_link_libraries(DataSelector NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(DataSelector PRIVATE include/)
set_target_properties(DataSelector PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(DistributionsSelector ${SOURCES_DistributionsSelector})
target_link_libraries(DistributionsSelector NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(DistributionsSelector PRIVATE include/)
set_target_properties(DistributionsSelector PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(McSelector ${SOURCES_McSelector})
target_link_libraries(McSelector NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(McSelector PRIVATE include/)
set_target_properties(McSelector PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)

add_executable(TemplatesSelector ${SOURCES_TemplatesSelector})
target_link_libraries(TemplatesSelector NAIA::NAIAChain NSL::NSLSelections)
target_include_directories(TemplatesSelector PRIVATE include/)
set_target_properties(TemplatesSelector PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./../bin/)