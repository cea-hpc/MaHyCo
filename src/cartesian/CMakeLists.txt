add_library(cartesian CartCartesianMeshT.h
        CartFaceDirectionMng.h
        interface/CartesianDefines.h
        CartCellDirectionMng.h
        CartItemEnumeratorT.h
        interface/CartesianNamespace.h
        CartConnectivityCellFaceNode.h
        CartItemGroup.h
        interface/CellDirectionMng.h
        CartConnectivityCellNode.h
        CartLocalIdNumberingT.h
        interface/FaceDirectionMng.h
        CartConnectivityNodeCell.h
        CartNeighCells.h
        interface/ICartesianMesh.h
        CartesianConnectivity.cc
        CartNodeDirectionMng.h
        interface/NodeDirectionMng.h
        CartesianConnectivity.h
        CartStencilDirItemT.h
        Interval3T.h
        CartesianFaceId.h
        CartTypes.h
        NodeDirectionMng.cc
        CartesianGlobal.h
        CellDirectionMng.cc
        NodeDirectionMng.h
        CartesianGridT.h
        CellDirectionMng.h
        NumberingConverterT.h
        CartesianItemSorter.cc
        CMakeLists.txt
        UnstructCartesianMeshT.h
        CartesianItemSorter.h
        FaceDirectionMng.cc
        UnstructConnectivityCellFaceNode.h
        CartesianMesh.cc
        FaceDirectionMng.h
        UnstructConnectivityCellNode.h
        CartesianMeshProperties.h
        FactCartDirectionMng.h
        UnstructNeighCells.h
        CartesianMeshT.h
        ICartesianMesh.h
        CartesianNumberingT.h
        interface/CartesianConnectivity.h)
target_include_directories(cartesian PRIVATE ..)
target_link_libraries(cartesian PUBLIC Arcane::arcane_core)
install(TARGETS cartesian EXPORT ${MAHYCO_EXPORT_TARGET})

arcane_accelerator_enable()
arcane_accelerator_add_source_files(CartesianConnectivity.cc)
arcane_accelerator_add_source_files(CartesianMesh.cc)
arcane_accelerator_add_source_files(CellDirectionMng.cc)
arcane_accelerator_add_source_files(FaceDirectionMng.cc)
arcane_accelerator_add_source_files(NodeDirectionMng.cc)
arcane_accelerator_add_source_files(CartesianItemSorter.cc)
arcane_accelerator_add_to_target(cartesian)
