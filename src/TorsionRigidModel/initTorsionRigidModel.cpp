/**
 * @file initTorsionRigidModel.cpp
 * @author Haitao Gao (haitao.gao@unsw.edu.au)
 * @brief Plugin entry point required by SOFA's PluginManager
 * @version 0.1
 * @date 2026-05-06
 *
 * @copyright Copyright (c) 2026
 */

#include <sofa/core/ObjectFactory.h>

namespace TRMCTR::engine   { extern void registerTRMForwardKinematicsEngine(sofa::core::ObjectFactory*); }
namespace TRMCTR::visual   { extern void registerTRMVisualModel(sofa::core::ObjectFactory*); }

extern "C" {

void initExternalModule() {}

const char* getModuleName()    { return "TorsionRigidModel"; }
const char* getModuleVersion() { return "0.1"; }
const char* getModuleLicense() { return "LGPL"; }
const char* getModuleDescription()
{
    return "CTR Torsion Rigid Model: Forward/Inverse Kinematics and Visual model "
           "for Concentric Tube Robots in SOFA.";
}

void registerObjects(sofa::core::ObjectFactory* factory)
{
    TRMCTR::engine::registerTRMForwardKinematicsEngine(factory);
    TRMCTR::visual::registerTRMVisualModel(factory);
}

} // extern "C"
