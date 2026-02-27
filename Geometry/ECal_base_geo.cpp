#include "XML/Utilities.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"

using namespace dd4hep;

static Ref_t create_element(Detector& theDetector, xml_h xml_ent, SensitiveDetector sens_det)
{
    xml_det_t x_det = xml_ent;
    std::string det_name = x_det.nameStr();
    sens_det.setType( "calorimeter" );

    xml_dim_t calo_dim = x_det.dimensions();
    const double ECaloX = calo_dim.x();
    const double ECaloY = calo_dim.y();
    const double ECaloZ = calo_dim.z();
    Box ECaloBox { ECaloX / 2., ECaloY / 2., ECaloZ / 2. };
    Volume ECaloVol { "ECalVol", ECaloBox, theDetector.material("Air") };

    DetElement subdet(det_name, x_det.id());
    Volume motherVolume = theDetector.pickMotherVolume(subdet);    
    PlacedVolume ECaloPlaced = motherVolume.placeVolume(ECaloVol, x_det.position());
    subdet.setPlacement(ECaloPlaced);
    return subdet;
}

DECLARE_DETELEMENT(NA64ECal_base, create_element)
