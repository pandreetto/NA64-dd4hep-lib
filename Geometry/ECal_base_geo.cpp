#include "XML/Utilities.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"

using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::Ref_t;
using dd4hep::SensitiveDetector;


class NA64EcalBuild : public dd4hep::xml::tools::VolumeBuilder
{
public:
    NA64EcalBuild(dd4hep::Detector& theDetector, xml_elt_t xml_ent, dd4hep::SensitiveDetector sens_det);
    virtual ~NA64EcalBuild() {}
    dd4hep::Volume build_ecal();
};

NA64EcalBuild::NA64EcalBuild(dd4hep::Detector& theDetector, xml_elt_t xml_ent, dd4hep::SensitiveDetector sens_det) :
    dd4hep::xml::tools::VolumeBuilder(theDetector, xml_ent, sens_det)
{}

dd4hep::Volume NA64EcalBuild::build_ecal()
{
    double ecal_width = dd4hep::_toDouble("ECALWidth");
    double ecal_height = dd4hep::_toDouble("ECALHeight");
    double ecal_layer_depth = dd4hep::_toDouble("ECALLayerDepth");
    int ecal_n_layers = dd4hep::_toInt("NLayersECAL");

    dd4hep::Box ecal_box { ecal_width / 2, ecal_height / 2, ecal_layer_depth * ecal_n_layers / 2 };
    dd4hep::Volume ecal_vol { "ECalVol", ecal_box, description.material("Air") };
    return ecal_vol;
}

static Ref_t create_element(Detector& theDetector, xml_h xml_ent, SensitiveDetector sens_det)
{
    xml_det_t x_det = xml_ent;
    std::string name = x_det.nameStr();
    sens_det.setType( "calorimeter" );    

    NA64EcalBuild builder { theDetector, x_det, sens_det };

    dd4hep::Volume lvEcalVol = builder.build_ecal();
    dd4hep::PlacedVolume pv = builder.placeDetector(lvEcalVol);
    pv.addPhysVolID("system", x_det.id());

    return builder.detector;
}

DECLARE_DETELEMENT(NA64ECal_base, create_element)
