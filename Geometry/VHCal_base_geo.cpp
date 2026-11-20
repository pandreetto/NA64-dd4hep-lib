#include "XML/Utilities.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"

#include <format>

using namespace dd4hep;

static Ref_t create_element(Detector& theDetector, xml_h xml_ent, SensitiveDetector sens_det)
{
    xml_det_t x_det = xml_ent;
    std::string det_name = x_det.nameStr();
    int det_id = x_det.id();
    sens_det.setType("calorimeter");

    const double VHCaloWidth = theDetector.constant<double>("VHCALWidth");
    const double ConvDepth = theDetector.constant<double>("VHCALConvDepth");
    const double CountDepth = theDetector.constant<double>("VHCALCountDepth");
    const double GapDepth = theDetector.constant<double>("VHCALGapDepth");
    const double VHCaloHoleWidth = theDetector.constant<double>("VHCALHoleWidth");
    const double VHCaloHoleHeight = theDetector.constant<double>("VHCALHoleHeight");
    const int numberOfLayers = theDetector.constant<int>("VHCALLayer");

    const double LayerDepth = ConvDepth + CountDepth + GapDepth;
    const double VHCaloDepth = numberOfLayers * LayerDepth;

    /* *********************************************************************
     * Converter
     * ********************************************************************* */
    xml_det_t x_caloconv = x_det.child(_Unicode(caloConverter));
	auto conv_material = x_caloconv.attr<std::string>(_Unicode(material));

	Box convBox { 0.5 * VHCaloWidth, 0.5 * VHCaloWidth, 0.5 * ConvDepth };
	Box convHole { 0.5 * VHCaloHoleWidth, 0.5 * VHCaloHoleHeight, 0.5 * ConvDepth + 0.1 * mm };
	SubtractionSolid convDiff { convBox, convHole, Position(0.0, 0.0, 0.0)};
    Volume convVol { "ConvVol", convDiff, theDetector.material(conv_material) };
    convVol.setVisAttributes(theDetector, x_caloconv.visStr());

    /* *********************************************************************
     * Counter
     * ********************************************************************* */
    xml_det_t x_calosc = x_det.child(_Unicode(caloCounter));
	auto sc_material = x_calosc.attr<std::string>(_Unicode(material));

	Box scBox { 0.5 * VHCaloWidth, 0.5 * VHCaloWidth, 0.5 * CountDepth };
	Box scHole { 0.5 * VHCaloHoleWidth, 0.5 * VHCaloHoleHeight, 0.5 * CountDepth + 0.1 * mm };
	SubtractionSolid scDiff { scBox, scHole, Position(0.0, 0.0, 0.0)};
    Volume scVol { "ScVol", scDiff, theDetector.material(sc_material) };
    scVol.setVisAttributes(theDetector, x_calosc.visStr());
    if (x_calosc.isSensitive()) scVol.setSensitiveDetector(sens_det);

    /* *********************************************************************
     * Layer
     * ********************************************************************* */
    Box layerBox { 0.5 * VHCaloWidth, 0.5 * VHCaloWidth, 0.5 * LayerDepth };
    Volume layerVol { "LayerVol", layerBox, theDetector.material("Air")};
    layerVol.setVisAttributes(theDetector, "VHCALDefaultVis");

    auto conv_pos = Position(0.0, 0.0, 0.5 * (ConvDepth - LayerDepth));
    layerVol.placeVolume(convVol, conv_pos);
    auto sc_pos = Position(0.0, 0.0, 0.5 * (LayerDepth - CountDepth));
    layerVol.placeVolume(scVol, sc_pos);

    /* *********************************************************************
     * VHCAL Modules
     * ********************************************************************* */
    DetElement subdet(det_name, det_id);
    Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, xml_ent, subdet);
	dd4hep::xml::setDetectorTypeFlag(xml_ent, subdet);

	for (xml_coll_t modItem(xml_ent, _U(module)); modItem; modItem++)
	{
    	xml_comp_t hModule(modItem);
		int module_id = hModule.id();
    	const double PosX = hModule.attr<double>(_Unicode(posx), 0.0);
    	const double PosY = hModule.attr<double>(_Unicode(posy), 0.0);
    	const double PosZ = hModule.attr<double>(_Unicode(posz), 0.0);
    	const double RotY = hModule.attr<double>(_Unicode(roty), 0.0);

    	Assembly modAssembly { std::format("VHCalModAssembly_{}", module_id) };
    	modAssembly.setVisAttributes(theDetector, "VHCALDefaultVis");

    	auto mc_pos = Position(0.0, 0.0, 0.0);
        for (int k = 0; k < numberOfLayers; k++)
        {
            auto dis = Position(0., 0., -0.5 * VHCaloDepth + LayerDepth * 0.5 + k * LayerDepth);
        	PlacedVolume layerPlaced = modAssembly.placeVolume(layerVol, k, mc_pos + dis);
            layerPlaced.addPhysVolID("layerid", k);
        }

        auto mod_pos = Position(PosX, PosY, PosZ);
        auto mod_rot = RotationZYX(0.0, RotY, 0.0);
        auto mod_transf = Transform3D(mod_rot, mod_pos);
    	PlacedVolume modPlaced = envelope.placeVolume(modAssembly, mod_transf);
		modPlaced.addPhysVolID("system", det_id).addPhysVolID("moduleid", module_id);
	}
	return subdet;
}

DECLARE_DETELEMENT(NA64VHCal_base, create_element)
