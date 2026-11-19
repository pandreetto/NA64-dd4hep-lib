#include "XML/Utilities.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"

#include <format>
#include <vector>

using namespace dd4hep;
using std::vector;

static Ref_t create_element(Detector& theDetector, xml_h xml_ent, SensitiveDetector sens_det)
{
    xml_det_t x_det = xml_ent;
    std::string det_name = x_det.nameStr();
    int det_id = x_det.id();
    sens_det.setType( "calorimeter" );

    const double HCaloWidth = theDetector.constant<double>("HCALWidth");
    const double HCaloHeight = theDetector.constant<double>("HCALHeight");
    const double HCaloDepth = theDetector.constant<double>("HCALDepth");

    const double HCaloModWidth = theDetector.constant<double>("HCALModWidth");
    const double HCaloModHeight = theDetector.constant<double>("HCALModHeight");

    const int numberOfLayers = theDetector.constant<int>("NLayersHCAL");
    const double converterDepth = theDetector.constant<double>("HCALConverterDepth");
    const double counterDepth = theDetector.constant<double>("HCALCounterDepth");
    const int cellNumX = theDetector.constant<int>("HCALCellNumX");
    const int cellNumY = theDetector.constant<int>("HCALCellNumY");

    const double LayerDepth = HCaloDepth / numberOfLayers;
    const double converterShellDepth = LayerDepth - (converterDepth + counterDepth);
    const double CellSizeX = HCaloWidth / cellNumX;
    const double CellSizeY = HCaloHeight / cellNumY;

    const double HCaloStart = theDetector.constant<double>("HCALStart");
    const double HCaloModDepth = HCaloDepth + converterDepth;

    Box HCaloBox { HCaloWidth / 2., HCaloHeight / 2., HCaloDepth / 2. };
    Volume HCaloVol { "HCalVol", HCaloBox, theDetector.material("Air") };
    HCaloVol.setVisAttributes(theDetector, "HCALDefaultVis");

    /* *********************************************************************
     * Layers
     * ********************************************************************* */
    Box layerBox { HCaloWidth / 2., HCaloHeight / 2., LayerDepth / 2 };
    Volume layerVol { "HCalLayerVol", layerBox, theDetector.material("Air") };
    layerVol.setVisAttributes(theDetector, "HCALDefaultVis");

    for (int k = 0; k < numberOfLayers; k++)
    {
        PlacedVolume layerPlaced = HCaloVol.placeVolume(layerVol, k,
            Position(0., 0., -HCaloDepth / 2. + LayerDepth / 2. + k * LayerDepth));
        layerPlaced.addPhysVolID("layerid", k);
    }

    /* *********************************************************************
     * HCAL converter
     * ********************************************************************* */
    xml_det_t x_caloconv = x_det.child(_Unicode(caloConverter));
    auto conv_material = x_caloconv.attr<std::string>(_Unicode(material));

    Box convBox { HCaloWidth / 2, HCaloHeight / 2, converterDepth / 2 };
    Volume convVol { "ConverterVol", convBox, theDetector.material(conv_material) };
    convVol.setVisAttributes(theDetector, x_caloconv.visStr());

    auto conv_pos = Position(0.0, 0.0, -0.5 * LayerDepth + 0.5 * converterDepth);
    PlacedVolume convPlaced = layerVol.placeVolume(convVol, 1, conv_pos);

    /* *********************************************************************
     * HCAL converter shell
     * ********************************************************************* */
    xml_det_t x_caloconvsh = x_det.child(_Unicode(caloConverterShell));
    auto convsh_material = x_caloconvsh.attr<std::string>(_Unicode(material));

    Box convShBox { HCaloWidth / 2, HCaloHeight / 2, converterShellDepth / 2 };
    Volume convShVol { "ConverterShellVol", convShBox, theDetector.material(convsh_material) };
    convShVol.setVisAttributes(theDetector, x_caloconvsh.visStr());

    double z_convsh = -0.5 * LayerDepth + converterDepth + counterDepth + 0.5 * converterShellDepth;
    auto convSh_pos = Position(0.0, 0.0, z_convsh);
    PlacedVolume convShPlaced = layerVol.placeVolume(convShVol, 1, convSh_pos);

    /* *********************************************************************
     * HCAL counter
     * ********************************************************************* */
    xml_det_t x_calosc = x_det.child(_Unicode(caloCounter));
    auto sc_material = x_calosc.attr<std::string>(_Unicode(material));

    Box scBox { HCaloWidth / 2, HCaloHeight / 2,  counterDepth / 2 };
    Volume scVol { "CounterVol", scBox, theDetector.material("Air") };
    scVol.setVisAttributes(theDetector, "HCALDefaultVis");

    auto sc_pos = Position(0.0, 0.0, -0.5 * LayerDepth + converterDepth + 0.5 * counterDepth);
    PlacedVolume scPlaced = layerVol.placeVolume(scVol, 1, sc_pos);

    Box cellBox { CellSizeX / 2, CellSizeY / 2, counterDepth };
    Volume cellVol { "CellVol", cellBox, theDetector.material(sc_material) };
    cellVol.setVisAttributes(theDetector, x_calosc.visStr());
    if (x_calosc.isSensitive()) cellVol.setSensitiveDetector(sens_det);

    for (int j = 0; j < cellNumX; j++)
    {
        for (int h = 0; h < cellNumY; h++)
        {
            auto cell_pos = Position(-HCaloWidth / 2. + CellSizeX / 2 + j * CellSizeX,
                -HCaloHeight / 2. + CellSizeY / 2 + h * CellSizeY, 0);
            int cell_idx = j * cellNumY + h;

            PlacedVolume cellPlaced = scVol.placeVolume(cellVol, cell_idx, cell_pos);
            cellPlaced.addPhysVolID("cellid", cell_idx);
        }
    }

    /* *********************************************************************
     * HCAL Modules
     * ********************************************************************* */
    // lateral shields
    xml_det_t x_calosh = x_det.child(_Unicode(caloShield));
    auto sh_material = x_calosh.attr<std::string>(_Unicode(material));

    Box shLatBox { (HCaloModWidth - HCaloWidth) / 4, HCaloModHeight / 2, HCaloModDepth / 2 };
    Volume shLatVol { "HCALShellX", shLatBox, theDetector.material(sh_material) };
    shLatVol.setVisAttributes(theDetector, x_calosh.visStr());
    
    // top/bottom shields
    Box shTBBox { HCaloWidth / 2, (HCaloModHeight - HCaloHeight) / 4, HCaloModDepth / 2 };
    Volume shTBVol { "HCALShellY", shTBBox, theDetector.material(sh_material) };
    shTBVol.setVisAttributes(theDetector, x_calosh.visStr());
    
    DetElement subdet(det_name, x_det.id());

    Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, xml_ent, subdet);
    dd4hep::xml::setDetectorTypeFlag(xml_ent, subdet);

    for (xml_coll_t modItem(xml_ent, _U(module)); modItem; modItem++)
    {
    	xml_comp_t hModule(modItem);
		int module_id = hModule.id();
    	const double PosX = hModule.attr<double>(_Unicode(posx), 0.0);
    	const double PosY = hModule.attr<double>(_Unicode(posy), 0.0);
    	const double PosZ = hModule.attr<double>(_Unicode(posz), 0.0);

    	Assembly modAssembly { std::format("HCalModAssembly_{}", module_id) };
    	modAssembly.setVisAttributes(theDetector, "HCALDefaultVis");

        auto mr_pos = Position(0.0, 0.0, -0.5 * HCaloModDepth + 0.5 * HCaloDepth);
        PlacedVolume cModPlaced = modAssembly.placeVolume(HCaloVol, mr_pos);
        auto bck_pos = Position(0.0, 0.0, 0.5 * HCaloModDepth - 0.5 * converterDepth);
        PlacedVolume bckPlaced = modAssembly.placeVolume(convVol, bck_pos);

        auto sh1_pos = Position(-0.25 * (HCaloModWidth + HCaloWidth), 0.0, 0.0);
        PlacedVolume sh1Placed = modAssembly.placeVolume(shLatVol, 0, sh1_pos);
        auto sh2_pos = Position(0.25 * (HCaloModWidth + HCaloWidth), 0.0, 0.0);
        PlacedVolume sh2Placed = modAssembly.placeVolume(shLatVol, 1, sh2_pos);

        auto sh3_pos = Position(0.0, -0.25 * (HCaloModHeight + HCaloHeight), 0.0);
        PlacedVolume sh3Placed = modAssembly.placeVolume(shTBVol, 2, sh3_pos);
        auto sh4_pos = Position(0.0, 0.25 * (HCaloModHeight + HCaloHeight), 0.0);
        PlacedVolume sh4Placed = modAssembly.placeVolume(shTBVol, 3, sh4_pos);

    	auto mod_pos = Position(PosX, PosY, PosZ);
    	PlacedVolume modPlaced = envelope.placeVolume(modAssembly, mod_pos);
		modPlaced.addPhysVolID("system", det_id).addPhysVolID("moduleid", module_id);
    }

    return subdet;
}

DECLARE_DETELEMENT(NA64HCal_base, create_element)
