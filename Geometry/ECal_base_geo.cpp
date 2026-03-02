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
    const double ECaloWidth = calo_dim.x();
    const double ECaloHeight = calo_dim.y();
    const double ECaloDepth = calo_dim.z();

    auto calo_pos = x_det.position();
    const double ECaloX = calo_pos.x();
    const double ECaloY = calo_pos.y();
    const double ECaloZ = calo_pos.z();

    const double FaceThickness = 1.5 * mm;
    const double ECaloStart = theDetector.constant<double>("ECALStart");
    const int numberOfLayers = theDetector.constant<int>("NLayersECAL");
    const int cellNumX = theDetector.constant<int>("CellNumX");
    const int cellNumY = theDetector.constant<int>("CellNumY");
    const double converterDepth = theDetector.constant<double>("ConverterDepth");

    const double LayerDepth = ECaloDepth / numberOfLayers;
    const double CellSizeX = ECaloWidth / cellNumX;
    const double CellSizeY = ECaloHeight / cellNumY;

    const double counterDepth = LayerDepth - converterDepth;
    if (counterDepth <= 0.0)
    {
        // TODO missing error handling
    }

    Box ECaloBox { ECaloWidth / 2., ECaloHeight / 2., ECaloDepth / 2. };
    Volume ECaloVol { "ECalVol", ECaloBox, theDetector.material("Air") };
    ECaloVol.setVisAttributes(theDetector, "DefaultVis");

    DetElement subdet(det_name, x_det.id());
    Volume motherVolume = theDetector.pickMotherVolume(subdet);    
    PlacedVolume ECaloPlaced = motherVolume.placeVolume(ECaloVol, calo_pos);
    subdet.setPlacement(ECaloPlaced);

//    Box FaceBox { ECaloWidth / 2, ECaloHeight / 2, FaceThickness / 2 };
//    Volume FaceVol { "ECALScreenZVol", FaceBox, theDetector.material("Al") };
//    motherVolume.placeVolume(FaceVol, Position(ECaloX, ECaloY, ECaloZ + FaceThickness));

    Box layerBox { ECaloWidth / 2., ECaloHeight / 2., LayerDepth / 2 };
    Volume layerVol { "ECalLayerVol", layerBox, theDetector.material("Air") };
    layerVol.setVisAttributes(theDetector, "DefaultVis");

    Box cellBox { CellSizeX / 2, CellSizeY / 2, LayerDepth / 2 };
    Volume cellVol { "ECalCellVol", cellBox, theDetector.material("Air") };
    cellVol.setVisAttributes(theDetector, "DefaultVis");

    /* *********************************************************************
     * ECAL converter
     * ********************************************************************* */
    xml_det_t x_caloconv = x_det.child(_Unicode(caloConverter));
    auto conv_material = x_caloconv.attr<std::string>(_Unicode(material));

    Box convBox { CellSizeX / 2, CellSizeY / 2, converterDepth / 2 };
    Volume convVol { "ConverterVol", convBox, theDetector.material(conv_material) };
    convVol.setVisAttributes(theDetector, x_caloconv.visStr());

    auto conv_pos = Position(0.0, 0.0, -0.5 * LayerDepth + 0.5 * converterDepth);
    cellVol.placeVolume(convVol, 1, conv_pos);

    /* *********************************************************************
     * ECAL counter
     * ********************************************************************* */
    xml_det_t x_calosc = x_det.child(_Unicode(caloCounter));
    auto sc_material = x_calosc.attr<std::string>(_Unicode(material));

    Box scBox { CellSizeX / 2, CellSizeY / 2, converterDepth / 2 };
    Volume scVol { "CounterVol", scBox, theDetector.material(sc_material) };
    convVol.setVisAttributes(theDetector, x_calosc.visStr());

    auto sc_pos = Position(0.0, 0.0, 0.5 * LayerDepth - 0.5 * counterDepth);
    cellVol.placeVolume(scVol, 1, sc_pos);

    /* *********************************************************************
     * Layers
     * ********************************************************************* */
    for (int k = 0; k < numberOfLayers; k++)
    {
        PlacedVolume layerPlaced = ECaloVol.placeVolume(layerVol, k,
            Position(0., 0., -ECaloDepth / 2. + LayerDepth / 2. + k * LayerDepth));
        layerPlaced.addPhysVolID("layer", k);

        for (int j = 0; j < cellNumX; j++)
        {
            for (int h = 0; h < cellNumY; h++)
            {
                auto cell_pos = Position(-ECaloWidth / 2. + CellSizeX / 2 + j * CellSizeX,
                    -ECaloHeight / 2. + CellSizeY / 2 + h * CellSizeY, 0);
                int cell_idx = (k * cellNumX + j) * cellNumY + h;

                PlacedVolume cellPlaced = layerVol.placeVolume(cellVol, cell_idx, cell_pos);
                cellPlaced.addPhysVolID("cell", cell_idx);
            }
        }
    }

    return subdet;
}

DECLARE_DETELEMENT(NA64ECal_base, create_element)
