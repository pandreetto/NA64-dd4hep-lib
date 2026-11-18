#include "XML/Utilities.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"

using namespace dd4hep;

static Ref_t create_element(Detector& theDetector, xml_h xml_ent, SensitiveDetector sens_det)
{
    xml_det_t x_det = xml_ent;
    std::string det_name = x_det.nameStr();
    sens_det.setType( "calorimeter" );

    const double ECaloWidth = theDetector.constant<double>("ECALWidth");
    const double ECaloHeight = theDetector.constant<double>("ECALHeight");
    const double ECaloDepth = theDetector.constant<double>("ECALDepth");

    const double ECaloX = theDetector.constant<double>("ECALPosX");
    const double ECaloY = theDetector.constant<double>("ECALPosY");
    const double ECaloZ = theDetector.constant<double>("ECALPosZ");

    const double FaceThickness = 1.5 * mm;
    const double ECaloStart = theDetector.constant<double>("ECALStart");
    const int numberOfLayers = theDetector.constant<int>("NLayersECAL");
    const int cellNumX = theDetector.constant<int>("ECALCellNumX");
    const int cellNumY = theDetector.constant<int>("ECALCellNumY");
    const double converterDepth = theDetector.constant<double>("ECALConverterDepth");
    const double holeRadius = theDetector.constant<double>("ECALHoleRadius");

    const double LayerDepth = ECaloDepth / numberOfLayers;
    const double CellSizeX = ECaloWidth / cellNumX;
    const double CellSizeY = ECaloHeight / cellNumY;

    const double counterDepth = LayerDepth - converterDepth;
    if (counterDepth <= 0.0)
    {
        // TODO missing error handling
    }
    const double holeDX = CellSizeX / 4;
    const double holeDY = CellSizeY / 4;

    DetElement subdet(det_name, x_det.id());
    Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, xml_ent, subdet);
    dd4hep::xml::setDetectorTypeFlag(xml_ent, subdet);

    Assembly detAss { "ECalDetAssembly" };
    detAss.setVisAttributes(theDetector, "ECALDefaultVis");

    Assembly cellAss { "ECalCellAssembly" };
    cellAss.setVisAttributes(theDetector, "ECALDefaultVis");

    /* *********************************************************************
     * Cover of ECAL fibers
     * ********************************************************************* */
    Box CoverBox { ECaloWidth / 2, ECaloHeight / 2, FaceThickness / 2 };
    Volume CoverVol { "ECALScreenZVol", CoverBox, theDetector.material("Al") };
    CoverVol.setVisAttributes(theDetector, "ECALCoverVis");
    auto coverPos = Position(0.0, 0.0, (ECaloDepth + FaceThickness) / 2);
    detAss.placeVolume(CoverVol, 0, coverPos);

    /* *********************************************************************
     * ECAL converter
     * ********************************************************************* */
    xml_det_t x_caloconv = x_det.child(_Unicode(caloConverter));
    auto conv_material = x_caloconv.attr<std::string>(_Unicode(material));

    Box convBox { CellSizeX / 2, CellSizeY / 2, converterDepth / 2 };
    Tube convHole { 0.0, holeRadius, converterDepth / 2 + 1.0 * mm };
    SubtractionSolid cSolid1 { convBox, convHole, Position(holeDX, holeDY, 0.0)};
    SubtractionSolid cSolid2 { cSolid1, convHole, Position(holeDX, -holeDY, 0.0)};
    SubtractionSolid cSolid3 { cSolid2, convHole, Position(-holeDX, holeDY, 0.0)};
    SubtractionSolid cSolid4 { cSolid3, convHole, Position(-holeDX, -holeDY, 0.0)};

    Volume convVol { "ConverterVol", cSolid4, theDetector.material(conv_material) };
    convVol.setVisAttributes(theDetector, x_caloconv.visStr());

    auto conv_pos = Position(0.0, 0.0, -0.5 * LayerDepth + 0.5 * converterDepth);
    cellAss.placeVolume(convVol, conv_pos);

    /* *********************************************************************
     * ECAL counter
     * ********************************************************************* */
    xml_det_t x_calosc = x_det.child(_Unicode(caloCounter));
    auto sc_material = x_calosc.attr<std::string>(_Unicode(material));

    Box scBox { CellSizeX / 2, CellSizeY / 2, counterDepth / 2 };
    Tube scHole { 0.0, holeRadius, counterDepth / 2 + 1.0 * mm };
    SubtractionSolid sSolid1 { scBox, scHole, Position(holeDX, holeDY, 0.0)};
    SubtractionSolid sSolid2 { sSolid1, scHole, Position(holeDX, -holeDY, 0.0)};
    SubtractionSolid sSolid3 { sSolid2, scHole, Position(-holeDX, holeDY, 0.0)};
    SubtractionSolid sSolid4 { sSolid3, scHole, Position(-holeDX, -holeDY, 0.0)};

    Volume scVol { "CounterVol", sSolid4, theDetector.material(sc_material) };
    scVol.setVisAttributes(theDetector, x_calosc.visStr());
    if (x_calosc.isSensitive()) scVol.setSensitiveDetector(sens_det);

    auto sc_pos = Position(0.0, 0.0, 0.5 * LayerDepth - 0.5 * counterDepth);
    cellAss.placeVolume(scVol, sc_pos);

    /* *********************************************************************
     * Iron rods
     * ********************************************************************* */
    xml_det_t x_rod = x_det.child(_Unicode(rod));
    auto rod_material = x_rod.attr<std::string>(_Unicode(material));

    Tube rodTube { 0.0, holeRadius, LayerDepth / 2 };
    Volume rodVol { "RodVol", rodTube, theDetector.material(rod_material) };
    rodVol.setVisAttributes(theDetector, x_rod.visStr());

    cellAss.placeVolume(rodVol, 1, Position(holeDX, holeDY, 0.0));
    cellAss.placeVolume(rodVol, 2, Position(holeDX, -holeDY, 0.0));
    cellAss.placeVolume(rodVol, 3, Position(-holeDX, holeDY, 0.0));
    cellAss.placeVolume(rodVol, 4, Position(-holeDX, -holeDY, 0.0));

    /* *********************************************************************
     * Layers
     * ********************************************************************* */
    Box layerBox { ECaloWidth / 2., ECaloHeight / 2., LayerDepth / 2 };
    Volume layerVol { "ECalLayerVol", layerBox, theDetector.material("Air") };
    layerVol.setVisAttributes(theDetector, "ECALDefaultVis");

    for (int k = 0; k < numberOfLayers; k++)
    {
        PlacedVolume layerPlaced = detAss.placeVolume(layerVol, k + 1,
            Position(0., 0., -ECaloDepth / 2. + LayerDepth / 2. + k * LayerDepth));
        layerPlaced.addPhysVolID("layerid", k);
    }

    for (int j = 0; j < cellNumX; j++)
    {
        for (int h = 0; h < cellNumY; h++)
        {
            auto cell_pos = Position(-ECaloWidth / 2. + CellSizeX / 2 + j * CellSizeX,
                -ECaloHeight / 2. + CellSizeY / 2 + h * CellSizeY, 0);
            int cell_idx = j * cellNumY + h;

            PlacedVolume cellPlaced = layerVol.placeVolume(cellAss, cell_idx, cell_pos);
            cellPlaced.addPhysVolID("cellx", j).addPhysVolID("celly", h);
        }
    }

    auto calo_pos = Position(ECaloX, ECaloY, ECaloZ);
    envelope.placeVolume(detAss, calo_pos);

    return subdet;
}

DECLARE_DETELEMENT(NA64ECal_base, create_element)
