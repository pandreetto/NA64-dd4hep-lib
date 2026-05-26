#include "XML/Utilities.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"

using namespace dd4hep;

static Ref_t create_element(Detector& theDetector, xml_h xml_ent, SensitiveDetector sens_det)
{
    xml_det_t x_det = xml_ent;
    std::string det_name = x_det.nameStr();
    sens_det.setType( "lyso" );

    const double LYSOModWidth = theDetector.constant<double>("LYSOModWidth");
    const double LYSOModHeight = theDetector.constant<double>("LYSOModHeight");
    const double LYSOModDepth = theDetector.constant<double>("LYSOModDepth");

    const int cellNumX = theDetector.constant<int>("LYSOCellNumX");
    const int cellNumY = theDetector.constant<int>("LYSOCellNumY");
    const int cellNumZ = theDetector.constant<int>("LYSOCellNumZ");

    const double frontWidth = theDetector.constant<double>("LYSOFrontWidth");
    const double frontHeight = theDetector.constant<double>("LYSOFrontHeight");
    const double frontDepth = theDetector.constant<double>("LYSOFrontDepth");
    const double sideWidth = theDetector.constant<double>("LYSOSideWidth");
    const double sideDepth = theDetector.constant<double>("LYSOSideDepth");
    const double backDepth = theDetector.constant<double>("LYSOBackDepth");
    const int numberOfModules = theDetector.constant<int>("NModulesLYSO");
    const double modDistance = theDetector.constant<double>("LYSOModDist");

    const double frontCenterX = theDetector.constant<double>("LYSOFrontCenterX");
    const double frontCenterY = theDetector.constant<double>("LYSOFrontCenterY");
    const double frontCenterZ = theDetector.constant<double>("LYSOFrontCenterZ");

    Box LYSOModBox { LYSOModWidth / 2., LYSOModHeight / 2., LYSOModDepth / 2. };
    Volume LYSOModVol { "LYSOModVol", LYSOModBox, theDetector.material("Air") };
    LYSOModVol.setVisAttributes(theDetector, "LYSODefaultVis");

    /* *********************************************************************
     * Carbon fiber front box
     * ********************************************************************* */
    xml_det_t x_lysofront = x_det.child(_Unicode(lysoFront));
    auto lfront_material = x_lysofront.attr<std::string>(_Unicode(material));

    Box FrontBox { frontWidth / 2., frontHeight / 2., frontDepth / 2. };
    Volume FrontVol { "FrontVol", FrontBox, theDetector.material(lfront_material) };

    auto front_pos = Position(0, 0, -0.5 * sideDepth);
    PlacedVolume frontPlaced = LYSOModVol.placeVolume(FrontVol, front_pos);

    /* *********************************************************************
     * Aluminum sides and back box
     * ********************************************************************* */
    xml_det_t x_lysoshell = x_det.child(_Unicode(lysoShell));
    auto lshell_material = x_lysoshell.attr<std::string>(_Unicode(material));

    Box SideBox { sideWidth / 2., frontHeight / 2., sideDepth / 2. };
    Volume SideVol { "SideVol", SideBox, theDetector.material(lshell_material) };

    auto rside_pos = Position((0.5 * frontWidth + 0.5 * frontDepth), 0, 0);
    PlacedVolume rsidePlaced = LYSOModVol.placeVolume(SideVol, 0, rside_pos);
    auto lside_pos = Position((-0.5 * frontWidth - 0.5 * frontDepth), 0, 0);
    PlacedVolume lsidePlaced = LYSOModVol.placeVolume(SideVol, 1, lside_pos);


    Box BackBox { frontWidth / 2., frontHeight / 2., backDepth / 2. };
    Volume BackVol { "BackVol", BackBox, theDetector.material(lshell_material) };
    auto back_pos = Position(0, 0, (0.5 * sideDepth - 0.5 * backDepth));
    PlacedVolume backPlaced = LYSOModVol.placeVolume(BackVol, back_pos);

    /* *********************************************************************
     * LYSO matrix
     * ********************************************************************* */
    xml_det_t x_lysomatrix = x_det.child(_Unicode(lysoMatrix));
    auto matrix_material = x_lysomatrix.attr<std::string>(_Unicode(material));

    Box matrixBox { 32.0 * mm, 45.0 * mm, 48.0 * mm };
    Volume matrixVol { "MatrixVol", matrixBox, theDetector.material(matrix_material) };
    auto matrix_pos = Position(0.0 * mm, -42.0 * mm, -5.0 * mm);
    PlacedVolume matrixPlaced = LYSOModVol.placeVolume(matrixVol, matrix_pos);
    matrixPlaced.addPhysVolID("cellid", 1);

    if (x_lysomatrix.isSensitive()) matrixVol.setSensitiveDetector(sens_det);

    /* *********************************************************************
     * Sub Detector
     * ********************************************************************* */
    DetElement subdet(det_name, x_det.id());

    Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, xml_ent, subdet);
    dd4hep::xml::setDetectorTypeFlag(xml_ent, subdet);

    DetElement mdet0(subdet, "LYSO_module_0", x_det.id());
    for (int k = 0; k < numberOfModules; k++)
    {
        auto lyso_pos = Position(frontCenterX - k * (LYSOModWidth + modDistance),
                                 frontCenterY + 42 * mm,
                                 frontCenterZ + LYSOModDepth);
        PlacedVolume modPlaced = envelope.placeVolume(LYSOModVol, lyso_pos);
        modPlaced.addPhysVolID("moduleid", k);
        if (k == 0)
        {
            mdet0.setPlacement(modPlaced);
            continue;
        }

        DetElement mdetk = mdet0.clone("LYSO_module_" + k, x_det.id());
        mdetk.setPlacement(modPlaced);
        subdet.add(mdetk);
    }

    return subdet;
}

DECLARE_DETELEMENT(NA64LYSO_base, create_element)
