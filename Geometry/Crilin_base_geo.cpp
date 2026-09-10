#include "XML/Utilities.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"

using namespace dd4hep;

static Ref_t create_element(Detector& theDetector, xml_h xml_ent, SensitiveDetector sens_det)
{
    xml_det_t x_det = xml_ent;
    std::string det_name = x_det.nameStr();
    sens_det.setType( "crilin" );

    const double PosX = theDetector.constant<double>("CrilinPosX");
    const double PosY = theDetector.constant<double>("CrilinPosY");
    const double PosZ = theDetector.constant<double>("CrilinPosZ");
    const double CrystalSize = theDetector.constant<double>("CrilinCrystalSize");
    const double CrystalLength = theDetector.constant<double>("CrilinCrystalLength");
    const double WrapThickness = theDetector.constant<double>("CrilinWrapThickness");
    const double HoneyCombThickness = theDetector.constant<double>("CrilinHoneyCombThickness");
    const double ElectroThickness = theDetector.constant<double>("CrilinElectroThickness");
    const double ShellThickness = theDetector.constant<double>("CrilinShellThickness");
    const double SiPMSize = theDetector.constant<double>("CrilinSiPMSize");
    const double SiPMThickness = theDetector.constant<double>("CrilinSiPMThickness");
    const int NumCellX = theDetector.constant<int>("CrilinNumCellX");
    const int NumCellY = theDetector.constant<int>("CrilinNumCellY");
    const int NumLayer = theDetector.constant<int>("CrilinNumLayer");

    /* *********************************************************************
     * Cell
     * ********************************************************************* */
	double WrapSize = CrystalSize + 2 * WrapThickness;
	double CellSize = WrapSize + 2 * HoneyCombThickness;

	xml_det_t x_crystal = x_det.child(_Unicode(crystal));
    auto crystal_material = x_crystal.attr<std::string>(_Unicode(material));

    Box crysBox { CrystalSize * 0.5, CrystalSize * 0.5, CrystalLength * 0.5 };
    Volume crysVol { "CrystalVol", crysBox, theDetector.material(crystal_material) };
    crysVol.setVisAttributes(theDetector, x_crystal.visStr());

    xml_det_t x_wrap = x_det.child(_Unicode(wrap));
    auto wrap_material = x_wrap.attr<std::string>(_Unicode(material));
	Box wrapBox { WrapSize * 0.5, WrapSize * 0.5, CrystalLength * 0.5 };
	SubtractionSolid wrapShell { wrapBox, crysBox, Position(0.0, 0.0, 0.0) };
	Volume wrapVol { "WrapVol", wrapShell, theDetector.material(wrap_material) };
	wrapVol.setVisAttributes(theDetector, x_wrap.visStr());

	xml_det_t x_hcomb = x_det.child(_Unicode(honeycomb));
	auto hcomb_material = x_hcomb.attr<std::string>(_Unicode(material));
	Box hcombBox { CellSize * 0.5, CellSize * 0.5, CrystalLength * 0.5 };
	SubtractionSolid hcombShell { hcombBox, wrapBox, Position(0.0, 0.0, 0.0) };
	Volume hcombVol { "hcombVol", hcombShell, theDetector.material(hcomb_material) };
	hcombVol.setVisAttributes(theDetector, x_hcomb.visStr());

    /* *********************************************************************
     * Electronics
     * ********************************************************************* */
	double electro_width = NumCellX * CellSize + 2 * ShellThickness;
	double electro_height = NumCellY * CellSize + 2 * ShellThickness;

	xml_det_t x_electro = x_det.child(_Unicode(electro));
	auto electro_material = x_electro.attr<std::string>(_Unicode(material));
	Box electroBox { electro_width * 0.5, electro_height * 0.5, ElectroThickness * 0.5 };
	Volume electroVol { "ElectroVol", electroBox, theDetector.material(electro_material) };
	electroVol.setVisAttributes(theDetector, x_electro.visStr());

    /* *********************************************************************
     * SiPM
     * ********************************************************************* */
	xml_det_t x_sipm = x_det.child(_Unicode(siPM));
	auto sipm_material = x_sipm.attr<std::string>(_Unicode(material));
	Box siPMBox { SiPMSize * 0.5, SiPMSize * 0.5, SiPMThickness * 0.5 };
	Volume siPMVol { "SiPMVol", siPMBox, theDetector.material(sipm_material) };
	siPMVol.setVisAttributes(theDetector, x_sipm.visStr());

    /* *********************************************************************
     * Shell (TODO missing holes)
     * ********************************************************************* */
	xml_det_t x_shell = x_det.child(_Unicode(shell));
	auto shell_material = x_shell.attr<std::string>(_Unicode(material));
	Box shBoxOut { electro_width * 0.5, electro_height * 0.5, CrystalLength * 0.5 };
	Box shBoxIn { NumCellX * CellSize * 0.5, NumCellY * CellSize * 0.5, CrystalLength * 0.5 };
	SubtractionSolid shShell { shBoxOut, shBoxIn, Position(0.0, 0.0, 0.0) };
	Volume shellVol { "ShellVol", shShell, theDetector.material(sipm_material) };
	shellVol.setVisAttributes(theDetector, x_shell.visStr());

	/* *********************************************************************
     * Sub Detector
     * ********************************************************************* */
    double layerSize = CrystalLength + SiPMThickness + ElectroThickness;
    double caloSizeX = electro_width;
    double caloSizeY = electro_height;      // TODO missing kapton tape
    double caloSizeZ = layerSize * NumLayer;
    Box caloBox { caloSizeX, caloSizeY, caloSizeZ };
    Volume caloVol { "CrilinVol", caloBox, theDetector.material("Air") };
    caloVol.setVisAttributes(theDetector, "CrilinDefaultVis");

    double x0 = (NumCellX - 1) * CellSize * -0.5;
    double y0 = (NumCellY - 1) * CellSize * -0.5;
    double z0 = caloSizeZ * -0.5;
    double siPMOffset = (CrystalSize - (3. * SiPMSize)) * 0.5;   // TODO check factor 3.0

    for (int k = 0; k < NumLayer; k++)
    {
    	double zOffset = k * layerSize + z0;

    	// Cells placement
        for (int j = 0; j < NumCellY; j++)
        {
        	for (int i = 0; i < NumCellX; i++)
        	{
        		int cellIdx = NumCellX * (NumCellY * k + j) + i;
        		auto cellCenter = Position(x0 + CellSize * i, y0 + CellSize * j,
        				                 zOffset + CrystalLength * 0.5);

        		PlacedVolume crystalPlaced = caloVol.placeVolume(crysVol, cellIdx, cellCenter);
        		crystalPlaced.addPhysVolID("cellid", cellIdx);
        		caloVol.placeVolume(wrapVol, cellIdx, cellCenter);
        		caloVol.placeVolume(hcombVol, cellIdx, cellCenter);

				auto sPos0 = Position(siPMOffset, siPMOffset, (CrystalLength + SiPMThickness) * 0.5);
				caloVol.placeVolume(siPMVol, 4 * cellIdx, cellCenter + sPos0);
				auto sPos1 = Position(siPMOffset, -1 * siPMOffset, (CrystalLength + SiPMThickness) * 0.5);
				caloVol.placeVolume(siPMVol, 4 * cellIdx + 1, cellCenter + sPos1);
				auto sPos2 = Position(-1 * siPMOffset, siPMOffset, (CrystalLength + SiPMThickness) * 0.5);
				caloVol.placeVolume(siPMVol, 4 * cellIdx + 2, cellCenter + sPos2);
				auto sPos3 = Position(-1 * siPMOffset, -1 * siPMOffset, (CrystalLength + SiPMThickness) * 0.5);
				caloVol.placeVolume(siPMVol, 4 * cellIdx + 3, cellCenter + sPos3);
        	}
        }
        // Electronics placement
		auto electroPos = Position(0, 0,
				zOffset + CrystalLength + SiPMThickness + ElectroThickness * 0.5);
        caloVol.placeVolume(electroVol, k, electroPos);

        // Shell placement
        auto shellPos = Position(0, 0, zOffset + CrystalLength * 0.5);
        caloVol.placeVolume(shellVol, k, shellPos);
    }

    DetElement subdet(det_name, x_det.id());
    Volume motherVolume = theDetector.pickMotherVolume(subdet);
    auto caloPos = Position(PosX, PosY, PosZ);
	PlacedVolume CrilinPlaced = motherVolume.placeVolume(caloVol, caloPos);
	subdet.setPlacement(CrilinPlaced);

    return subdet;
}

DECLARE_DETELEMENT(NA64Crilin_base, create_element)
