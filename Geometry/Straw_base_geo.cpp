#include "XML/Utilities.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"

#include <format>

using namespace dd4hep;

static Ref_t create_element(Detector& theDetector, xml_h xml_ent, SensitiveDetector sens_det)
{
    xml_det_t x_det = xml_ent;
    std::string det_name = x_det.nameStr();
    sens_det.setType( "Straw" );

    DetElement subdet(det_name, x_det.id());
    Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, xml_ent, subdet);
    dd4hep::xml::setDetectorTypeFlag(xml_ent, subdet);

    const double outerRadius = theDetector.constant<double>("StrawOuterRadius");
    const double shellThick = theDetector.constant<double>("ShellThickness");
    const double kapThick = theDetector.constant<double>("KaptonThickness");
    const double wireThick = theDetector.constant<double>("WireThickness");

    for (xml_coll_t statItem(xml_ent, _U(station)); statItem; statItem++)
    {
    	xml_comp_t station(statItem);
    	int stationId = station.id();
    	const double PosX = station.attr<double>(_Unicode(posx), 0.0);
    	const double PosY = station.attr<double>(_Unicode(posy), 0.0);
    	const double PosZ = station.attr<double>(_Unicode(posz), 0.0);
    	const double RotX = station.attr<double>(_Unicode(rotx), 0.0);
    	const double RotY = station.attr<double>(_Unicode(roty), 0.0);
    	const double RotZ = station.attr<double>(_Unicode(rotz), 0.0);
    	const double sLength = station.attr<double>(_Unicode(length));
    	const int nWires = station.attr<int>(_Unicode(wires));

    	/* *********************************************************************
		 * Straw tube
		 * ********************************************************************* */
		auto r1 = outerRadius - shellThick;
		auto r2 = r1 - kapThick;

		Tube shTube { r1, outerRadius, sLength * 0.5 };
		Volume shTVol { "AlVol", shTube, theDetector.material("Aluminum") };
		shTVol.setVisAttributes(theDetector, "StrawShellVis");

		Tube kapTube { r2 , r1, sLength * 0.5 };
		Volume kapTVol { "KapVol", kapTube, theDetector.material("Kapton") };
		kapTVol.setVisAttributes(theDetector, "StrawKaptonVis");

		Tube gasTube { wireThick, r2, sLength * 0.5 };
		Volume gasTVol { "GasVol", gasTube, theDetector.material("StrawGas") };
		gasTVol.setVisAttributes(theDetector, "StrawGasVis");

		Tube wireTube { 0, wireThick, sLength * 0.5 };
		Volume wireTVol { "WireVol", wireTube, theDetector.material("Tungsten") };
		wireTVol.setVisAttributes(theDetector, "StrawWireVis");

    	/* *********************************************************************
		 * Straw assembly
		 * ********************************************************************* */
		Assembly stAssembly { std::format("StrawAssembly_{}", stationId) };
		stAssembly.setVisAttributes(theDetector, "StrawDefaultVis");

		for (int idx = 0; idx < nWires; idx++)
		{
			auto t1_center = Position(outerRadius, (2 * idx - nWires + 1)  * outerRadius, 0.0);
			auto t2_center = Position(-1 * outerRadius, (2 * idx - nWires)  * outerRadius, 0.0);

			stAssembly.placeVolume(shTVol, t1_center);
			stAssembly.placeVolume(kapTVol, t1_center);
			auto gas1PV = stAssembly.placeVolume(gasTVol, t1_center);
			stAssembly.placeVolume(wireTVol, t1_center);

			stAssembly.placeVolume(shTVol, t2_center);
			stAssembly.placeVolume(kapTVol, t2_center);
			auto gas2PV = stAssembly.placeVolume(gasTVol, t2_center);
			stAssembly.placeVolume(wireTVol, t2_center);

			gas1PV.addPhysVolID("layerid", 0).addPhysVolID("tubeid", idx);
			gas2PV.addPhysVolID("layerid", 1).addPhysVolID("tubeid", idx);
		}
		/* *********************************************************************
		 * Straw station
		 * ********************************************************************* */
		auto st_pos = Position(PosX, PosY, PosZ);
		auto st_rot = RotationZYX(RotZ, RotY, RotX);
		auto st_transf = Transform3D(st_rot, st_pos);
		PlacedVolume modPlaced = envelope.placeVolume(stAssembly, st_transf);
		modPlaced.addPhysVolID("stationid", stationId);
    }
    return subdet;
}

DECLARE_DETELEMENT(NA64Straw_base, create_element)
