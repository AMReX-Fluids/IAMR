
#include <NavierStokes.H>
#include <NS_error_F.H>
#include <AMReX_ErrorList.H>
#include <AMReX_ParmParse.H>

using std::string;

using namespace amrex;

void
NavierStokes::error_setup()
{
    //
    // Dynamically generated error tagging functions
    //
    std::string amr_prefix = "amr";
    ParmParse ppamr(amr_prefix);
    Vector<std::string> refinement_indicators;
    ppamr.queryarr("refinement_indicators",refinement_indicators,0,ppamr.countval("refinement_indicators"));
    for (int i=0; i<refinement_indicators.size(); ++i)
    {
        std::string ref_prefix = amr_prefix + "." + refinement_indicators[i];

        ParmParse ppr(ref_prefix);
        RealBox realbox;
        if (ppr.countval("in_box_lo")) {
            std::vector<Real> box_lo(BL_SPACEDIM), box_hi(BL_SPACEDIM);
            ppr.getarr("in_box_lo",box_lo,0,box_lo.size());
            ppr.getarr("in_box_hi",box_hi,0,box_hi.size());
            realbox = RealBox(&(box_lo[0]),&(box_hi[0]));
        }

        AMRErrorTagInfo info;

        if (realbox.ok()) {
            info.SetRealBox(realbox);
        }
        if (ppr.countval("start_time") > 0) {
            Real min_time; ppr.get("start_time",min_time);
            info.SetMinTime(min_time);
        }
        if (ppr.countval("end_time") > 0) {
            Real max_time; ppr.get("end_time",max_time);
            info.SetMaxTime(max_time);
        }
        if (ppr.countval("max_level") > 0) {
            int max_level; ppr.get("max_level",max_level);
            info.SetMaxLevel(max_level);
        }

        if (ppr.countval("value_greater")) {
            Real value; ppr.get("value_greater",value);
            std::string field; ppr.get("field_name",field);
            errtags.push_back(AMRErrorTag(value,AMRErrorTag::GREATER,field,info));
        }
        else if (ppr.countval("value_less")) {
            Real value; ppr.get("value_less",value);
            std::string field; ppr.get("field_name",field);
            errtags.push_back(AMRErrorTag(value,AMRErrorTag::LESS,field,info));
        }
        else if (ppr.countval("vorticity_greater")) {
            Real value; ppr.get("vorticity_greater",value);
            const std::string field="mag_vort";
            errtags.push_back(AMRErrorTag(value,AMRErrorTag::VORT,field,info));
        }
        else if (ppr.countval("adjacent_difference_greater")) {
            Real value; ppr.get("adjacent_difference_greater",value);
            std::string field; ppr.get("field_name",field);
            errtags.push_back(AMRErrorTag(value,AMRErrorTag::GRAD,field,info));
        }
        else if (realbox.ok())
        {
            errtags.push_back(AMRErrorTag(info));
        }
	// //
	// // Could create a user defined function here as outlined below.
	// // However, this only allows you to use one "field", i.e. one
	// // component of State or a derived value (as defined in NS_setup.cpp).
	// // For all cases I can think of, a better option is to create a
	// // derived value in NS_setup.cpp and use one of the comparisons
	// // defined above (eg. value_greater will tag based on
	// // derived_value > value).
	// //
        // else if (ppr.countval("value")) {
        //     Real value; ppr.get("value",value);
        //     std::string field; ppr.get("field_name",field);

	//     // set ngrow for "field" based on what errFunc needs
	//     int ngrow = ;
	//     AMRErrorTag::UserFunc* errFunc;
	//     //
	//     // define error estimation function
	//     //
	    
	//     errtags.push_back(AMRErrorTag(errFunc,field,ngrow,info));
        // }
        else {
            Abort(std::string("Unrecognized refinement indicator for " + refinement_indicators[i]).c_str());
        }
    }

    //
    // User-defined error estimation functions 
    //
}
